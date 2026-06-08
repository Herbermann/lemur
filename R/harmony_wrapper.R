find_zcorr_offsets <- function() {
  so_path <- system.file("libs", package = "harmony")
  so_file <- list.files(so_path, 
                        pattern = "harmony\\.so$|harmony\\.dll$|harmony\\.dylib$",
                        full.names = TRUE)[1]
  
  if (is.na(so_file)) stop("Cannot find harmony shared library")
  
  sysname <- Sys.info()["sysname"]
  
  if (sysname %in% c("Darwin", "Linux")) {
    asm <- system(
      paste("objdump -d", 
            if (sysname == "Darwin") "--macho" else "",
            shQuote(so_file),
            "| grep -A 80 '_ZN7harmony8getZcorrEv'"),
      intern = TRUE
    )
    
    # n_rows offset: first ldp or mov loading from 'this' register
    nrows_line <- grep("ldp.*\\[x[0-9]+,.*#0x|mov.*0x.*\\(%r", asm, value = TRUE)[1]
    if (is.na(nrows_line)) stop("Could not find n_rows load in getZcorr disassembly")
    nrows_offset <- strtoi(
      sub(".*#(0x[0-9a-f]+)\\].*|.*0x([0-9a-f]+)\\(%r.*", "\\1\\2", nrows_line), 
      16L
    )
    
    # data ptr: last load into argument register before memcpy call
    memcpy_idx <- grep("memcpy", asm)[1]
    if (is.na(memcpy_idx)) stop("Could not find memcpy in getZcorr disassembly")
    
    if (sysname == "Darwin") {
      # ARM64: ldr x1, [x19, #0xNNN]
      ptr_lines <- grep("ldr.*x1.*\\[x[0-9]+,.*#0x", asm[seq_len(memcpy_idx)])
      ptr_line  <- asm[tail(ptr_lines, 1)]
      ptr_offset <- strtoi(sub(".*#(0x[0-9a-f]+)\\].*", "\\1", ptr_line), 16L)
    } else {
      # x86_64: mov 0xNNN(%rXX), %rsi
      ptr_lines <- grep("mov.*0x.*\\(%r.*%rsi", asm[seq_len(memcpy_idx)])
      ptr_line  <- asm[tail(ptr_lines, 1)]
      ptr_offset <- strtoi(sub(".*0x([0-9a-f]+)\\(%r.*", "\\1", ptr_line), 16L)
    }
    
  } else {
    stop("Unsupported platform: ", sysname, 
         ". Please report this at https://github.com/const-ae/lemur/issues")
  }
  
  if (is.na(nrows_offset) || is.na(ptr_offset)) {
    stop("Could not parse offsets from disassembly")
  }
  
  list(nrows_offset = as.integer(nrows_offset), 
       ptr_offset   = as.integer(ptr_offset))
}

verify_harmony_zcorr_offset <- function() {
  ho <- harmony::RunHarmony(
    harmony::cell_lines$scaled_pcs,
    harmony::cell_lines$meta_data,
    "dataset", return_object = TRUE, verbose = FALSE
  )
  Z    <- ho$getZcorr()
  test <- matrix(99.0, nrow = nrow(Z), ncol = ncol(Z))
  
  ok <- tryCatch({
    set_harmony_Zcorr(ho$.pointer, test)
    ho$getZcorr()[1,1] == 99.0
  }, error = function(e) FALSE)
  
  if (!ok) {
    stop(
      "set_harmony_Zcorr() verification failed with harmony ", 
      packageVersion("harmony"), " on ", Sys.info()["sysname"], ".\n",
      "Please report this at https://github.com/const-ae/lemur/issues"
    )
  }
  invisible(TRUE)
}

# Modified harmony_init that also returns phi and setup args
harmony_init <- function(embedding, design_matrix,
                         theta = 2, lambda = 1, sigma = 0.1, 
                         nclust = min(round(ncol(embedding) / 30), 100),
                         tau = 0, block.size = 0.05, max.iter.cluster = 200,
                         epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                         verbose = TRUE) {

  mm_groups <- get_groups(design_matrix)
  n_groups <- length(unique(mm_groups))
  phi <- matrix(0, nrow = n_groups, ncol = ncol(embedding))
  phi[mm_groups + n_groups * (seq_along(mm_groups)-1)] <- 1
  phi <- as(phi, "sparseMatrix")

  N <- ncol(embedding)
  N_b <- MatrixGenerics::rowSums2(phi)
  B_vec <- rep(1, n_groups)

  theta <- rep_len(theta, n_groups)
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))
  lambda <- rep_len(lambda, n_groups)
  lambda_vec <- c(0, lambda)
  sigma <- rep_len(sigma, nclust)
  alpha <- 0.2
  batch.prop.cutoff <- 1e-5

  harmonyObj <- harmony::RunHarmony(embedding, mm_groups, nclust = nclust, 
                                     max_iter = 0, return_object = TRUE, 
                                     verbose = FALSE)
  harmonyObj$setup(
    embedding, phi, sigma, theta, lambda_vec, alpha, max.iter.cluster, 
    epsilon.cluster, epsilon.harmony, nclust, block.size, B_vec, 
    batch.prop.cutoff, verbose
  )
  
  harmony_init_clustering(harmonyObj)
  
  # Return object AND the setup args needed to reinitialize embedding later
  list(
    obj = harmonyObj,
    phi = phi,
    sigma = sigma,
    theta = theta,
    lambda_vec = lambda_vec,
    alpha = alpha,
    max.iter.cluster = max.iter.cluster,
    epsilon.cluster = epsilon.cluster,
    epsilon.harmony = epsilon.harmony,
    nclust = nclust,
    block.size = block.size,
    B_vec = B_vec,
    batch.prop.cutoff = batch.prop.cutoff
  )
}

#' Create an arbitrary Harmony object so that I can modify it later
#'
#' @returns The full [`harmony`] object (R6 reference class type).
#'
#' @keywords internal
harmony_new_object <- function(){
  Y <- randn(3, 100)
  harmony::RunHarmony(Y, rep(c("a", "b"), length.out = 100), nclust = 2, max.iter = 0, return_object = TRUE, verbose = FALSE)
  harmonyObj
}

harmony_init_clustering <- function(harmonyObj, iter.max = 25, nstart = 10){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  Z_corr <- harmonyObj$getZcorr()
  Z_cos <- sweep(Z_corr, 2, sqrt(colSums(Z_corr^2)), "/")
  harmonyObj$Y <- t(stats::kmeans(t(Z_cos), centers = harmonyObj$K, 
                                   iter.max = iter.max, nstart = nstart)$centers)
  harmonyObj$init_cluster_cpp()
  harmonyObj
}

harmony_max_div_clustering <- function(harmonyObj){
  stopifnot(is(harmonyObj, "Rcpp_harmony"))
  err_status <- harmonyObj$cluster_cpp()
  if (err_status == -1) {
    stop('terminated by user')
  } else if (err_status != 0) {
    stop(gettextf('Harmony exited with non-zero exit status: %d',
                  err_status))
  }
  harmonyObj
}

