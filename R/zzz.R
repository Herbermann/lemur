.onLoad <- function(libname, pkgname) {
  tryCatch({
    offsets <- find_zcorr_offsets()
    init_harmony_offsets(offsets$nrows_offset, offsets$ptr_offset)
    verify_harmony_zcorr_offset()
  }, error = function(e) {
    packageStartupMessage(
      "lemur: Could not initialize harmony Z_corr setter: ", conditionMessage(e), "\n",
      "align_harmony() will not work on this system.\n",
      "Please report this at https://github.com/const-ae/lemur/issues"
    )
  })
}