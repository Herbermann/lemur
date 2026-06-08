// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cstdint>
#include <cstring>

static size_t g_nrows_offset = 0;
static size_t g_ptr_offset   = 0;

// [[Rcpp::export]]
void init_harmony_offsets(int nrows_offset, int ptr_offset) {
    g_nrows_offset = (size_t) nrows_offset;
    g_ptr_offset   = (size_t) ptr_offset;
}

// [[Rcpp::export]]
int sizeof_uword() { return (int)sizeof(arma::uword); }

// [[Rcpp::export]]
void set_harmony_Zcorr(SEXP harmonyPtr, const arma::mat& Z_new) {
    if (g_nrows_offset == 0) 
        Rcpp::stop("Harmony offsets not initialized.");
    
    uint8_t* base = (uint8_t*) R_ExternalPtrAddr(harmonyPtr);

    // Read as uint64 — harmony is always compiled 64-bit
    uint64_t nr, nc;
    memcpy(&nr, base + g_nrows_offset,     sizeof(uint64_t));
    memcpy(&nc, base + g_nrows_offset + 8, sizeof(uint64_t));

    if (nr != (uint64_t)Z_new.n_rows || nc != (uint64_t)Z_new.n_cols) {
        Rcpp::stop("Dimension mismatch: Z_corr is %d x %d, Z_new is %d x %d",
                   (int)nr, (int)nc, (int)Z_new.n_rows, (int)Z_new.n_cols);
    }

    uintptr_t data_ptr;
    memcpy(&data_ptr, base + g_ptr_offset, sizeof(uintptr_t));

    if (data_ptr == 0) Rcpp::stop("Null data pointer");
    if (data_ptr < 0x10000ULL) Rcpp::stop("Implausible data pointer");

    std::memcpy((void*)data_ptr, Z_new.memptr(), Z_new.n_elem * sizeof(double));
}