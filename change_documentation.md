# Harmony 2.0.3 Compatibility Fix for `align_harmony`

## Background

`align_harmony` implements a custom iterative batch correction loop where LEMUR alternates between:
1. **Harmony** computing maximally diverse soft cluster assignments (`R` matrix)
2. **LEMUR** computing a geometrically-aware embedding correction (`align_impl`)

The corrected embedding is then fed back into harmony for the next clustering iteration.
This requires writing LEMUR's corrected embedding into harmony's internal `Z_corr` field
between iterations.

---

## What broke

In harmony **1.2.4**, `Z_corr` was registered as a readable/writable field in the Rcpp module:

```cpp
// harmony 1.2.4 src/harmony.cpp
RCPP_MODULE(harmony_module) {
    class_<harmony>("harmony")
        .field("Z_corr", &harmony::Z_corr)  // ← read/write
        .field("Z_cos",  &harmony::Z_cos)
        ...
}
```

So `harm_obj$Z_corr <- my_matrix` worked directly from R.

In harmony **2.0.3**, these were replaced with getter methods only:

```cpp
// harmony 2.0.3
.method("getZcorr", &harmony::getZcorr)  // ← read only, no setter
```

The old code called `harm_obj$setZcorr(alignment$embedding)` which no longer exists,
causing `align_harmony` to fail entirely.

---

## Why a simple fix wasn't possible

The obvious fixes were all blocked:

- **`harm_obj$Z_corr <- x`**: fails — not registered as writable field in 2.0.3
- **`harm_obj$setZcorr(x)`**: doesn't exist — maintainer declined to add it
- **`LinkingTo: harmony`**: harmony's `inst/include/` is empty — no headers shipped
- **Shipping a modified harmony**: not acceptable for CRAN distribution
- **Reinitializing harmony each iteration**: breaks convergence tracking since
  `objective_harmony` history is reset, and discards accumulated clustering state

---

## What we did instead

We inject a setter at runtime by directly writing to harmony's C++ object memory.
This works because:

1. The harmony R object (`Rcpp_harmony`) is an S4 wrapper around a raw C++ pointer
   (`ho$.pointer`)
2. `Z_corr` is a public `arma::mat` field in the C++ object — there are no access
   restrictions at the binary level, only at the Rcpp module registration level
3. An `arma::mat` stores its data in a heap-allocated buffer pointed to by an internal
   pointer. We can overwrite that buffer in place without touching `arma::mat`'s
   bookkeeping

### Finding the correct memory offsets

This is the core technical challenge. We cannot use `offsetof()` on the harmony class
because harmony ships no headers. Instead we disassemble harmony's own `getZcorr()`
function at runtime to read the exact byte offsets the compiler used:

```r
# Disassemble getZcorr to find:
# - offset of Z_corr's n_rows/n_cols within the harmony object
# - offset of Z_corr's data pointer within the harmony object
asm <- system(paste("objdump -d [--macho on macOS]", so_file,
                    "| grep -A 80 '_ZN7harmony8getZcorrEv'"), intern = TRUE)
```

From the ARM64 disassembly of harmony 2.0.3: