// src/harmony_defines.h
#pragma once
#include <RcppArmadillo.h>
#include <vector>

typedef arma::mat    MATTYPE;
typedef arma::sp_mat SPMAT;
typedef arma::vec    VECTYPE;

class harmony_layout {
public:
    // 32 bytes of fields present in harmony 2.0.3 but not in the
    // header we have access to. Empirically determined: offsetof(Z_corr)
    // should be 0x180=384, our layout gives 352, difference = 32 bytes.
    uint8_t _padding[32];
    
    // Fields from harmony 2.0.3 src/harmony.h in original order
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y;
    SPMAT Phi, Phi_moe, Phi_moe_t, Phi_t, Rk;
    VECTYPE Pr_b, theta, N_b, sigma, lambda;
    std::vector<float> objective_kmeans, objective_kmeans_dist,
                       objective_kmeans_entropy, objective_kmeans_cross,
                       objective_harmony;
    std::vector<int> kmeans_rounds, B_vec, covariate_bounds;
    std::vector<arma::uvec> index;
    arma::uvec batch_sizes, new_index, original_index, batch_indptr;
    float block_size, epsilon_kmeans, epsilon_harmony, alpha,
          batch_proportion_cutoff;
    unsigned int N, K, B, d, max_iter_kmeans, window_size;
    MATTYPE W, dist_mat, O, E, dir_prior;
    arma::uvec update_order, cells_update;
    bool ran_setup, ran_init, l