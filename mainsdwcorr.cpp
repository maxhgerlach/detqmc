#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/program_options.hpp"
#pragma GCC diagnostic pop
#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include "exceptions.h"
#include "metadata.h"
#include "git-revision.h"
#include "tools.h"


typedef double num;
typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;

// slice indexes timeslice, column indexes order parameter
// dimension, row indexes site    
// [Index order: row, col, slice]
// [Data layout; slice after slice, matrices then column major] 
typedef CubeNum PhiConfig;

// < phi_(iy_delta, ix_delta)(k_delta) phi_(0,0)(0) >
// row indexes site-y-delta, column indexes site-x-delta,
// slice indexes timeslice-delta
// [Index order: row, col, slice]
typedef CubeNum PhiCorrelations;

struct ConfigParameters {
    uint32_t L;                 // linear system size
    // uint32_t d;                 // lattice dimension, should be 2
    uint32_t N;                 // total lattice volume
    uint32_t m;                 // number of timeslices
    uint32_t opdim;             // 1, 2, or 3
};

// return size of sample in units of doubles
uint32_t get_size_of_one_sample(const ConfigParameters& conf_params) {
    // uint32_t N = uint_pow(conf_params.L, conf_params.d);
    return conf_params.N * conf_params.m * conf_params.opdim;
}


// returns false on failure
bool readSystemConfiguration(PhiConfig& phi_target, std::ifstream& binary_float_input,
                             const ConfigParameters& conf_params) {
    bool return_success = false;
    uint32_t sample_size = get_size_of_one_sample(conf_params);
    if (binary_float_input) {
        std::vector<double> one_sample;
        one_sample.resize(sample_size, 0.0);
        binary_float_input.read(reinterpret_cast<char*>(&(*one_sample.begin())), sizeof(double) * sample_size);
        if (binary_float_input) {
            // no failure: store configuration
            phi_target.set_size(conf_params.N,     // rows ... sites
                                conf_params.opdim, // cols ... op dim
                                conf_params.m + 1  // slcs ... time slices
                );
            uint32_t n = 0;
            for (uint32_t ix = 0; ix < conf_params.L; ++ix) {
                for (uint32_t iy = 0; iy < conf_params.L; ++iy) {
                    uint32_t i = iy*conf_params.L + ix;
                    for (uint32_t k = 1; k <= conf_params.m; ++k) {
                        for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
                            phi_target(i, dim, k) = one_sample[n];
                            ++n;
                        }
                    }
                }
            }
            phi_target.slice(0).zeros();
            return_success  = true;
        }
    }
    return return_success;
}


// compute delta_x = x2 - x1 on a periodic ring of length ring_length,
// -L/2 < delta_x <= +L/2
inline int32_t pbc_diff(int32_t x1, int32_t x2, int32_t ring_length) {
    // difference restricted to pbc ring
    int32_t delta_x = (x2 - x1) % ring_length;
    if (delta_x < 0) {
        delta_x += ring_length;
    }
    assert(delta_x >= 0);
    assert(delta_x < ring_length);
    // restrict to range -L/2 < delta_x <= +L/2
    if (delta_x > ring_length / 2) {
        delta_x -= ring_length;
    }
    assert(-ring_length/2 < delta_x);
    assert(delta_x <= ring_length/2);
    return delta_x;
}

inline uint32_t pbc_diff_to_array_index(int32_t pbc_diff, int32_t ring_length) {
    assert(-ring_length/2 < pbc_diff);
    assert(pbc_diff <= ring_length/2);
    int32_t index = pbc_diff + ring_length / 2 - 1;
    assert(index >= 0);
    assert(index <  ring_length);
    return static_cast<uint32_t>(index);
}


// to speed things up do not use translational invariance
void computeCorrelations_reduced(PhiCorrelations& corr_target, const PhiConfig& conf,
                                 const ConfigParameters& conf_params) {
    corr_target.zeros(conf_params.L, // rows ... site y delta
                      conf_params.L, // cols ... site x delta
                      conf_params.m  // slcs ... time slice delta
        );
    uint32_t i2y = 0;
    uint32_t i2x = 0;
    uint32_t i2 = i2y * conf_params.L + i2x;
    uint32_t k2 = 1;
    for (uint32_t k1 = 1; k1 <= conf_params.m; ++k1) {
        for (uint32_t i1y = 0; i1y < conf_params.L; ++i1y) {
            for (uint32_t i1x = 0; i1x < conf_params.L; ++i1x) {
                uint32_t i1 = i1y * conf_params.L + i1x;
                int32_t ix_delta = pbc_diff(i1x, i2x, conf_params.L);
                int32_t iy_delta = pbc_diff(i1y, i2y, conf_params.L);
                int32_t k_delta  = pbc_diff(k1,  k2,  conf_params.m);
                            
                num dot_product = 0.0;
                for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
                    dot_product += conf(i1, dim, k1) * conf(i2, dim, k2);
                }

                corr_target(pbc_diff_to_array_index(iy_delta, conf_params.L),
                            pbc_diff_to_array_index(ix_delta, conf_params.L),
                            pbc_diff_to_array_index(k_delta,  conf_params.m))
                    += dot_product;
            }
        }
    }
}



void computeCorrelations_full(PhiCorrelations& corr_target, const PhiConfig& conf,
                              const ConfigParameters& conf_params) {
    corr_target.zeros(conf_params.L, // rows ... site y delta
                      conf_params.L, // cols ... site x delta
                      conf_params.m  // slcs ... time slice delta
        );
    for (uint32_t k1 = 1; k1 <= conf_params.m; ++k1) {
        for (uint32_t k2 = 1; k2 <= conf_params.m; ++k2) {
            for (uint32_t i1y = 0; i1y < conf_params.L; ++i1y) {
                for (uint32_t i1x = 0; i1x < conf_params.L; ++i1x) {
                    uint32_t i1 = i1y * conf_params.L + i1x;
                    for (uint32_t i2y = 0; i2y < conf_params.L; ++i2y) {
                        for (uint32_t i2x = 0; i2x < conf_params.L; ++i2x) {
                            uint32_t i2 = i2y * conf_params.L + i2x;

                            // int32_t ix_delta = (i1x - i2x) % conf_params.L;
                            // int32_t iy_delta = (i1y - i2y) % conf_params.L;
                            // int32_t k_delta  = (k1  - k2)  % conf_params.m;;

                            // int32_t ix_delta = i2x - i1x;
                            // int32_t iy_delta = i2y - i1y;
                            // int32_t k_delta  = k2  - k1;

                            int32_t ix_delta = pbc_diff(i1x, i2x, conf_params.L);
                            int32_t iy_delta = pbc_diff(i1y, i2y, conf_params.L);
                            int32_t k_delta  = pbc_diff(k1,  k2,  conf_params.m);
                            
                            num dot_product = 0.0;
                            for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
                                dot_product += conf(i1, dim, k1) * conf(i2, dim, k2);
                            }

                            corr_target(pbc_diff_to_array_index(iy_delta, conf_params.L),
                                        pbc_diff_to_array_index(ix_delta, conf_params.L),
                                        pbc_diff_to_array_index(k_delta,  conf_params.m))
                                += dot_product;
                        }
                    }
                }                
            }
        }
    }

    corr_target /= num(conf_params.N * conf_params.m);
}



int main(int argc, char *argv[]) {
    // test
    ConfigParameters params;
//    assert(2 > 3); // assert really are ignored in our debug builds
    params.L = 10;
    params.N = 100;
    params.m = 200;
    params.opdim = 2;
    PhiConfig config = arma::randu<PhiConfig>(params.N, params.opdim, params.m + 1);
    PhiCorrelations corr_full = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    PhiCorrelations corr_reduced = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    computeCorrelations_full(corr_full, config, params);
    computeCorrelations_reduced(corr_reduced, config, params);
    PhiCorrelations diff = arma::abs((corr_full - corr_reduced) / corr_full);
    std::cout << "max: " << diff.max() << ", mean: " << arma::accu(diff) / diff.n_elem << std::endl;
    // for a single purely random sample there will of course be a
    // significant difference between _reduced and _full
    return 0;
}
