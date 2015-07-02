#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#include <chrono>
#include <memory>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <armadillo>
#include <fftw3.h>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/program_options.hpp"
#pragma GCC diagnostic pop
#include "exceptions.h"
#include "metadata.h"
#include "git-revision.h"
#include "tools.h"


typedef double num;
typedef std::complex<double> cpx;
typedef arma::Col<num> VecNum;
typedef arma::Col<cpx> VecCpx;
typedef arma::Mat<num> MatNum;
typedef arma::Mat<cpx> MatCpx;
typedef arma::Cube<num> CubeNum;
typedef arma::Cube<cpx> CubeCpx;


// slice indexes timeslice, column indexes order parameter
// dimension, row indexes site    
// [Index order: row, col, slice]
// [Data layout; slice after slice, matrices then column major] 
typedef CubeNum PhiConfig;


// TODO: clarify index values

// for real space:
// < phi_(iy_delta, ix_delta)(nt_delta) phi_(0,0)(0) >
// row indexes site-y-delta, column indexes site-x-delta,
// slice indexes timeslice-delta

// for Fourier space:
// row indexes ky, column indexes kx, slice indexes omega

// [Index order: row, col, slice]
typedef CubeNum PhiCorrelations;

struct ConfigParameters {
    uint32_t L;                 // linear system size
    // uint32_t d;                 // lattice dimension, should be 2
    uint32_t N;                 // total lattice volume
    uint32_t m;                 // number of timeslices
    num dtau;                   // beta = m * dtau
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
            // in the files on disk we stored time slices 1 ... m,
            // for our calculations here it is easier to consider time slices,
            // 0 ... m-1, where, due to periodic boundary conditions, slice m ==
            // slice 0
            phi_target.slice(0) = phi_target.slice(conf_params.m);
            phi_target.shed_slice(conf_params.m);
            
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


// compute positive delta_x = x2 - x1 on a periodic ring of length
// ring_length, 0 <= delta_x < ring_length -- so this gives 'too long'
// distances.  This can directly be used as real space index.
inline int32_t pos_pbc_dist(int32_t x1, int32_t x2, int32_t ring_length) {
    // difference restricted to pbc ring
    int32_t delta_x = (x2 - x1) % ring_length;
    if (delta_x < 0) {
        delta_x += ring_length;
    }
    assert(delta_x >= 0);
    assert(delta_x < ring_length);
    return delta_x;
}


// // real space index
// inline uint32_t pbc_diff_to_array_index(int32_t pbc_diff, int32_t ring_length) {
//     assert(-ring_length/2 < pbc_diff);
//     assert(pbc_diff <= ring_length/2);
//     int32_t index = pbc_diff + ring_length / 2 - 1;
//     assert(index >= 0);
//     assert(index <  ring_length);
//     return static_cast<uint32_t>(index);
// }


// give meaning to array indices (1)
// 0 ... (L-1) * spacing
VecNum get_pos_pbc_dist_values(uint32_t ring_length, num spacing=1.0) {
    assert(ring_length > 0);
    assert(spacing > 0);
    VecNum dists(ring_length);
    for (uint32_t counter = 0; counter < ring_length; ++counter) {
        dists[counter] = counter * spacing;
    }
    return dists;
}

// VecNum get_pbc_diff_values(uint32_t ring_length, num spacing=1.0) {
//     assert(ring_length > 0);
//     VecNum pbc_diffs(ring_length);
//     // this was intended for even ring_length, but seems to work fine
//     // with odd ring_length too
//     pbc_diffs[ring_length - 1] = spacing * (ring_length / 2);
//     for (int32_t counter = ring_length - 2; counter >= 0; --counter) {
//         pbc_diffs[counter] = pbc_diffs[counter + 1] - spacing;
//     }
//     return pbc_diffs;
// }



// give meaning to array indices (2)
// k values: 0 ... (L-1)*2pi / (L*spacing)
VecNum get_k_values(uint32_t ring_length, num spacing=1.0) {
    assert(ring_length > 0);
    const auto pi = arma::datum::pi;
    VecNum k_values(ring_length);
    k_values[0] = 0.0;
    for (uint32_t counter = 1; counter < ring_length; ++counter) {
        k_values[counter] = k_values[counter - 1] + 2*pi / (ring_length * spacing);
    }
    return k_values;
}




// to speed things up do not use translational invariance -- needs way
// more samples
void computeCorrelations_reduced(PhiCorrelations& corr_target, const PhiConfig& conf,
                                 const ConfigParameters& conf_params) {
    corr_target.zeros(conf_params.L, // rows ... site y delta
                      conf_params.L, // cols ... site x delta
                      conf_params.m  // slcs ... time slice delta
        );
    uint32_t i2y = 0;
    uint32_t i2x = 0;
    uint32_t i2 = i2y * conf_params.L + i2x;
    uint32_t nt2 = 1;
    for (uint32_t nt1 = 0; nt1 < conf_params.m; ++nt1) {
        for (uint32_t i1y = 0; i1y < conf_params.L; ++i1y) {
            for (uint32_t i1x = 0; i1x < conf_params.L; ++i1x) {
                uint32_t i1 = i1y * conf_params.L + i1x;
                int32_t ix_delta = pos_pbc_dist(i1x, i2x, conf_params.L);
                int32_t iy_delta = pos_pbc_dist(i1y, i2y, conf_params.L);
                int32_t nt_delta = pos_pbc_dist(nt1, nt2, conf_params.m);
                            
                num dot_product = 0.0;
                for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
                    dot_product += conf(i1, dim, nt1) * conf(i2, dim, nt2);
                }

                // corr_target(pbc_diff_to_array_index(iy_delta, conf_params.L),
                //             pbc_diff_to_array_index(ix_delta, conf_params.L),
                //             pbc_diff_to_array_index(nt_delta,  conf_params.m))
                //     += dot_product;
                corr_target(iy_delta, ix_delta, nt_delta) += dot_product;
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
    for (uint32_t nt1 = 0; nt1 < conf_params.m; ++nt1) {
        for (uint32_t nt2 = 0; nt2 < conf_params.m; ++nt2) {
            for (uint32_t i1y = 0; i1y < conf_params.L; ++i1y) {
                for (uint32_t i1x = 0; i1x < conf_params.L; ++i1x) {
                    uint32_t i1 = i1y * conf_params.L + i1x;
                    for (uint32_t i2y = 0; i2y < conf_params.L; ++i2y) {
                        for (uint32_t i2x = 0; i2x < conf_params.L; ++i2x) {
                            uint32_t i2 = i2y * conf_params.L + i2x;

                            int32_t ix_delta = pos_pbc_dist(i1x, i2x, conf_params.L);
                            int32_t iy_delta = pos_pbc_dist(i1y, i2y, conf_params.L);
                            int32_t nt_delta = pos_pbc_dist(nt1, nt2, conf_params.m);
                            
                            num dot_product = 0.0;
                            for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
                                dot_product += conf(i1, dim, nt1) * conf(i2, dim, nt2);
                            }

                            // corr_target(pbc_diff_to_array_index(iy_delta, conf_params.L),
                            //             pbc_diff_to_array_index(ix_delta, conf_params.L),
                            //             pbc_diff_to_array_index(nt_delta, conf_params.m))
                            //     += dot_product;
                            corr_target(iy_delta, ix_delta, nt_delta) += dot_product;
                        }
                    }
                }                
            }
        }
    }

    corr_target /= num(conf_params.N * conf_params.m);
}

// this assumes corr(r) = corr(-r), so its Fourier transform is real
void slow_ft(PhiCorrelations& corr_ft, const PhiCorrelations& corr, const ConfigParameters& conf_params) {
    corr_ft.zeros(conf_params.L,
                  conf_params.L,
                  conf_params.m
        );

    VecNum k_values = get_k_values(conf_params.L, 1.0);
    VecNum omega_values = get_k_values(conf_params.m, conf_params.dtau);

    VecNum xy_values = get_pos_pbc_dist_values(conf_params.L, 1.0);
    VecNum t_values  = get_pos_pbc_dist_values(conf_params.m, conf_params.dtau);
    
    for (uint32_t omega_index = 0; omega_index < conf_params.m; ++omega_index) {
        num omega = omega_values[omega_index];
        for (uint32_t ky_index = 0; ky_index < conf_params.L; ++ky_index) {
            num ky = k_values[ky_index];
            for (uint32_t kx_index = 0; kx_index < conf_params.L; ++kx_index) {
                num kx = k_values[kx_index];
                for (uint32_t y_index = 0; y_index < conf_params.L; ++y_index) {
                    num y = xy_values[y_index];
                    for (uint32_t x_index = 0; x_index < conf_params.L; ++x_index) {
                        num x = xy_values[x_index];
                        for (uint32_t t_index = 0; t_index < conf_params.m; ++t_index) {
                            num t = t_values[t_index];

                            num argument = -(kx * x + ky * y) + t * omega;
                            cpx phase = std::exp(cpx(0, argument));

                            cpx contrib = phase * corr(y_index, x_index, t_index);
                            corr_ft(ky_index, kx_index, omega_index) += contrib.real();
                        }
                    }
                }
            }
        }
    }

    corr_ft /= (conf_params.L * conf_params.L * conf_params.m);
}

// use FFTW3 on Armadillo vector
class FFT1d {
    int n;
    fftw_complex *in, *out;
    fftw_plan my_plan;
public:
    // We use pointers to make sure the vectors use the same memory
    // all the time. Pass pointers to vectors of equal length.
    // Keep using the same vectors.    
    // sign == -1: FFTW_FORWARD
    // sign == +1: FFTW_BACKWARD
    FFT1d(std::shared_ptr<VecCpx> vec_in, std::shared_ptr<VecCpx> vec_out,
          int sign) {
        assert(vec_in->n_elem == vec_out->n_elem);
        n = vec_in->n_elem;
        in  = reinterpret_cast<fftw_complex*>(vec_in->memptr());
        out = reinterpret_cast<fftw_complex*>(vec_out->memptr());
        my_plan = fftw_plan_dft_1d(n, in, out, sign, FFTW_MEASURE);
    }
    
    void execute() {
        fftw_execute(my_plan);
    }
    
    ~FFT1d() {
        fftw_destroy_plan(my_plan);
    }    
};

class FFT2d {
    int rows, cols;
    int n;
    fftw_complex *in, *out;     // we will store column-major data in here (like BLAS/Lapack/Armadillo)
    fftw_plan my_plan;
public:
    // We use pointers to make sure the matrices use the same memory
    // all the time. Pass pointers to matrices of equal dimensions.
    // Keep using the same matrices.
    // sign == -1: FFTW_FORWARD
    // sign == +1: FFTW_BACKWARD
    FFT2d(std::shared_ptr<MatCpx> mat_in, std::shared_ptr<MatCpx> mat_out,
          int sign) {
        assert(mat_in->n_rows == mat_out->n_rows);
        assert(mat_in->n_cols == mat_out->n_cols);
        rows = mat_in->n_rows;
        cols = mat_in->n_cols;
        n = rows * cols;
        in  = reinterpret_cast<fftw_complex*>(mat_in->memptr());
        out = reinterpret_cast<fftw_complex*>(mat_out->memptr());
        // switch row & column dimensions: understand our column-major data
        // as row-major
        my_plan = fftw_plan_dft_2d(cols, rows, in, out, sign, FFTW_MEASURE);    
    }
    
    void execute() {
        fftw_execute(my_plan);
    }
    
    ~FFT2d() {
        fftw_destroy_plan(my_plan);
    }    
};


// Armadillo seems to lack the necessary overload
void vecCpx_set_real_from_cubeNum_tube(VecCpx& vec, const CubeNum& cube,
                                       uint32_t row, uint32_t col) {
    for (uint32_t slice = 0; slice < cube.n_slices; ++slice) {
        vec[slice].real(cube(row, col, slice));
    }
}


// temporal and spatial FFT set up
struct FFT_workspace {
    std::shared_ptr<VecCpx> vec_in, vec_out;
    std::shared_ptr<MatCpx> mat_in, mat_out;
    FFT1d fft_temporal;
    FFT2d fft_spatial;

    FFT_workspace(const ConfigParameters& conf_params)
        : vec_in(new VecCpx(conf_params.m)),
          vec_out(new VecCpx(conf_params.m)),
          mat_in (new MatCpx(conf_params.L, conf_params.L)),          
          mat_out(new MatCpx(conf_params.L, conf_params.L)),
          fft_temporal(vec_in, vec_out, +1),
          fft_spatial(mat_in, mat_out, -1)
    {
        vec_in ->zeros();
        mat_in ->zeros();
    }
    
};


// compute Fourier transformed correlations for one system configuration
void computeCorrelations_fft(PhiCorrelations& corr_ft, const PhiConfig& conf,
                             const ConfigParameters& conf_params,
                             FFT_workspace& fft) {
    // complex intermediary
    CubeCpx phi_ft(conf_params.L, conf_params.L, conf_params.m);

    // real result
    corr_ft.zeros(conf_params.L, conf_params.L, conf_params.m);
    
    // FFTs in time and space for each order parameter dimension
    for (uint32_t dim = 0; dim < conf_params.opdim; ++dim) {
        phi_ft.zeros();

        // temporal
        for (uint32_t site = 0; site < conf_params.N; ++site) {
            // instead of
            //    vec_in->set_real( conf.tube(site, dim) );
            // do this:
            vecCpx_set_real_from_cubeNum_tube(*(fft.vec_in), conf, site, dim);
            fft.fft_temporal.execute();
            uint32_t y_index = site / conf_params.L;
            uint32_t x_index = site % conf_params.L;
            phi_ft.tube(y_index, x_index) = *(fft.vec_out);
        }

        // spatial
        for (uint32_t nt = 0; nt < conf_params.m; ++nt) {
            *(fft.mat_in) = phi_ft.slice(nt); // copies into mat_in's memory
            fft.fft_spatial.execute();
            phi_ft.slice(nt) = *(fft.mat_out);
        }

        // the FT'ed correlation function is given by the squared modulus of the 
        // FT'ed spin correlation function
        corr_ft += arma::real(phi_ft % arma::conj(phi_ft));
    }
    
    // normalization of Fourier transforms was missing, correct for that:
    corr_ft /= (conf_params.N * conf_params.N * conf_params.m * conf_params.m); 
}


void test_corr_ft() {
    ConfigParameters params;
    params.L = 10;
    params.N = params.L * params.L;
    params.m = 200;
    params.dtau = 0.1;
    params.opdim = 2;
    PhiConfig config = arma::randu<PhiConfig>(params.N, params.opdim, params.m);

    PhiCorrelations corr_fft = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    FFT_workspace fft(params);
    
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    computeCorrelations_fft(corr_fft, config, params, fft);
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "corr_fft: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " microseconds" << std::endl;

    PhiCorrelations corr_full = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    PhiCorrelations corr_full_ft = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    
    start = std::chrono::steady_clock::now();
    computeCorrelations_full(corr_full, config, params);
    slow_ft(corr_full_ft, corr_full, params);
    end = std::chrono::steady_clock::now();
    std::cout << "corr_full & corr_full_ft: "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
              << " microseconds" << std::endl;    

    PhiCorrelations diff_corr_ft = arma::abs((corr_full_ft - corr_fft) / corr_full_ft);

    std::cout << "full_ft vs fft --- max: " << diff_corr_ft.max()
              << ", mean: " << arma::accu(diff_corr_ft) / diff_corr_ft.n_elem << std::endl;


    // PhiCorrelations corr_reduced = arma::zeros<PhiCorrelations>(params.L, params.L, params.m);
    // computeCorrelations_reduced(corr_reduced, config, params);
    // PhiCorrelations diff = arma::abs((corr_full - corr_reduced) / corr_full);
    // std::cout << "full vs reduced --- max: " << diff.max() << ", mean: " << arma::accu(diff) / diff.n_elem << std::endl;
    // // for a single purely random sample there will of course be a
    // // significant difference between _reduced and _full
    
}


int main(int argc, char *argv[]) {
    test_corr_ft();

    return 0;
}

