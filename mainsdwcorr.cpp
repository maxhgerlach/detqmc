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
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#include "cnpy.h"
#pragma GCC diagnostic pop
#include "exceptions.h"
#include "metadata.h"
#include "git-revision.h"
#include "statistics.h"
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

    bool operator==(const ConfigParameters& other) const {
        return (L == other.L) and
               (N == other.N) and
               (m == other.m) and
               (opdim == opdim) and
               (std::abs(dtau - other.dtau) < 1E-7);
    }
    bool operator!=(const ConfigParameters& other) const {
        return not operator==(other);
    }

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
        : vec_in (new VecCpx(conf_params.m)),
          vec_out(new VecCpx(conf_params.m)),
          mat_in (new MatCpx(conf_params.L, conf_params.L)),          
          mat_out(new MatCpx(conf_params.L, conf_params.L)),
          fft_temporal(vec_in, vec_out, +1),
          fft_spatial(mat_in, mat_out, -1)
    {
        vec_in->zeros();
        mat_in->zeros();
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


uintmax_t get_file_size(const std::string& filename) {
    return boost::filesystem::file_size(filename);
}

ConfigParameters get_conf_params(const std::string& metadata_filename) {
    ConfigParameters conf_params;
    MetadataMap meta = readOnlyMetadata(metadata_filename);
    getMeta(meta, "L", conf_params.L);
    getMeta(meta, "m", conf_params.m);
    getMeta(meta, "dtau", conf_params.dtau);
    getMeta(meta, "opdim", conf_params.opdim);
    conf_params.N = conf_params.L * conf_params.L;
    return conf_params;
}

ConfigParameters get_conf_params_for_directories(const std::vector< std::string >& input_directories) {
    namespace fs = boost::filesystem;
    ConfigParameters params;
    bool first = false;
    for (const auto& d : input_directories) {
        ConfigParameters d_params = get_conf_params((fs::path(d) / "info.dat").string());
        if (not first) {
            params = d_params;
            first = true;
        } else {
            if (params != d_params) {
                throw_GeneralError("configuration parameters do not agree");
            }
        }
    }
    return params;
}

std::vector<uintmax_t> get_sample_counts_for_directories(const std::vector< std::string >& input_directories,
                                                         const ConfigParameters& params) {
    namespace fs = boost::filesystem;
    // get number of samples for each directory
    std::vector<uintmax_t> sample_counts;
    uintmax_t sample_size = get_size_of_one_sample(params);
    for (const auto& d : input_directories) {
        std::string f = (fs::path(d) / "configs-phi.binarystream").string();
        uintmax_t file_size = get_file_size(f);
        if (file_size % (sample_size * sizeof(double)) != 0) {
            throw_GeneralError("unexpected binarystream file size");
        }
        uintmax_t this_sample_count = file_size / (sample_size * sizeof(double));
        sample_counts.push_back(this_sample_count);
    }
    return sample_counts;
}


void debug_print_slice(const PhiCorrelations& cube, uint32_t slc) {
    cube.slice(slc).eval().print();
}


// main entry for work
void process(const std::vector< std::string >& input_directories,
             const std::string& output_directory,
             uint32_t discard = 0, uint32_t jkblocks = 1) {
    if (jkblocks == 0) jkblocks = 1;

    namespace fs = boost::filesystem;

    uintmax_t dir_count = input_directories.size();

    ConfigParameters params = get_conf_params_for_directories(input_directories);

    std::vector<uintmax_t> sample_counts = get_sample_counts_for_directories(
        input_directories, params);

    // we leave out directories where we do not have more than $discard samples
    std::vector<uintmax_t> effective_sample_counts;
    uintmax_t total_sample_count = 0;
    for (uint32_t i = 0; i < dir_count; ++i) {
        uintmax_t effective_count = 0;
        if (sample_counts[i] > discard) {
            effective_count = sample_counts[i] - discard;
        }
        effective_sample_counts.push_back(effective_count);
        total_sample_count += effective_count;
    }

    // set up fft workspace
    FFT_workspace fft(params);

    // go through all input directories one after the other and
    // compute correlations for one sample after the other.  Find out,
    // what overall jackknife block we are at, and do the averages
    // correspondingly.
    uintmax_t jkblock_size = total_sample_count / jkblocks;
    // if jkblocks is not a divisor or total_sample_count, some data at the end will be discarded (just a small bit)
    uintmax_t total_sample_count_jk = jkblocks * jkblock_size;
    PhiConfig cur_config;
    PhiCorrelations cur_corr_ft;
    std::vector<PhiCorrelations> jkblock_corr_ft(
        jkblocks, arma::zeros<PhiCorrelations>(params.L, params.L, params.m));
    uintmax_t effective_sample_counter = 0;
    
    for (uint32_t i = 0; i < dir_count; ++i) {
        if (effective_sample_counts[i] == 0) {
            continue;
        }
        
        std::string d = input_directories[i];
        std::string f = (fs::path(d) / "configs-phi.binarystream").string();
        std::ifstream binary_float_input(f.c_str(),
                                         std::ios::in | std::ios::binary);
        if (not binary_float_input) {
            throw_ReadError(f);
        }
        uintmax_t this_directory_sample_counter = 0;
        while (readSystemConfiguration(cur_config, binary_float_input, params)) {
            if (effective_sample_counter > total_sample_count_jk) {
                // too many samples, skip the remaining
                break;
            }
            if (this_directory_sample_counter < discard) {
                // discard some initial configurations
                continue;
            } else {
                computeCorrelations_fft(cur_corr_ft, cur_config,
                                        params, fft);
                if (jkblocks > 1) {
                    uintmax_t cur_block = effective_sample_counter / jkblock_size;
                    for (uint32_t jb = 0; jb < jkblocks; ++jb) {
                        if (jb != cur_block) {
                            jkblock_corr_ft[jb] += cur_corr_ft;
                        }
                    }
                } else {
                    jkblock_corr_ft[0] += cur_corr_ft;
                }

                ++effective_sample_counter;
            }
            ++this_directory_sample_counter;
        }
    }
    
    // average and error bars
    PhiCorrelations avg_corr_ft, err_corr_ft;
    
    if (jkblocks > 1) {
        for (auto& single_block_corr_ft : jkblock_corr_ft) {
            single_block_corr_ft /= double((jkblocks - 1) * jkblock_size);
        }
        jackknife(avg_corr_ft, err_corr_ft,
                  jkblock_corr_ft,
                  arma::zeros<PhiCorrelations>(params.L, params.L, params.m).eval());
    } else {
        avg_corr_ft = jkblock_corr_ft[0] / double(total_sample_count_jk);
        err_corr_ft.zeros(params.L, params.L, params.m);
    }

    // wavevector and frequency values matching the Fourier transforms
    VecNum k_values = get_k_values(params.L, 1.0);
    unsigned int k_values_shape[] = {k_values.n_elem};
    VecNum omega_values = get_k_values(params.m, params.dtau);
    unsigned int omega_values_shape[] = {omega_values.n_elem};

    
    // Save results to a Numpy npz file, first convert cubes to
    // C-order.  The package Cnpy does not readily support writing
    // C-ordered data.
    PhiCorrelations avg_corr_ft_c_ordered = transpose_3d(avg_corr_ft);
    PhiCorrelations err_corr_ft_c_ordered = transpose_3d(err_corr_ft);
    unsigned int corr_ft_c_ordered_shape[] = {avg_corr_ft.n_slices, avg_corr_ft.n_cols, avg_corr_ft.n_rows};
    fs::path od(output_directory);
    fs::create_directories(od);
    std::string f = (od / "corr_ft.npz").string(); 
    cnpy::npz_save(f, "k_values",
                   k_values.memptr(), k_values_shape, 1, "w");
    cnpy::npz_save(f, "omega_values",
                   omega_values.memptr(), omega_values_shape, 1, "a");
    cnpy::npz_save(f, "avg_corr_ft__ky_kx_omega", 
                   avg_corr_ft_c_ordered.memptr(), corr_ft_c_ordered_shape, 3, "a");
    cnpy::npz_save(f, "err_corr_ft__ky_kx_omega", 
                   err_corr_ft_c_ordered.memptr(), corr_ft_c_ordered_shape, 3, "a");

    // writeout some info
    std::string finfo = (od / "corr_ft-info.dat").string();
    MetadataMap info;
    info["sdwcorr-discard"] = numToString(discard);
    info["sdwcorr-jkblocks"] = numToString(jkblocks);
    info["sdwcorr-totalsamples"] = numToString(total_sample_count_jk);
    writeOnlyMetaData(finfo, info, "sdwcorr info");
}


int main(int argc, char *argv[]) {
    std::vector< std::string > input_directories;
    std::string output_directory;
    uint32_t discard = 0;
    uint32_t jkblocks = 1;
    
    //parse command line options
    namespace po = boost::program_options;
    po::options_description options("Options for extraction of data from config binarystream");
    options.add_options()
        ("help", "print help on allowed options and exit")
        ("version,v", "print version information (git hash, build date) and exit")
        ("test", "run a simple test routine")
        ("input", po::value< std::vector< std::string > >(&input_directories)->multitoken(),
         "list of directories cotaining input data (multiple simindex for the same data point)")
        ("output", po::value< std::string >(&output_directory),
         "directory where to put results")
        ("discard,d", po::value<uint32_t>(&discard)->default_value(0),
         "number of initial configuration samples to discard (additional thermalization)")
        ("jkblocks,j", po::value<uint32_t>(&jkblocks)->default_value(1),
         "number jackknife blocks for estimating error bars")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    //handle simple options
    if (vm.count("help")) {
        std::cout << options << std::endl;
        return 0;
    }
    if (vm.count("version")) {
        std::cout << "Build info:\n"
                  << metadataToString(collectVersionInfo())
                  << std::endl;
        return 0;
    }
    if (vm.count("test")) {
        test_corr_ft();
        return 0;
    }

    if (not (vm.count("input") or vm.count("output"))) {
        std::cerr << "pass --input and --output options\n";
        return 1;
    }
         
    std::cout << "processing input: ";
    for (const auto& s : input_directories) {
        std::cout << s << " ";
    }
    std::cout << std::endl;

    process(input_directories, output_directory, discard, jkblocks);
    
    return 0;
}

