/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

// Compute bosonic equal-time, q=0 susceptibility from stored system
// configurations, write out as a time series


#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#include <memory>
#include <complex>
#include <fstream>
#include <iostream>
#include <string>
#include <armadillo>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#include "boost/filesystem.hpp"
#include "boost/program_options.hpp"
#pragma GCC diagnostic pop
#include "dataserieswritersucc.h"
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

//return path to configs-phi.binarystream or extracted-configs-phi.binarystream
//if more appropriate
boost::filesystem::path get_input_file_path(const std::string& input_directory) {
    namespace fs = boost::filesystem;
    fs::path p_result;
    fs::path p_1 = fs::path(input_directory) / "configs-phi.binarystream";
    fs::path p_2 = fs::path(input_directory) / "extracted-configs-phi.binarystream";
    if (fs::exists(p_1)) {
        p_result = p_1;
    } else {
        if (not fs::exists(p_2)) {
            throw_GeneralError("No binary configuration stream file found");
        }
        p_result = p_2;
    }
    return p_result;
}

ConfigParameters get_conf_params_for_directory(const std::string& input_directory) {
    namespace fs = boost::filesystem;
    ConfigParameters params = get_conf_params((fs::path(input_directory) / "info.dat").string());
    return params;
}


MetadataMap get_metadata_for_directory(const std::string& input_directory) {
    namespace fs = boost::filesystem;
    return readOnlyMetadata((fs::path(input_directory) / "info.dat").string());
}


std::size_t get_sample_count_for_directory(const std::string& input_directory,
                                         const ConfigParameters& params) {
    namespace fs = boost::filesystem;
    std::size_t sample_size = get_size_of_one_sample(params);
    std::string f = get_input_file_path(input_directory).string();
    std::size_t file_size = static_cast<std::size_t>(get_file_size(f));
    
    if (file_size % (sample_size * sizeof(double)) != 0) {
        throw_GeneralError("unexpected binarystream file size");
    }
    std::size_t sample_count = file_size / (sample_size * sizeof(double));
    
    return sample_count;
}

num computePhiSusceptibilityEqualTime(const PhiConfig& config,
                                      const ConfigParameters& params) {
    num eqtime_susc = 0.0;
    for (uint32_t timeslice = 0; timeslice < params.m; ++timeslice) {
        VecNum eqtime_mag(params.opdim, arma::fill::zeros);
        for (uint32_t dim = 0; dim < params.opdim; ++dim) {
            for (uint32_t site = 0; site < params.N; ++site) {
                eqtime_mag[dim] += config(site, dim, timeslice);
            }
        }
        eqtime_mag /= num(params.N);
        eqtime_susc += arma::dot(eqtime_mag, eqtime_mag);
    }
    eqtime_susc *= (num(params.N) / num(params.m));
    return eqtime_susc;
}


// main entry for work
void process() {
    namespace fs = boost::filesystem;

    MetadataMap meta = get_metadata_for_directory(".");
    ConfigParameters params = get_conf_params_for_directory(".");

    std::size_t sample_count = get_sample_count_for_directory(".", params);

    std::vector<num> timeseries;
    timeseries.reserve(sample_count);

    // read in sample after sample and buffer in time series
    std::string f = get_input_file_path(".").string();
    std::ifstream binary_float_input(f.c_str(), std::ios::in | std::ios::binary);
    if (not binary_float_input) {
        throw_ReadError(f);
    }
    PhiConfig cur_config;
    uintmax_t sample_counter = 0;
    while (readSystemConfiguration(cur_config, binary_float_input, params)) {
        ++sample_counter;

        timeseries.push_back(computePhiSusceptibilityEqualTime(cur_config, params));
    }

    // write out everything using DataSeriesWriter
    const std::string observable = "phiSusceptibilityEqualTime";
    DoubleVectorWriterSuccessive writer(observable + ".series", false);
    writer.addMetadataMap(meta);
    writer.addMeta("observable", observable);
    writer.addHeaderText("Computed using sdweqtimesusc");
    writer.writeHeader();
    writer.writeData(timeseries);
}



int main(int argc, char *argv[]) {
    //parse command line options
    namespace po = boost::program_options;
    po::options_description options("Options for extraction of data from config binarystream to compute bosonic equal-time q=0 susceptibility (input and output in current working directory)");
    options.add_options()
        ("help", "print help on allowed options and exit")
        ("version,v", "print version information (git hash, build date) and exit")
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
         
    std::cout << "processing input" << std::endl;

    process();
    
    return 0;
}
