#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

// discard some initial samples, then discard all but every n'th
// sample of system configurations in a binarystream, write that to a
// new binarystream

// input expected in working directory:
//    configs-phi.infoheader    ... metadata about this timeseries
//    configs-phi.binarystream  ... phi configurations in binary doubles
//
// output written into working directory:
//    extracted-configs-phi.infoheader
//    extracted-configs-phi.binarystream
//
// commandline parameters
//    -d M   ...  M: number of samples to discard at beginning
//    -s N   ...  then take only every N'th sample


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/program_options.hpp"
#pragma GCC diagnostic pop
#include <fstream>
#include <iostream>
#include <string>
#include "exceptions.h"
#include "metadata.h"
#include "git-revision.h"
#include "tools.h"

// returns number of samples written to output file
uint32_t extract(const std::string& input_data_file,
                 const std::string& output_data_file,
                 uint32_t one_sample_size, // number of doubles in one sample
                 uint32_t discard, uint32_t subsample_interval) {
    
    std::ifstream binary_float_input(input_data_file.c_str(), std::ios::in | std::ios::binary);
    if (not binary_float_input) {
        throw_ReadError(input_data_file);
    }

    std::ofstream binary_float_output(output_data_file.c_str(), std::ios::out | std::ios::binary);

    uint32_t read_length = one_sample_size * sizeof(double);
    std::string sample_buffer;
    sample_buffer.resize(read_length, '\0');
    char* sample_buffer_begin = &(*sample_buffer.begin());

    uint32_t sample_counter = 0;
    uint32_t extracted_counter = 0;
    while (binary_float_input) {
        binary_float_input.read(sample_buffer_begin, read_length);
        uint32_t actually_read = uint32_t(binary_float_input.gcount());
        if (actually_read == 0) {
            // nothing left
            continue;
        }
        else if (actually_read != read_length) {
            throw_GeneralError("Could not read " + numToString(read_length) +
                               " characters from file: actually read " +
                               numToString(actually_read) + " characters.");
        }
        ++sample_counter;
        if (sample_counter <= discard) {
            continue;
        } else if (((sample_counter - discard - 1) % subsample_interval) != 0) {
            continue;
        } else {
            binary_float_output.write(sample_buffer_begin, read_length);
            ++extracted_counter;
        }
    }

    return extracted_counter;    
}

// return size of sample in units of doubles
uint32_t get_size_of_one_sample(const MetadataMap& meta) {
    uint32_t N = fromString<uint32_t>(meta.at("N"));
    uint32_t m = fromString<uint32_t>(meta.at("m"));
    uint32_t opdim = fromString<uint32_t>(meta.at("opdim"));
    return N * m * opdim;
}

// write to extracted-configs-phi.infoheader
void write_info(const std::string& out_infoheader_filename,
                const MetadataMap& initial_meta,
                uint32_t discard, uint32_t subsample_interval, uint32_t extracted_samples) {
    MetadataMap meta = initial_meta;
    meta["extract-discard"] = numToString(discard);
    meta["extract-subsample"] = numToString(subsample_interval);
    meta["extract-samples"] = numToString(extracted_samples);
    writeOnlyMetaData(out_infoheader_filename, meta);
}

int main(int argc, char *argv[]) {
    uint32_t discard = 0;
    uint32_t subsample_interval = 1;

    //parse command line options
    namespace po = boost::program_options;
    po::options_description extractOptions("Options for extraction of data from config binarystream");
    extractOptions.add_options()
        ("help", "print help on allowed options and exit")
        ("version,v", "print version information (git hash, build date) and exit")
        ("discard,d", po::value<uint32_t>(&discard)->default_value(0),
         "number of initial configuration samples to discard (additional thermalization)")
        ("subsample,s", po::value<uint32_t>(&subsample_interval)->default_value(1),
         "only write out every s'th sample of the input")
        ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, extractOptions), vm);
    po::notify(vm);

    //handle simple options
    if (vm.count("help")) {
        std::cout << " discard some initial samples, then discard all but every n'th \n"
                  << " sample of system configurations in a binarystream, write that to a \n"
                  << " new binarystream \n"
                  << " \n"
                  << " input expected in working directory: \n"
                  << "    configs-phi.infoheader    ... metadata about this timeseries \n"
                  << "    configs-phi.binarystream  ... phi configurations in binary doubles \n"
                  << " \n"
                  << " output written into working directory: \n"
                  << "    extracted-configs-phi.infoheader \n"
                  << "    extracted-configs-phi.binarystream \n"
                  << " \n"
                  << extractOptions << std::endl;
        return 0;
    }
    if (vm.count("version")) {
        std::cout << "Build info:\n"
                  << metadataToString(collectVersionInfo())
                  << std::endl;
        return 0;
    }
    
    // read meta information about configuration time series
    MetadataMap meta = readOnlyMetadata("configs-phi.infoheader");
    if (meta.empty()) {
        std::cerr << "Did not obtain meta information about configuration time series\n";
        return 1;
    }

    // process
    uint32_t sample_size = get_size_of_one_sample(meta);
    uint32_t extracted_samples = 0;
    extracted_samples = extract("configs-phi.binarystream",
                                "extracted-configs-phi.binarystream",
                                sample_size,
                                discard, subsample_interval);
    write_info("extracted-configs-phi.infoheader", meta,
               discard, subsample_interval, extracted_samples);    
    
    return 0;
}
