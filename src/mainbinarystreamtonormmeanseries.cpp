#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include "metadata.h"
#include "tools.h"

// convert a sdw model binary stream to a normMeanPhi time series in stdout


// return size of sample in units of doubles
uint32_t get_size_of_one_sample(const MetadataMap& meta) {
    uint32_t N = fromString<uint32_t>(meta.at("N"));
    uint32_t m = fromString<uint32_t>(meta.at("m"));
    uint32_t opdim = fromString<uint32_t>(meta.at("opdim"));
    return N * m * opdim;
}


int main(int argc, char *argv[])
{
    MetadataMap meta = readOnlyMetadata("configs-phi.infoheader");
    uint32_t opdim = fromString<uint32_t>(meta.at("opdim"));
    uint32_t sample_size = get_size_of_one_sample(meta);

    std::ifstream binary_float_input("configs-phi.binarystream", std::ios::in | std::ios::binary);

    std::cout.precision(14);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);

    std::vector<double> one_sample;
    one_sample.resize(sample_size, 0.0);
    
    while (binary_float_input) {
        binary_float_input.read(reinterpret_cast<char*>(&(*one_sample.begin())), sizeof(double) * sample_size);
        if (binary_float_input) {
            //no failure

            double squared_norm = 0.0;
            for (uint32_t offset = 0; offset < opdim; ++offset) {
                double component_sum = 0;
                for (uint32_t index = offset; index < sample_size; index += opdim) {
                    component_sum += one_sample[index];
                }
                squared_norm += std::pow((component_sum / (sample_size / opdim)), 2);
            }
            double norm = std::sqrt(squared_norm);

            std::cout << norm << "\n";
        }
    }
    
    return 0;
}
