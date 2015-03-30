#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/iostreams/stream.hpp"
#include "boost/iostreams/device/back_inserter.hpp"
#include "boost/iostreams/device/array.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#pragma GCC diagnostic pop
#include "boost_serialize_armadillo.h"
#include "detsdwsystemconfig.h"
#include "exceptions.h"

DetSDW_SystemConfig::DetSDW_SystemConfig()
    : L(0), m(0), opdim(0), phi(), cdwl()
{
}

DetSDW_SystemConfig::DetSDW_SystemConfig(const ModelParamsDetSDW& pars, const CubeNum& phi_current, const MatInt& cdwl_current)
    : L(pars.L), m(pars.m), opdim(pars.opdim), phi(phi_current), cdwl(cdwl_current)
{
}

DetSDW_SystemConfig::DetSDW_SystemConfig(const ModelParamsDetSDW& pars, const CubeNum& phi_current)
    : L(pars.L), m(pars.m), opdim(pars.opdim), phi(phi_current), cdwl()
{
}

void DetSDW_SystemConfig::write_to_disk(DetSDW_SystemConfig_FileHandle& file_handle) const {
    if (not phi.is_empty()) {
        // we have a phi config to save
        bool could_save_phi = false;
        if (file_handle.phi_output_text) {
            // the pointer is not NULL
            if (*(file_handle.phi_output_text)) {
                // the ofstream is valid
                write_to_disk_phi_text(file_handle);
                could_save_phi = true;
            } else {
                throw_GeneralError("std::ofstream *file_handle.phi_output_text is not ready");
            }
        }
        if (file_handle.phi_output_binary) {
            // the pointer is not NULL
            if (*(file_handle.phi_output_binary)) {
                write_to_disk_phi_binary(file_handle);
                could_save_phi = true;
            } else {
                throw_GeneralError("std::ofstream *file_handle.phi_output_binary is not ready");
            }
        }
        if (not could_save_phi) {
            throw_GeneralError("file_handle does not have a valid handle to save phi system configuration");
        }
    } else {
        throw_GeneralError("No phi system configuration available to save");
    }

    if (not cdwl.is_empty()) {
        // we have a cdwl config to save
        bool could_save_cdwl = false;
        if (file_handle.cdwl_output_text) {
            // the pointer is not NULL
            if (*(file_handle.cdwl_output_text)) {
                // the ofstream is valid
                write_to_disk_cdwl_text(file_handle);
                could_save_cdwl = true;
            } else {
                throw_GeneralError("std::ofstream *file_handle.cdwl_output_text is not ready");
            }
        }
        if (file_handle.cdwl_output_binary) {
            // the pointer is not NULL
            if (*(file_handle.cdwl_output_binary)) {
                write_to_disk_cdwl_binary(file_handle);
                could_save_cdwl = true;
            } else {
                throw_GeneralError("std::ofstream *file_handle.cdwl_output_binary is not ready");
            }
        }
        if (not could_save_cdwl) {
            throw_GeneralError("file_handle does not have a valid handle to save cdwl system configuration");
        }
        
    }
}

void DetSDW_SystemConfig::write_to_disk_phi_text(DetSDW_SystemConfig_FileHandle& file_handle) const {
    file_handle.phi_output_text->precision(14);
    file_handle.phi_output_text->setf(std::ios::scientific, std::ios::floatfield);

    for (uint32_t ix = 0; ix < L; ++ix) {
        for (uint32_t iy = 0; iy < L; ++iy) {
            uint32_t i = iy*L + ix;
            for (uint32_t k = 1; k <= m; ++k) {
                for (uint32_t dim = 0; dim < opdim; ++dim) {
                    *(file_handle.phi_output_text) << phi(i, dim, k) << "\n";
                }
            }
        }
    }
}

void DetSDW_SystemConfig::write_to_disk_cdwl_text(DetSDW_SystemConfig_FileHandle& file_handle) const {
    for (uint32_t ix = 0; ix < L; ++ix) {
        for (uint32_t iy = 0; iy < L; ++iy) {
            uint32_t i = iy*L + ix;
            for (uint32_t k = 1; k <= m; ++k) {
                *(file_handle.cdwl_output_text) << cdwl(i, k) << "\n";
            }
        }
    }    
}

void DetSDW_SystemConfig::write_to_disk_phi_binary(DetSDW_SystemConfig_FileHandle& file_handle) const {
    for (uint32_t ix = 0; ix < L; ++ix) {
        for (uint32_t iy = 0; iy < L; ++iy) {
            uint32_t i = iy*L + ix;
            for (uint32_t k = 1; k <= m; ++k) {
                for (uint32_t dim = 0; dim < opdim; ++dim) {
                    file_handle.phi_output_binary->write(reinterpret_cast<const char*>(&(phi(i, dim, k))),
                                                         sizeof(phi(i, dim, k)));
                }
            }
        }
    }
}

void DetSDW_SystemConfig::write_to_disk_cdwl_binary(DetSDW_SystemConfig_FileHandle& file_handle) const {
    for (uint32_t ix = 0; ix < L; ++ix) {
        for (uint32_t iy = 0; iy < L; ++iy) {
            uint32_t i = iy*L + ix;
            for (uint32_t k = 1; k <= m; ++k) {
                file_handle.cdwl_output_binary->write(reinterpret_cast<const char*>(&(cdwl(i, k))),
                                                      sizeof(cdwl(i, k)));
            }
        }
    }
}




void serialize_systemConfig_to_buffer(std::string& buffer, const DetSDW_SystemConfig& systemConfig) {
    //serialize object into a std::basic_string
    //cf. http://stackoverflow.com/questions/3015582/direct-boost-serialization-to-char-array
    namespace ios = boost::iostreams;
    ios::back_insert_device<std::string> inserter(buffer);
    ios::stream<ios::back_insert_device<std::string> > s(inserter);
    boost::archive::binary_oarchive oa(s);

    oa << systemConfig;

    s.flush();    
}

void deserialize_systemConfig_from_buffer(DetSDW_SystemConfig& systemConfig, const std::string& buffer) {
    //wrap buffer inside a stream and deserialize into object
    //cf. http://stackoverflow.com/questions/3015582/direct-boost-serialization-to-char-array
    namespace ios = boost::iostreams;
    ios::basic_array_source<char> device(buffer.data(), buffer.size());
    ios::stream<ios::basic_array_source<char> > s(device);
    boost::archive::binary_iarchive ia(s);

    ia >> systemConfig;    
}
