#ifndef DETSDWSYSTEMCONFIG_H_
#define DETSDWSYSTEMCONFIG_H_

#include <armadillo>
#include "detsdwsystemconfigfilehandle.h"
#include "detsdwparams.h"

typedef arma::Col<num> VecNum;
typedef arma::Cube<num> CubeNum;
typedef arma::Mat<int32_t> MatInt;



// class to pass around and save system configurations for the SDW model

class DetSDW_SystemConfig {
public:
    DetSDW_SystemConfig(const ModelParamsDetSDW& pars, const CubeNum& phi_current, const MatInt& cdwl_current);
    DetSDW_SystemConfig(const ModelParamsDetSDW& pars, const CubeNum& phi_current);

    void write_to_disk(DetSDW_SystemConfig_FileHandle& file_handle) const;
private:
    uint32_t L;
    uint32_t m;
    uint32_t opdim;
    //slice indexes timeslice, column indexes order parameter
    // dimension, row indexes site    
    // [Index order: row, col, slice]
    // [Data layout; slice after slice, matrices then column major] 
    CubeNum phi;

    //discrete field for cdwU: l_i(\tau_k) == cdwl(i,k).
    //row indexes site, column indexes timeslice
    MatInt cdwl;

    void write_to_disk_phi_text(DetSDW_SystemConfig_FileHandle& file_handle) const;
    void write_to_disk_cdwl_text(DetSDW_SystemConfig_FileHandle& file_handle) const;
    void write_to_disk_phi_binary(DetSDW_SystemConfig_FileHandle& file_handle) const;
    void write_to_disk_cdwl_binary(DetSDW_SystemConfig_FileHandle& file_handle) const;

private:
    friend class boost::serialization::access;
    
    // serialization to pass around system configurations by MPI
    template<class Archive>
    void serialize(Archive & ar, const unsigned int /*version*/) {
        ar & L & m & opdim
           & phi & cdwl;
    }
};


void serialize_systemConfig_to_buffer(std::string& buffer, const DetSDW_SystemConfig& systemConfig);

void deserialize_systemConfig_from_buffer(DetSDW_SystemConfig& systemConfig, const std::string& buffer);


#endif//DETSDWSYSTEMCONFIG_H_
