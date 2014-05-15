/*
 * parameters.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"
#pragma GCC diagnostic pop

#include "metadata.h"

// Collect various structs defining various parameters.

// The set specified included in each struct contains string representations
// of all parameters actually specified.  This allows throwing an exception
// at the appropriate point in the program if a parameter is missing.


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
typedef double num;     //possibility to switch to single precision if ever desired
#pragma GCC diagnostic pop

// Representation of Monte Carlo simulation parameters for DetQMC
struct DetQMCParams {
    //public data members, call check() after setting these
    uint32_t sweeps;            // number of sweeps used for measurements
    uint32_t thermalization;    // number of warm-up sweeps allowed before equilibrium is assumed
    uint32_t jkBlocks;          // number of jackknife blocks for error estimation
    bool timeseries;            // if true, write time series of individual measurements to disk
    uint32_t measureInterval;   // take measurements every measureInterval sweeps
    uint32_t saveInterval;      // write measurements to disk every saveInterval sweeps
    uint32_t rngSeed;           // seed for random number generator

    std::string greenUpdateType_string;    //"simple" or "stabilized"
    enum GreenUpdateType {GreenUpdateTypeSimple, GreenUpdateTypeStabilized};
    GreenUpdateType greenUpdateType;    

    std::string stateFileName;      //for serialization dumps
    bool sweepsHasChanged;          //true, if the number of target sweeps has changed after resuming

    std::set<std::string> specified; // used to record names of specified parameters

    DetQMCParams() :
        sweeps(), thermalization(), jkBlocks(), timeseries(false), measureInterval(), saveInterval(),
        rngSeed(), greenUpdateType_string(), stateFileName(), sweepsHasChanged(false), specified()
    { }

    // check consistency, convert strings to enums
    void check();

    // to export human readable form of parameters
    MetadataMap prepareMetadataMap() const;
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar  & sweeps & thermalization & jkBlocks & timeseries
            & measureInterval & saveInterval & rngSeed
            & greenUpdateType_string & greenUpdateType
            & stateFileName
            & sweepsHasChanged
            & specified;
    }
};



#endif /* PARAMETERS_H_ */
