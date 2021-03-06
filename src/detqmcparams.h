/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

#ifndef DETQMCPARAMS_H_
#define DETQMCPARAMS_H_

#include <string>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"
#pragma GCC diagnostic pop

#include "metadata.h"

// Collect various structs defining various parameters. [more structs
// in the other ...params.h headers]

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
    uint32_t simindex;          // an index to discriminate between multiple simulation instances for the same set of parameters
    uint32_t sweeps;            // number of sweeps used for measurements
    uint32_t thermalization;    // number of warm-up sweeps allowed before equilibrium is assumed
    uint32_t jkBlocks;          // number of jackknife blocks for error estimation
    bool timeseries;            // if true, write time series of individual measurements to disk
    uint32_t measureInterval;   // take measurements every measureInterval sweeps
    uint32_t saveInterval;      // write measurements to disk every saveInterval sweeps
    uint32_t saveConfigurationStreamInterval; // interval in sweeps where full system configurations are buffered and saved to disk.  Must be an integer multiple of measureInterval.  This is only effective if one of the boolean flags for saving configurations is set to true.
    uint32_t rngSeed;           // seed for random number generator

    std::string greenUpdateType_string;    //"simple" or "stabilized"
    enum GreenUpdateType {GreenUpdateTypeSimple, GreenUpdateTypeStabilized};
    GreenUpdateType greenUpdateType;    

    bool saveConfigurationStreamText;
    bool saveConfigurationStreamBinary;
    
    std::string stateFileName;      //for serialization dumps
    bool sweepsHasChanged;          //true, if the number of target sweeps has changed after resuming

    std::set<std::string> specified; // used to record names of specified parameters

    DetQMCParams() :
        simindex(0), sweeps(), thermalization(), jkBlocks(), timeseries(false), measureInterval(), saveInterval(),
        saveConfigurationStreamInterval(0),
        rngSeed(), greenUpdateType_string(), saveConfigurationStreamText(false), saveConfigurationStreamBinary(false),
        stateFileName(), sweepsHasChanged(false), specified()
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
        ar  & simindex
            & sweeps & thermalization & jkBlocks & timeseries
            & measureInterval & saveInterval & saveConfigurationStreamInterval
            & rngSeed
            & greenUpdateType_string & greenUpdateType
            & saveConfigurationStreamText & saveConfigurationStreamBinary
            & stateFileName
            & sweepsHasChanged
            & specified;
    }
};



#endif /* DETQMCPARAMS_H_ */
