/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"
#pragma GCC diagnostic pop

#include "tools.h"
#include "detqmcparams.h"
#include "exceptions.h"

void DetQMCParams::check() {
    //check parameters
    using namespace boost::assign;
    std::vector<std::string> neededMCPars;
    neededMCPars += "sweeps", "thermalization", "jkBlocks", "timeseries", "measureInterval";
    for (auto p = neededMCPars.cbegin(); p != neededMCPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw_ParameterMissing(*p);
        }
    }
    if (specified.count("saveInterval") == 0) {
        saveInterval = sweeps;        //only save at end
    }

    if (sweeps % 2 != 0) {
        throw_ParameterWrong_message("Parameter sweeps must be even [else serialization consistency cannot be guaranteed]");
    }
    if (thermalization % 2 != 0) {
        throw_ParameterWrong_message("Parameter thermalization must be even [else serialization consistency cannot be guaranteed]");
    }
    if (saveInterval % 2 != 0) {
        throw_ParameterWrong_message("Parameter saveInterval must be even [else serialization consistency cannot be guaranteed]");
    }

    if (greenUpdateType_string == "simple") {
        greenUpdateType = GreenUpdateTypeSimple;
    } else if (greenUpdateType_string == "stabilized") {
        greenUpdateType = GreenUpdateTypeStabilized;
    } else {
        throw_ParameterWrong("greenUpdateType", greenUpdateType_string);
    }

    //some parameter consistency checking:
    if (sweeps % jkBlocks != 0) {
        throw_ParameterWrong_message("Number of jackknife blocks " + numToString(jkBlocks)
                                     + " does not match number of sweeps " + numToString(sweeps));
    }
    if ((measureInterval > sweeps) or
        (sweeps % measureInterval != 0)) {
        throw_ParameterWrong_message("Measurement interval " + numToString(measureInterval)
                                     + " ill-chosen for number of sweeps " + numToString(sweeps));
    }
    if (sweeps % saveInterval != 0) {
        throw_ParameterWrong_message("saveInterval (" + numToString(saveInterval) +
                                     ") needs to be a divisor of sweeps (" + numToString(sweeps) + ")");
    }

    if (specified.count("saveConfigurationStreamInterval") and
        (saveConfigurationStreamInterval % measureInterval != 0)) {
        throw_ParameterWrong_message("saveConfigurationStreamInterval must be an integer multiple of measureInterval");
    }

    if (not specified.count("saveConfigurationStreamInterval")) {
        saveConfigurationStreamInterval = measureInterval;
    }

}

MetadataMap DetQMCParams::prepareMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(VAR)
    META_INSERT(greenUpdateType_string);
    META_INSERT(simindex);
    META_INSERT(sweeps);
    META_INSERT(thermalization);
    META_INSERT(jkBlocks);
    META_INSERT(measureInterval);
    META_INSERT(saveInterval);
    META_INSERT(saveConfigurationStreamInterval);
    META_INSERT(rngSeed);
#undef META_INSERT
    meta["timeseries"] = (timeseries ? "true" : "false");
    meta["saveConfigurationStreamText"]   = (saveConfigurationStreamText   ? "true" : "false");
    meta["saveConfigurationStreamBinary"] = (saveConfigurationStreamBinary ? "true" : "false");    
    meta["stateFileName"] = stateFileName;
    return meta;
}
