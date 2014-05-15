#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"
#pragma GCC diagnostic pop

#include "tools.h"
#include "parameters.h"

void DetQMCParams::check() {
    //check parameters
    if (parsmodel.specified.count("model") == 0) {
        throw ParameterMissing("model");
    }
    using namespace boost::assign;
    std::vector<std::string> neededMCPars;
    neededMCPars += "sweeps", "thermalization", "jkBlocks", "timeseries", "measureInterval";
    for (auto p = neededMCPars.cbegin(); p != neededMCPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }
    if (specified.count("saveInterval") == 0) {
        saveInterval = sweeps;        //only save at end
    }

    if (sweeps % 2 != 0) {
        throw ParameterWrong("Parameter sweeps must be even [else serialization consistency cannot be guaranteed]");
    }
    if (thermalization % 2 != 0) {
        throw ParameterWrong("Parameter thermalization must be even [else serialization consistency cannot be guaranteed]");
    }
    if (saveInterval % 2 != 0) {
        throw ParameterWrong("Parameter saveInterval must be even [else serialization consistency cannot be guaranteed]");
    }

    if (greenUpdateType_string == "simple") {
        greenUpdateType = GreenUpdateTypeSimple;
    } else if (greenUpdateType_string == "stabilized") {
        greenUpdateType = GreenUpdateTypeStabilized;
    } else {
        throw ParameterWrong("greenUpdateType", greenUpdateType_string);
    }

    //some parameter consistency checking:
    if (parsmc.sweeps % parsmc.jkBlocks != 0) {
        throw ParameterWrong("Number of jackknife blocks " + numToString(parsmc.jkBlocks)
                             + " does not match number of sweeps " + numToString(parsmc.sweeps));
    }
    if ((parsmc.measureInterval > parsmc.sweeps) or
        (parsmc.sweeps % parsmc.measureInterval != 0)) {
        throw ParameterWrong("Measurement interval " + numToString(parsmc.measureInterval)
                             + " ill-chosen for number of sweeps " + numToString(parsmc.sweeps));
    }
    if (parsmc.sweeps % parsmc.saveInterval != 0) {
        throw ParameterWrong("saveInterval (" + numToString(parsmc.saveInterval) +
                             ") needs to be a divisor of sweeps (" + numToString(parsmc.sweeps) + ")");
    }

}

MetadataMap prepareMetaDataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(VAR)
    META_INSERT(greenUpdateType_string);
    META_INSERT(sweeps);
    META_INSERT(thermalization);
    META_INSERT(jkBlocks);
    META_INSERT(measureInterval);
    META_INSERT(saveInterval);
    META_INSERT(rngSeed);
#undef META_INSERT
    meta["timeseries"] = (timeseries ? "true" : "false");
    meta["stateFileName"] = stateFileName;
    return meta;
}
