/*
 * detqmc.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include <cstdlib>
#include <limits>
#include <ctime>
#include <functional>
#include <fstream>
#include <armadillo>
#include "boost/assign/std/vector.hpp"

#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"

#include "tools.h"
#include "detqmc.h"
#include "dethubbard.h"
#include "detsdw.h"
#include "tools.h"
#include "git-revision.h"
#include "exceptions.h"
#include "timing.h"

using std::cout;
using std::endl;

void DetQMC::initFromParameters(const ModelParams& parsmodel_, const MCParams& parsmc_) {
    parsmodel = parsmodel_;
    parsmc = parsmc_;
    //check parameters
    if (parsmodel.specified.count("model") == 0) {
        throw ParameterMissing("model");
    }
    using namespace boost::assign;
    std::vector<std::string> neededMCPars;
    neededMCPars += "sweeps", "thermalization", "jkBlocks", "timeseries", "measureInterval";
    for (auto p = neededMCPars.cbegin(); p != neededMCPars.cend(); ++p) {
        if (parsmc.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }
    if (parsmc.specified.count("saveInterval") == 0) {
        parsmc.saveInterval = parsmc.sweeps;        //only save at end
    }

    if (parsmc.sweeps % 2 != 0) {
        throw ParameterWrong("Parameter sweeps must be even [else serialization consistency cannot be guaranteed]");
    }
    if (parsmc.thermalization % 2 != 0) {
        throw ParameterWrong("Parameter thermalization must be even [else serialization consistency cannot be guaranteed]");
    }
    if (parsmc.saveInterval % 2 != 0) {
        throw ParameterWrong("Parameter saveInterval must be even [else serialization consistency cannot be guaranteed]");
    }

    if (parsmc.specified.count("rngSeed") == 0) {
        cout << "No rng seed specified, will use std::time(0)" << endl;
        parsmc.rngSeed = (uint32_t) std::time(0);
    }
    rng = RngWrapper(parsmc.rngSeed);

    if (parsmodel.model == "hubbard") {
        replica = createDetHubbard(rng, parsmodel);
    } else if (parsmodel.model == "sdw") {
        replica = createDetSDW(rng, parsmodel);
    } else {
        throw ParameterWrong("model", parsmodel.model);
    }

    if (parsmc.greenUpdateType == "simple") {
        greenUpdateType = GreenUpdateTypeSimple;
    } else if (parsmc.greenUpdateType == "stabilized") {
        greenUpdateType = GreenUpdateTypeStabilized;
    } else {
        throw ParameterWrong("greenUpdateType", parsmc.greenUpdateType);
    }

//  if (greenUpdateType == GreenUpdateType::Simple) {
//      sweepFunc = [this]() {replica->sweepSimple();};
//      sweepThermalizationFunc= [this]() {replica->sweepSimpleThermalization();};
//  } else if (greenUpdateType == GreenUpdateType::Stabilized) {
//      sweepFunc = [this]() {replica->sweep();};
//      sweepThermalizationFunc = [this]() {replica->sweepThermalization();};
//  } else {
//      throw GeneralError("greenUpdateType not defined!");
//  }

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

    //prepare metadata
    modelMeta = replica->prepareModelMetadataMap();
    mcMeta = prepareMCMetadataMap();

    //prepare observable handlers
    auto scalarObs = replica->getScalarObservables();
    for (auto obsP = scalarObs.cbegin(); obsP != scalarObs.cend(); ++obsP) {
        obsHandlers.push_back(ObsPtr(
                new ScalarObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
    }
    auto vectorObs = replica->getVectorObservables();
    for (auto obsP = vectorObs.cbegin(); obsP != vectorObs.cend(); ++obsP) {
        vecObsHandlers.push_back(VecObsPtr(
                new VectorObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
    }
    auto keyValueObs = replica->getKeyValueObservables();
    for (auto obsP = keyValueObs.cbegin(); obsP != keyValueObs.cend(); ++obsP) {
        vecObsHandlers.push_back(VecObsPtr(
                new KeyValueObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
    }

    //query allowed walltime
    const char* pbs_walltime = std::getenv("PBS_WALLTIME");
    if (pbs_walltime) {
        grantedWalltimeSecs = fromString<decltype(grantedWalltimeSecs)>(pbs_walltime);
    } else {
        grantedWalltimeSecs = std::numeric_limits<decltype(grantedWalltimeSecs)>::max();
    }
    cout << "Granted walltime: " << grantedWalltimeSecs << " seconds.\n";

    cout << "\nSimulation initialized, parameters: " << endl;
    cout << metadataToString(mcMeta, " ") << metadataToString(modelMeta, " ") << endl;
}

DetQMC::DetQMC(const ModelParams& parsmodel_, const MCParams& parsmc_) :
        parsmodel(), parsmc(),
        //proper initialization of default initialized members done in initFromParameters
        greenUpdateType(), //sweepFunc(), sweepThermalizationFunc(),
        modelMeta(), mcMeta(), rng(), replica(),
        obsHandlers(), vecObsHandlers(),
        sweepsDone(0), sweepsDoneThermalization(),
        swCounter(0),
        elapsedTimer(),     // start timing
        totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
        grantedWalltimeSecs(0)
{
    initFromParameters(parsmodel_, parsmc_);
}

DetQMC::DetQMC(const std::string& stateFileName, const MCParams& newParsmc) :
        parsmodel(), parsmc(),
        //proper initialization of default initialized members done by loading from archive
        greenUpdateType(), //sweepFunc(), sweepThermalizationFunc(),
        modelMeta(), mcMeta(), rng(), replica(),
        obsHandlers(), vecObsHandlers(),
        sweepsDone(), sweepsDoneThermalization(),
        swCounter(0),
        elapsedTimer(),     // start timing
        totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
        grantedWalltimeSecs(0)
{
    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit | std::ifstream::failbit);
    ifs.open(stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ModelParams parsmodel_;
    MCParams parsmc_;
    ia >> parsmodel_ >> parsmc_;

    if (newParsmc.sweeps > parsmc_.sweeps) {
        std::cout << "Target sweeps will be changed from " << parsmc_.sweeps
                  << " to " << newParsmc.sweeps << std::endl;
        parsmc_.sweeps = newParsmc.sweeps;
        parsmc_.sweepsHasChanged = true;
    }
    if (newParsmc.saveInterval > 0 and newParsmc.saveInterval != parsmc_.saveInterval) {
        std::cout << "saveInterval will be changed from " << parsmc_.saveInterval
                  << " to " << newParsmc.saveInterval << std::endl;
        parsmc_.saveInterval = newParsmc.saveInterval;
    }
    parsmc_.stateFileName = stateFileName;

    //make sure mcparams are set correctly as "specified"
#define SPECIFIED_INSERT_VAL(x) if (parsmc_.x) { parsmc_.specified.insert(#x); }
#define SPECIFIED_INSERT_STR(x) if (not parsmc_.x.empty()) { parsmc_.specified.insert(#x); }
    SPECIFIED_INSERT_VAL(sweeps);
    SPECIFIED_INSERT_VAL(thermalization);
    SPECIFIED_INSERT_VAL(jkBlocks);
    SPECIFIED_INSERT_VAL(measureInterval);
    SPECIFIED_INSERT_VAL(saveInterval);
    SPECIFIED_INSERT_STR(greenUpdateType);
    SPECIFIED_INSERT_STR(stateFileName);
#undef SPECIFIED_INSERT_VAL
#undef SPECIFIED_INSERT_STR

    initFromParameters(parsmodel_, parsmc_);
    loadContents(ia);

    std::cout << "\n"
    		  << "State of previous simulation has been loaded.\n"
    		  << "  sweepsDoneThermalization: " << sweepsDoneThermalization << "\n"
    		  << "  sweepsDone: " << sweepsDone << std::endl;
}

void DetQMC::saveState() {
    timing.start("saveState");
    std::ofstream ofs;
    ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    ofs.open(parsmc.stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << parsmodel << parsmc;
    saveContents(oa);
    timing.stop("saveState");
}

DetQMC::~DetQMC() {
}


void DetQMC::run() {
    enum Stage { T, M, F };     //Thermalization, Measurement, Finished
    Stage stage = T;

    //local helper functions to initialize a "stage" of the big loop
    auto thermalizationStage = [&stage, this]() {
        stage = T;
        cout << "Thermalization for " << parsmc.thermalization << " sweeps..." << endl;
    };
    auto measurementsStage = [&stage, this]() {
        stage = M;
        cout << "Measurements for " << parsmc.sweeps << " sweeps..." << endl;
    };
    auto finishedStage = [&stage]() {
        stage = F;
        cout << "Measurements finished\n" << endl;
    };

    if (sweepsDoneThermalization < parsmc.thermalization) {
        thermalizationStage();
    } else if (sweepsDone < parsmc.sweeps) {
        measurementsStage();
    } else {
        finishedStage();
    }

    const uint32_t SavetyMinutes = 35;

    while (stage != F) {                //big loop
        if (curWalltimeSecs() > grantedWalltimeSecs - SavetyMinutes*60) {
            //close to exceeded walltime, but only save state and exit if we have done an even
            //number of sweeps for ("economic") serialization guarantee [else do one sweep more]
            if (swCounter % 2 == 0) {
                cout << "Granted walltime will be exceeded in less than " << SavetyMinutes << " minutes.\n"
                        << "Save state / results and exit gracefully."
                        << endl;
                if (stage == Stage::M) {
                    saveResults();
                }
                saveState();
                break;  //while
            }
        }

        //thermalization & measurement stages
        switch (stage) {

        case T:
            switch(greenUpdateType) {
            case GreenUpdateType::GreenUpdateTypeSimple:
                replica->sweepSimpleThermalization();
                break;
            case GreenUpdateType::GreenUpdateTypeStabilized:
                replica->sweepThermalization();
                break;
            }
            ++sweepsDoneThermalization;
            ++swCounter;
            if (swCounter == parsmc.saveInterval) {
                cout  << "  " << sweepsDoneThermalization << " ... saving state...";
                swCounter = 0;
                saveState();
                cout << endl;
            }
            if (sweepsDoneThermalization == parsmc.thermalization) {
                cout << "Thermalization finished\n" << endl;
                replica->thermalizationOver();
                swCounter = 0;
                measurementsStage();
            }
            break;  //case

        case M:
            switch(greenUpdateType) {
            case GreenUpdateTypeSimple:
                replica->sweepSimple();
                break;
            case GreenUpdateTypeStabilized:
                replica->sweep();
                break;
            }
            ++swCounter;
            if (swCounter % parsmc.measureInterval == 0) {
                replica->measure();
                for (auto ph = obsHandlers.begin(); ph != obsHandlers.end(); ++ph) {
                    (*ph)->insertValue(sweepsDone);
                }
                for (auto ph = vecObsHandlers.begin(); ph != vecObsHandlers.end(); ++ph) {
                    (*ph)->insertValue(sweepsDone);
                }
            }
            ++sweepsDone;
            if (swCounter == parsmc.saveInterval) {
                cout << "  " << sweepsDone << " ... saving results and state ...";
                swCounter = 0;
                saveResults();
                saveState();
                cout << endl;
            }
            if (sweepsDone == parsmc.sweeps) {
                swCounter = 0;
                finishedStage();
            }
            break;  //case

        case F:
            break;  //case

        }
    }
}


MetadataMap DetQMC::prepareMCMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(parsmc.VAR)
    META_INSERT(greenUpdateType);
    META_INSERT(sweeps);
    META_INSERT(thermalization);
    META_INSERT(jkBlocks);
    META_INSERT(measureInterval);
    META_INSERT(saveInterval);
    META_INSERT(rngSeed);
#undef META_INSERT
    meta["timeseries"] = (parsmc.timeseries ? "true" : "false");
    return meta;
}


void DetQMC::saveResults() {
    timing.start("saveResults");

    outputResults(obsHandlers);
    for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
        (*p)->outputTimeseries();
    }
    outputResults(vecObsHandlers);
    std::string commonInfoFilename = "info.dat";
    writeOnlyMetaData(commonInfoFilename, collectVersionInfo(),
            "Collected innformation about this determinantal quantum Monte Carlo simulation",
            false);
    writeOnlyMetaData(commonInfoFilename, modelMeta,
            "Model parameters:",
            true);
    writeOnlyMetaData(commonInfoFilename, mcMeta,
            "Monte Carlo parameters:",
            true);

    MetadataMap currentState;
    currentState["sweepsDoneThermalization"] = numToString(sweepsDoneThermalization);
    currentState["sweepsDone"] = numToString(sweepsDone);

    uint32_t cwts = curWalltimeSecs();
    totalWalltimeSecs += (cwts - walltimeSecsLastSaveResults);
    walltimeSecsLastSaveResults = cwts;

    currentState["totalWallTimeSecs"] = numToString(totalWalltimeSecs);
    writeOnlyMetaData(commonInfoFilename, currentState,
            "Current state of simulation:",
            true);

    timing.stop("saveResults");
}
