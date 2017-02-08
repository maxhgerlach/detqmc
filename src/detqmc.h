/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * detqmc.h
 *
 * Handling of determinantal quantum Monte Carlo simulations
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */

#ifndef DETQMC_H_
#define DETQMC_H_

#include <vector>
#include <functional>
#include <memory>
#include <cstdlib>
#include <limits>
#include <ctime>
#include <functional>
#include <fstream>
#include <armadillo>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/preprocessor/comma.hpp"
#include "boost/timer/timer.hpp"
#include "boost/serialization/split_member.hpp"
#include "boost/assign/std/vector.hpp"
#include "boost/filesystem.hpp"
#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"
#pragma GCC diagnostic pop
#include "metadata.h"
#include "detqmcparams.h"
#include "detmodelloggingparams.h"
#include "detmodelparams.h"
#include "detmodel.h"
#include "observablehandler.h"
#include "rngwrapper.h"
#include "exceptions.h"
#include "tools.h"
#include "git-revision.h"
#include "timing.h"




// Class handling the simulation
template<class Model, class ModelParams = ModelParams<Model> >
class DetQMC {
public:
    //constructor to init a new simulation:
    DetQMC(const ModelParams& parsmodel, const DetQMCParams& parsmc,
           const DetModelLoggingParams& loggingParams = DetModelLoggingParams());

    //constructor to resume a simulation from a dumped state file:
    //we allow to change some MC parameters at this point:
    //  sweeps & saveInterval
    //if values > than the old values are specified, change them
    DetQMC(const std::string& stateFileName, const DetQMCParams& newParsmc);


    //carry out simulation determined by parsmc given in construction,
    //- handle thermalization & measurement stages as necessary
    //- save state and results periodically
    //- if granted walltime is almost over, save state & results
    //  and exit gracefully
    void run();

    // update results stored on disk
    void saveResults();
    // dump simulation parameters and the current state to a Boost::S11n archive,
    // also write out information about the current simulation state to info.dat
    void saveState();

    virtual ~DetQMC();
protected:
    //helper for constructors -- set all parameters and initialize contained objects
    void initFromParameters(const ModelParams& parsmodel, const DetQMCParams& parsmc,
                            const DetModelLoggingParams& loggingParams = DetModelLoggingParams());

    ModelParams parsmodel;
    DetQMCParams parsmc;
    DetModelLoggingParams parslogging;
    typedef DetQMCParams::GreenUpdateType GreenUpdateType;
    
    MetadataMap modelMeta;
    MetadataMap mcMeta;
    RngWrapper rng;
    std::unique_ptr<Model> replica;
    typedef std::unique_ptr<ScalarObservableHandler> ObsPtr;
    typedef std::unique_ptr<VectorObservableHandler> VecObsPtr;
    std::vector<ObsPtr> obsHandlers;
    std::vector<VecObsPtr> vecObsHandlers;      //need to be pointers: holds both KeyValueObservableHandlers and VectorObservableHandlers
    uint32_t sweepsDone;                        //Measurement sweeps done
    uint32_t sweepsDoneThermalization;          //thermalization sweeps done

    uint32_t swCounter; //helper counter in run() -- e.g. sweeps between measurements -- should also be serialized

    boost::timer::cpu_timer elapsedTimer;           //during this simulation run
    uint32_t curWalltimeSecs() {
        return static_cast<uint32_t>(elapsedTimer.elapsed().wall / 1000 / 1000 / 1000); // ns->mus->ms->s
    }
    uint32_t totalWalltimeSecs;             //this is serialized and carries the elapsed walltime in seconds
                                                    //accumulated over all runs, updated on call of saveResults()
    uint32_t walltimeSecsLastSaveResults;       //timer seconds at previous saveResults() call --> used to update totalWalltimeSecs
    uint32_t grantedWalltimeSecs;               //walltime the simulation is allowed to run
    std::string jobid;							//id string from the job scheduling system, or "nojobid"

private:
    //Serialize only the content data that has changed after construction.
    //Only call for deserialization after DetQMC has already been constructed and initialized!

    //separate functions loadContents, saveContents; both employ serializeContentsCommon
    template<class Archive>
    void loadContents(Archive& ar) {
        serializeContentsCommon(ar);

        replica->loadContents(ar);
    }

    template<class Archive>
    void saveContents(Archive& ar) {
        serializeContentsCommon(ar);

        replica->saveContents(ar);
    }


    template<class Archive>
    void serializeContentsCommon(Archive& ar) {
        ar & rng;                   //serialize completely

        for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
            //ATM no further derived classes of ScalarObservableHandler have a method serializeContents
            (*p)->serializeContents(ar);
        }
        for (auto p = vecObsHandlers.begin(); p != vecObsHandlers.end(); ++p) {
            //ATM no further derived classes of VectorObservableHandler have a method serializeContents
            (*p)->serializeContents(ar);
        }
        ar & sweepsDone & sweepsDoneThermalization;
        ar & swCounter;

        ar & totalWalltimeSecs;
    }
};








template<class Model, class ModelParams>
void DetQMC<Model,ModelParams>::initFromParameters(const ModelParams& parsmodel_, const DetQMCParams& parsmc_,
                                                   const DetModelLoggingParams& loggingParams /*default argument*/) {
    parsmodel = parsmodel_;
    parsmodel = updateTemperatureParameters(parsmodel);
    parsmc = parsmc_;
    parslogging = loggingParams;

    parsmc.check();
    parsmodel.check();
    parslogging.check();

    if (parsmc.specified.count("rngSeed") == 0) {
        std::cout << "No rng seed specified, will use std::time(0)" << std::endl;
        parsmc.rngSeed = (uint32_t) std::time(0);
    }
    rng = RngWrapper(parsmc.rngSeed, (parsmc.simindex + 1));

    createReplica(replica, rng, parsmodel, parslogging);    

    //prepare metadata
    modelMeta = parsmodel.prepareMetadataMap();
    mcMeta = parsmc.prepareMetadataMap();

    //prepare observable handlers
    auto scalarObs = replica->getScalarObservables();
    for (auto obsP = scalarObs.cbegin(); obsP != scalarObs.cend(); ++obsP) {
        obsHandlers.push_back(
            ObsPtr(new ScalarObservableHandler(*obsP, parsmc, modelMeta, mcMeta))
        );
    }
    auto vectorObs = replica->getVectorObservables();
    for (auto obsP = vectorObs.cbegin(); obsP != vectorObs.cend(); ++obsP) {
        vecObsHandlers.push_back(
            VecObsPtr(new VectorObservableHandler(*obsP, parsmc, modelMeta, mcMeta))
        );
    }
    auto keyValueObs = replica->getKeyValueObservables();
    for (auto obsP = keyValueObs.cbegin(); obsP != keyValueObs.cend(); ++obsP) {
        vecObsHandlers.push_back(
            VecObsPtr(new KeyValueObservableHandler(*obsP, parsmc, modelMeta, mcMeta))
        );
    }

    // setup files for system configuration streams [if files do not exist already]
    // [each process for its local replica]
    if (parsmc.saveConfigurationStreamText or parsmc.saveConfigurationStreamBinary) {
        std::string headerInfoText = metadataToString(modelMeta, "#")
            + metadataToString(mcMeta, "#");
        if (parsmc.saveConfigurationStreamText) {
            replica->saveConfigurationStreamTextHeader(headerInfoText);
        }
        if (parsmc.saveConfigurationStreamBinary) {
            replica->saveConfigurationStreamBinaryHeaderfile(headerInfoText);
        }
    }

    
    //query allowed walltime
    const char* pbs_walltime = std::getenv("PBS_WALLTIME");
    if (pbs_walltime) {
        grantedWalltimeSecs = fromString<decltype(grantedWalltimeSecs)>(pbs_walltime);
    } else {
        grantedWalltimeSecs = std::numeric_limits<decltype(grantedWalltimeSecs)>::max();
    }
    std::cout << "Granted walltime: " << grantedWalltimeSecs << " seconds.\n";

    //query SLURM Jobid
    const char* jobid_env = std::getenv("SLURM_JOBID");
    if (jobid_env) {
    	jobid = jobid_env;
    } else {
    	jobid = "nojobid";
    }
    std::cout << "Job ID: " << jobid << "\n";

    std::cout << "\nSimulation initialized, parameters: " << std::endl;
    std::cout << metadataToString(mcMeta, " ") << metadataToString(modelMeta, " ") << std::endl;
}



template<class Model, class ModelParams>
DetQMC<Model, ModelParams>::DetQMC(const ModelParams& parsmodel_, const DetQMCParams& parsmc_,
                                   const DetModelLoggingParams& parslogging_ /* default argument */) :
    parsmodel(), parsmc(), parslogging(),
    //proper initialization of default initialized members done in initFromParameters
    modelMeta(), mcMeta(), rng(), replica(),
    obsHandlers(), vecObsHandlers(),
    sweepsDone(0), sweepsDoneThermalization(),
    swCounter(0),
    elapsedTimer(),     // start timing
    totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
    grantedWalltimeSecs(0)
{
    initFromParameters(parsmodel_, parsmc_, parslogging_);
}

template<class Model, class ModelParams>
DetQMC<Model, ModelParams>::DetQMC(const std::string& stateFileName, const DetQMCParams& newParsmc) :
    parsmodel(), parsmc(), parslogging(),
    //proper initialization of default initialized members done by loading from archive
    modelMeta(), mcMeta(), rng(), replica(),
    obsHandlers(), vecObsHandlers(),
    sweepsDone(), sweepsDoneThermalization(),
    swCounter(0),
    elapsedTimer(),     // start timing
    totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
    grantedWalltimeSecs(0), jobid("")
{
    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit | std::ifstream::failbit);
    ifs.open(stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    DetModelLoggingParams parslogging_;
    ModelParams parsmodel_;
    DetQMCParams parsmc_;
    ia >> parslogging_ >> parsmodel_ >> parsmc_;

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
    SPECIFIED_INSERT_STR(stateFileName);
#undef SPECIFIED_INSERT_VAL
#undef SPECIFIED_INSERT_STR
    if (not parsmc_.greenUpdateType_string.empty()) {
        parsmc_.specified.insert("greenUpdateType");
    }
    
    initFromParameters(parsmodel_, parsmc_, parslogging_);
    loadContents(ia);

    std::cout << "\n"
              << "State of previous simulation has been loaded.\n"
              << "  sweepsDoneThermalization: " << sweepsDoneThermalization << "\n"
              << "  sweepsDone: " << sweepsDone << std::endl;
}

template<class Model, class ModelParams>
void DetQMC<Model, ModelParams>::saveState() {
    timing.start("saveState");

    //serialize state to file
    std::ofstream ofs;
    ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    ofs.open(parsmc.stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << parslogging << parsmodel << parsmc;
    saveContents(oa);

    //write out info about state of simulation to "info.dat"
    std::string commonInfoFilename = "info.dat";
    writeOnlyMetaData(commonInfoFilename, collectVersionInfo(),
                      "Collected information about this determinantal quantum Monte Carlo simulation",
                      false);
    MetadataMap modelMeta_current = replica->prepareModelMetadataMap();
    writeOnlyMetaData(commonInfoFilename, modelMeta_current,
                      "Model parameters and some data:",
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

    std::cout << "State has been saved." << std::endl;

    timing.stop("saveState");
}

template<class Model, class ModelParams>
DetQMC<Model, ModelParams>::~DetQMC() {
}


template<class Model, class ModelParams>
void DetQMC<Model, ModelParams>::run() {
    enum Stage { T, M, F };     //Thermalization, Measurement, Finished
    Stage stage = T;

    //local helper functions to initialize a "stage" of the big loop
    auto thermalizationStage = [&stage, this]() {
        stage = T;
        std::cout << "Thermalization for " << parsmc.thermalization << " sweeps..." << std::endl;
    };
    auto measurementsStage = [&stage, this]() {
        stage = M;
        std::cout << "Measurements for " << parsmc.sweeps << " sweeps..." << std::endl;
    };
    auto finishedStage = [&stage]() {
        stage = F;
        std::cout << "Measurements finished\n" << std::endl;
    };

    if (sweepsDoneThermalization < parsmc.thermalization) {
        thermalizationStage();
    } else if (sweepsDone < parsmc.sweeps) {
        measurementsStage();
    } else {
        finishedStage();
    }

    const uint32_t SavetyMinutes = 35;

    const std::string abortFilenames[] = { "ABORT." + jobid,
                                           "../ABORT." + jobid,
                                           "ABORT.all",
                                           "../ABORT.all" };

    while (stage != F) {                //big loop
    	if (swCounter % 2 == 0) {
            bool stop_now = false;
            if (curWalltimeSecs() > grantedWalltimeSecs - SavetyMinutes*60) {
                std::cout << "Granted walltime will be exceeded in less than " << SavetyMinutes << " minutes.\n";
                stop_now = true;
            } else {
                for (auto abortfn : abortFilenames) {
                    if (boost::filesystem::exists(abortfn)) {
                        std::cout << "Found file " << abortfn << ".\n";
                        stop_now = true;
                    }
                }
            }
            if (stop_now) {
                //close to exceeded walltime or we find that a file has been placed,
                //which signals us to abort this run for some other reason.
                //but only save state and exit if we have done an even
                //number of sweeps for ("economic") serialization guarantee [else do one sweep more]
                std::cout << "Current stage:\n"
                          << " sweeps done thermalization: " << sweepsDoneThermalization << "\n"
                          << " sweeps done measurements:   " << sweepsDone << "\n";
                std::cout << "Save state / results and exit gracefully... ";
                if (stage == Stage::M) {
                    saveResults();
                }
                saveState();
                std::cout << " OK " << std::endl;
                break;  //while
            }
    	}

        //thermalization & measurement stages
        switch (stage) {

        case T:
            switch(parsmc.greenUpdateType) {
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
                std::cout  << "  " << sweepsDoneThermalization << " ... saving state...";
                swCounter = 0;
                saveState();
                std::cout << " OK" << std::endl;
            }
            if (sweepsDoneThermalization == parsmc.thermalization) {
                std::cout << "Thermalization finished\n" << std::endl;
                replica->thermalizationOver();
                swCounter = 0;
                measurementsStage();
            }
            break;  //case

        case M: {
            ++swCounter;
            bool takeMeasurementNow = (swCounter % parsmc.measureInterval == 0);
            
            switch(parsmc.greenUpdateType) {
            case GreenUpdateType::GreenUpdateTypeSimple:
                replica->sweepSimple(takeMeasurementNow);
                break;
            case GreenUpdateType::GreenUpdateTypeStabilized:
                replica->sweep(takeMeasurementNow);
                break;
            }

            if (takeMeasurementNow) {
                for (auto ph = obsHandlers.begin(); ph != obsHandlers.end(); ++ph) {
                    (*ph)->insertValue(sweepsDone);
                }
                for (auto ph = vecObsHandlers.begin(); ph != vecObsHandlers.end(); ++ph) {
                    (*ph)->insertValue(sweepsDone);
                }

                if (swCounter % parsmc.saveConfigurationStreamInterval == 0) {
                    // This is a good time to write the current system configuration to disk
                    if (parsmc.saveConfigurationStreamText) {
                        replica->saveConfigurationStreamText();
                    }
                    if (parsmc.saveConfigurationStreamBinary) {
                        replica->saveConfigurationStreamBinary();
                    }
                }
            }
            ++sweepsDone;
            if (swCounter == parsmc.saveInterval) {
                std::cout << "  " << sweepsDone << " ... saving results and state ...";
                swCounter = 0;
                saveResults();
                saveState();
                std::cout << " OK" << std::endl;
            }
            if (sweepsDone == parsmc.sweeps) {
                swCounter = 0;
                finishedStage();
            }
            break;  //case
        }

        case F:
            break;  //case

        }  //switch
    }
}






template<class Model, class ModelParams>
void DetQMC<Model, ModelParams>::saveResults() {
    timing.start("saveResults");

    outputResults(obsHandlers);
    for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
        (*p)->outputTimeseries();
    }
    outputResults(vecObsHandlers);

    timing.stop("saveResults");
}




#endif /* DETQMC_H_ */
