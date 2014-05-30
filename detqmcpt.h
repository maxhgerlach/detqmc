// Parallel Tempering Determinantal QMC simulation handling

#ifndef DETQMC_H_
#define DETQMC_H_

#include <mpi.h>
#include <vector>
#include <functional>
#include <numeric>              // std::iota
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
#include "detmodelparams.h"
#include "detmodel.h"
#include "observablehandler.h"
#include "rngwrapper.h"
#include "exceptions.h"
#include "tools.h"
#include "git-revision.h"
#include "timing.h"



class SerializeContentsKey;

// Class handling the simulation
template<class Model, class ModelParams = ModelParams<Model> >
class DetQMCPT {
public:
    //constructor to init a new simulation:
    DetQMCPT(const ModelParams& parsmodel, const DetQMCParams& parsmc,
             const DetQMCPTParams& parspt);

    //constructor to resume a simulation from a dumped state file:
    //we allow to change some MC parameters at this point:
    //  sweeps & saveInterval
    //if values > than the old values are specified, change them
    DetQMCPT(const std::string& stateFileName, const DetQMCParams& newParsmc,
             const DetQMCPTParams& newParspt);


    //carry out simulation determined by parsmc and parspt given in construction,
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

    virtual ~DetQMCPT();
protected:
    //helper for constructors -- set all parameters and initialize contained objects
    void initFromParameters(const ModelParams& parsmodel, const DetQMCParams& parsmc,
                            const DetQMCPTParams& parspt);

    ModelParams parsmodel;
    DetQMCParams parsmc;
    DetQMCPTParams parspt;
    typedef DetQMCParams::GreenUpdateType GreenUpdateType;
    
    MetadataMap modelMeta;
    MetadataMap mcMeta;
    MetadataMap ptMeta;    
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

    //MPI specifics:
    uint32_t numProcesses;      //total number of parallel processes
    uint32_t processIndex;      //MPI-rank of the current process

    // specific to the root process
    std::vector<uint32_t> current_process_par; // indexed by process rank number, giving control parameter index
                                               // currently associated to the replica at that process
    std::vector<uint32_t> current_par_process; // the reverse association
    std::vector<double> exchange_action;          // for each replica: its locally measured exchange action
private:
    //Serialize only the content data that has changed after construction.
    //Only call for deserialization after DetQMCPT has already been constructed and initialized!

    //separate functions loadContents, saveContents; both employ serializeContentsCommon
    template<class Archive>
    void loadContents(Archive& ar) {
        serializeContentsCommon(ar);

        replica->loadContents(SerializeContentsKey(), ar);
    }

    template<class Archive>
    void saveContents(Archive& ar) {
        serializeContentsCommon(ar);

        replica->saveContents(SerializeContentsKey(), ar);
    }


    template<class Archive>
    void serializeContentsCommon(Archive& ar) {
        ar & rng;                   //serialize completely

        for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
            //ATM no further derived classes of ScalarObservableHandler have a method serializeContents
            (*p)->serializeContents(SerializeContentsKey(), ar);
        }
        for (auto p = vecObsHandlers.begin(); p != vecObsHandlers.end(); ++p) {
            //ATM no further derived classes of VectorObservableHandler have a method serializeContents
            (*p)->serializeContents(SerializeContentsKey(), ar);
        }
        ar & sweepsDone & sweepsDoneThermalization;
        ar & swCounter;

        ar & totalWalltimeSecs;

        ar & current_process_par;

        ar & current_par_process;
    }
};



//Only few member functions of DetQMCPT are allowed to make instances of
//this class.  In this way access to the member functions
//serializeContents(), saveContents(), loadContents() of other classes
//is restricted.  Compare to
//http://stackoverflow.com/questions/6310720/declare-a-member-function-of-a-forward-declared-class-as-friend
class SerializeContentsKey {
  SerializeContentsKey() {} // default ctor private
  SerializeContentsKey(const SerializeContentsKey&) {} // copy ctor private

  // grant access to few methods
  template<class Model, class ModelParams>  
  template<class Archive>
  friend void DetQMCPT<Model,ModelParams>::saveContents(Archive& ar);
  template<class Model, class ModelParams>  
  template<class Archive>
  friend void DetQMCPT<Model,ModelParams>::loadContents(Archive& ar);
  template<class Model, class ModelParams>  
  template<class Archive>
  friend void DetQMCPT<Model,ModelParams>::serializeContentsCommon(Archive& ar);
};





template<class Model, class ModelParams>
void DetQMCPT<Model,ModelParams>::initFromParameters(const ModelParams& parsmodel_, const DetQMCParams& parsmc_,
                                                     const DetQMCPTParams& parspt_) {
    parsmodel = parsmodel_;
    parsmc = parsmc_;
    parspt = parspt_;

    parsmc.check();
    parspt.check();

    // Set up MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &processIndex);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    if (numProcesses != parspt.controlParameterValues.size()) {
        throw ConfigurationError("Number of processes " + numToString(numProcesses) +
                                 " does not match number of control parameter values " +
                                 numToString(parspt.controlParameterValues.size()));
    }

    // set up RNG
    if (parsmc.specified.count("rngSeed") == 0) {
        if (processIndex == 0) {
            std::cout << "No rng seed specified, will use std::time(0) determined at root process" << std::endl;
            parsmc.rngSeed = (uint32_t) std::time(0);
        }
        MPI_Bcast(&rngSeed, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
    }
    rng = RngWrapper(parsmc.rngSeed, processIndex);

    // set up control parameters for current process replica parameters
    parsmodel.set_exchange_parameter_value(parspt.controlParameterValues[processIndex]);

    // at rank 0 keep track of which process has which control parameter currently
    // and track exchange action contributions
    if (processIndex == 0) {
        // fill with 0, 1, 2, ... numProcesses-1
        current_process_par.resize(numProcesses);
        std::iota(current_process_par.begin(), current_process_par.end(), 0);
        current_par_process.resize(numProcesses);
        std::iota(current_par_process.begin(), current_par_process.end(), 0);
        exchange_action.resize(numProcesses, 0);
    }
    
    replica = createReplica<Model>(rng, parsmodel);

    //prepare metadata
    modelMeta = replica->prepareModelMetadataMap();
    mcMeta = parsmc.prepareMetadataMap();
    ptMeta = parspt.prepareMetadataMap();

    //prepare observable handlers
    auto scalarObs = replica->getScalarObservables();
    for (auto obsP = scalarObs.cbegin(); obsP != scalarObs.cend(); ++obsP) {
        obsHandlers.push_back(
            ObsPtr(new ScalarObservableHandler(*obsP, current_process_par, parsmc, parspt, modelMeta, mcMeta))
        );
    }

    auto vectorObs = replica->getVectorObservables();
    for (auto obsP = vectorObs.cbegin(); obsP != vectorObs.cend(); ++obsP) {
        vecObsHandlers.push_back(
            VecObsPtr(new VectorObservableHandler(*obsP, current_process_par, parsmc, parspt, modelMeta, mcMeta))
        );
    }
    auto keyValueObs = replica->getKeyValueObservables();
    for (auto obsP = keyValueObs.cbegin(); obsP != keyValueObs.cend(); ++obsP) {
        vecObsHandlers.push_back(
            VecObsPtr(new KeyValueObservableHandler(*obsP, current_process_par, parsmc, parspt, modelMeta, mcMeta))
        );
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
DetQMCPT<Model, ModelParams>::DetQMCPT(const ModelParams& parsmodel_, const DetQMCParams& parsmc_,
                                       const DetQMCPTParams& parspt_) :
    parsmodel(), parsmc(), parspt(),
    //proper initialization of default initialized members done in initFromParameters
    modelMeta(), mcMeta(), ptMeta(), rng(), replica(),
    obsHandlers(), vecObsHandlers(),
    sweepsDone(0), sweepsDoneThermalization(),
    swCounter(0),
    elapsedTimer(),     // start timing
    totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
    grantedWalltimeSecs(0), jobid(""),
    numProcesses(1), processIndex(0),
    current_process_par(1, 0),
    current_par_process(1, 0),
    exchange_action(1, 0)
{
    initFromParameters(parsmodel_, parsmc_, parspt_);
}

template<class Model, class ModelParams>
DetQMCPT<Model, ModelParams>::DetQMCPT(const std::string& stateFileName, const DetQMCParams& newParsmc,
                                       const DetQMCPTParams& newParspt) :
    parsmodel(), parsmc(), parspt(),
    //proper initialization of default initialized members done by loading from archive
    modelMeta(), mcMeta(), rng(), replica(),
    obsHandlers(), vecObsHandlers(),
    sweepsDone(), sweepsDoneThermalization(),
    swCounter(0),
    elapsedTimer(),     // start timing
    totalWalltimeSecs(0), walltimeSecsLastSaveResults(0),
    grantedWalltimeSecs(0), jobid(""),
    numProcesses(1), processIndex(0),
    current_process_par(1, 0),
    current_par_process(1, 0),
    exchange_action(1, 0)
{
    std::ifstream ifs;
    ifs.exceptions(std::ifstream::badbit | std::ifstream::failbit);
    ifs.open(stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ModelParams parsmodel_;
    DetQMCParams parsmc_;
    DetQMCPTParams parspt_;
    ia >> parsmodel_ >> parsmc_ >> parspt_;

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
    
    initFromParameters(parsmodel_, parsmc_, parspt_);
    loadContents(ia);

    std::cout << "\n"
              << "State of previous simulation has been loaded.\n"
              << "  sweepsDoneThermalization: " << sweepsDoneThermalization << "\n"
              << "  sweepsDone: " << sweepsDone << std::endl;
}

template<class Model, class ModelParams>
void DetQMCPT<Model, ModelParams>::saveState() {
    timing.start("saveState");

    //serialize state to file
    std::ofstream ofs;
    ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    ofs.open(parsmc.stateFileName.c_str(), std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << parsmodel << parsmc << parspt;
    saveContents(oa);

    //write out info about state of simulation to "info.dat"
    std::string commonInfoFilename = "info.dat";
    writeOnlyMetaData(commonInfoFilename, collectVersionInfo(),
                      "Collected information about this determinantal quantum Monte Carlo simulation",
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

    std::cout << "State has been saved." << std::endl;

    timing.stop("saveState");
}

template<class Model, class ModelParams>
DetQMCPT<Model, ModelParams>::~DetQMCPT() {
}


template<class Model, class ModelParams>
void DetQMCPT<Model, ModelParams>::run() {
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

    const std::string abortFilename1 = "ABORT." + jobid;
    const std::string abortFilename2 = "../" + abortFilename1;

    while (stage != F) {                //big loop
        // do we need to quit?
    	if (swCounter % 2 == 0) {
            bool stop_now = false;
            if (curWalltimeSecs() > grantedWalltimeSecs - SavetyMinutes*60) {
                std::cout << "Granted walltime will be exceeded in less than " << SavetyMinutes << " minutes.\n";
                stop_now = true;
            } else if (boost::filesystem::exists(abortFilename1) or
                       boost::filesystem::exists(abortFilename2)) {
                std::cout << "Found file " << abortFilename1 << ".\n";
                stop_now = true;
            }
            if (stop_now) {
                //close to exceeded walltime or we find that a file has been placed,
                //which signals us to abort this run for some other reason.
                //but only save state and exit if we have done an even
                //number of sweeps for ("economic") serialization guarantee [else do one sweep more]
                std::cout << "Save state / results and exit gracefully." << std::endl;
                if (stage == Stage::M) {
                    saveResults();
                }
                saveState();
                break;  //while
            }
    	}

        //thermalization & measurement stages | main work
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
                std::cout << std::endl;
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
            }
            ++sweepsDone;
            if (swCounter == parsmc.saveInterval) {
                std::cout << "  " << sweepsDone << " ... saving results and state ...";
                swCounter = 0;
                saveResults();
                saveState();
                std::cout << std::endl;
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

        //replica exchange
        if (stage == T or stage == M) {
            if (swCounter % parspt.exchangeInterval == 0) {
                // Gather exchange action contribution from replicas
                double localAction = replica->get_exchange_action_contribution();
                MPI_Gather( &localAction,           // send buf
                            1,
                            MPI_DOUBLE,
                            exchange_action.data(), // recv buf
                            1,
                            MPI_DOUBLE,
                            0,                      // root process
                            MPI_COMM_WORLD
                    );
                // serially walk through control parameters and propose exchange
                if (processIndex == 0) {
                    for (uint32_t cpi1 = 0; cpi1 < numReplicas - 1; ++cpi1) {
                        uint32_t cpi2 = cpi1 + 1;
                        double par1 = parspt.controlParameterValues[cpi1];
                        double par2 = parspt.controlParameterValues[cpi2];
                        uint32_t indexProc1 = current_par_process[par1];
                        uint32_t indexProc2 = current_par_process[par2];
                        double action1 = exchange_action[indexProc1];                    
                        double action2 = exchange_action[indexProc2];

                        num exchange_prob = get_replica_exchange_probability<Model>(par1, action1,
                                                                                    par2, action2);
                        if (exchange_prob >= 1 or rng.rand01() <= exchange_prob) {
                            // swap control parameters
                            current_process_par[indexProc1] = cpi2;
                            current_process_par[indexProc2] = cpi1;
                            current_par_process[cpi1] = indexProc2;
                            current_par_process[cpi2] = indexProc1;
                        }                    
                    }
                } // if (processIndex == 0)
                // distribute and update control parameter
                uint32_t new_param_index = 0;
                MPI_Scatter( current_process_par.data(), // send buf
                             1,
                             MPI_UINT32_T,
                             &new_param_index,           // recv buf
                             1
                             MPI_UINT32_T,
                             0,                          // root process
                             MPI_COMM_WORLD
                    );
                replica->set_exchange_parameter_value(
                    parspt.controlParameterValues[new_param_index]
                    );
                
            }
        } //replica exchange
        
    } // while (stage != F)
}






template<class Model, class ModelParams>
void DetQMCPT<Model, ModelParams>::saveResults() {
    timing.start("saveResults");

    outputResults(obsHandlers);
    for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
        (*p)->outputTimeseries();
    }
    outputResults(vecObsHandlers);

    timing.stop("saveResults");
}




#endif /* DETQMC_H_ */
