#ifndef MPIOBSERVABLEHANDLERPT_H
#define MPIOBSERVABLEHANDLERPT_H

// manage measurements of an observable, gather measurement values
// from various replicas, calculate expectation values and jackknife
// error bars; optionally store time series

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <armadillo>
#include "tools.h"
#include "detqmcparams.h"
#include "detqmcptparams.h"
#include "observable.h"
#include "metadata.h"
#include "dataserieswritersucc.h"
#include "datamapwriter.h"
#include "statistics.h"

class SerializeContentsKey;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/vector.hpp"
#include "boost/serialization/export.hpp"
#include "boost_serialize_uniqueptr.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "boost/mpi.hpp"
#pragma GCC diagnostic pop


template <typename ObsType>
class ObservableHandlerPTCommon {
public:
    ObservableHandlerPTCommon(
        const Observable<ObsType>& localObservable,
        const std::vector<int>& current_process_par,
        const DetQMCParams& simulationParameters,
        const DetQMCPTParams& ptParams,
        const MetadataMap& metadataToStoreModel,
        const MetadataMap& metadataToStoreMC,
        const MetadataMap& metadataToStorePT,
        ObsType zeroValue = ObsType());
    virtual ~ObservableHandlerPTCommon() { }
    //To be called from rank 0:
    //return [mean value, error] at end of simulation
    //if jkBlockCount <= 1, only estimate an error (using variance()) if the whole
    //timeseries is in memory
    //
    //return [mean value, 0] if this is called earlier.
    //return [0, 0] if this is called from rank != 0.
    std::tuple<ObsType,ObsType> evaluateJackknife(int control_parameter_index) const;    
protected:
    // this is to be called at rank0 by method insertValue() of a
    // derived class
    void handleValues(uint32_t curSweep);

    Observable<ObsType> localObs; // Handle containing the last measured value for
                                  // the observable, at the local replica    
    
    const std::string& name;    //reference to name in obs
    ObsType zero;               //an instance of ObsType that works like the number zero
                                //for addition -- this is not totally trivial for vector
                                //valued observables

    DetQMCParams mcparams;
    DetQMCPTParams ptparams;
    MetadataMap metaModel;
    MetadataMap metaMC;
    MetadataMap metaPT;    
    uint32_t jkBlockCount;
    uint32_t jkBlockSizeSweeps;

    uint32_t lastSweepLogged;
    uint32_t countValues;

    //MPI specifics
    int numProcesses;      //total number of parallel processes
    int processIndex;      //MPI-rank of the current process
    
    // root process specifics:
    // -----------------------
    // This is indexed by process number, and contains the control
    // parameter value currently associated to each process. It is a reference
    // to the vector held and kept uptodate by DetQMCPT. It is only sensible
    // at the root process.
    const std::vector<int>& process_par;
    //This is the receive buffer of most recently measured observable
    //values for each replica, ordered by process.  The method
    //insertValue() of the derived classes below fill this, then
    //call handleValues() defined above in this base class.
    std::vector<ObsType> process_cur_value;
    // These vectors contain one entry per control parameter value,
    // they are indexed by the control parameter index.
    std::vector<std::vector<ObsType>> par_jkBlockValues; // running counts of jackknife block values
    std::vector<ObsType> par_total; // running accumulation regardless of jackknife block
    // for each control parameter value [sorted by index]: store a
    // separate MetadataMap with just that entry replaced
    std::vector<MetadataMap> par_metaModel;
public:
    // only functions that can pass the key to this function have access
    // -- in this way access is granted only to DetQMC::serializeContents
    template<class Archive>
    void serializeContents(SerializeContentsKey const &, Archive &ar) {
        ar & lastSweepLogged;
        ar & countValues;
        ar & par_jkBlockValues;
        ar & par_total;
        // par_metaModel does not need to be serialized -- is reset
        // upon initialization
    }
};

template <typename ObsType>
ObservableHandlerPTCommon<ObsType>::ObservableHandlerPTCommon(
    const Observable<ObsType>& localObservable,
    const std::vector<int>& current_process_par,
    const DetQMCParams& simulationParameters,
    const DetQMCPTParams& ptParams,
    const MetadataMap& metadataToStoreModel,
    const MetadataMap& metadataToStoreMC,
    const MetadataMap& metadataToStorePT,
    ObsType zeroValue)
    : localObs(localObservable), name(localObs.name),
      zero(zeroValue),          //ObsType() may not be a valid choice!
      mcparams(simulationParameters),
      ptparams(ptParams),
      metaModel(metadataToStoreModel),
      metaMC(metadataToStoreMC),
      metaPT(metadataToStorePT),
      jkBlockCount(mcparams.jkBlocks),
      jkBlockSizeSweeps(mcparams.sweeps / jkBlockCount),
      lastSweepLogged(0),
      countValues(0),
      numProcesses(1),
      processIndex(0),
      process_par(current_process_par),
      process_cur_value(),
      par_jkBlockValues(),
      par_total(),
      par_metaModel()
{
    boost::mpi::communicator world;
    processIndex = world.rank();
    numProcesses = world.size();
    assert(int(ptparams.controlParameterValues.size()) == numProcesses);
    if (processIndex == 0) {
        process_cur_value.resize(numProcesses, zero);
        par_jkBlockValues.resize(numProcesses, std::vector<ObsType>(jkBlockCount, zero));
        par_total.resize(numProcesses, zero);
        par_metaModel.resize(numProcesses, metaModel);
        for (int cpi = 0; cpi < numProcesses; ++cpi) {
            par_metaModel[cpi][ptparams.controlParameterName] =
                numToString(ptparams.controlParameterValues[cpi]);
        }
    }
}


template <typename ObsType>
void ObservableHandlerPTCommon<ObsType>::handleValues(uint32_t curSweep) {
    if (processIndex == 0) {
        uint32_t curJkBlock = curSweep / jkBlockSizeSweeps;
        for (int p_i = 0; p_i < numProcesses; ++p_i) {
            int controlParameterIndex = process_par[p_i];
            for (uint32_t jb = 0; jb < jkBlockCount; ++jb) {
                if (jb != curJkBlock) {
                    par_jkBlockValues[controlParameterIndex][jb] += process_cur_value[p_i];
                }
            }
            par_total[controlParameterIndex] += process_cur_value[p_i];
        }
    }
    ++countValues;
    lastSweepLogged = curSweep;
}

template <typename ObsType>
std::tuple<ObsType,ObsType> ObservableHandlerPTCommon<ObsType>::evaluateJackknife(
    int cpi) const {
    if (processIndex != 0) {
        return std::make_tuple(zero, zero);
    } else {
        ObsType mean = par_total[cpi] / countValues;
        ObsType error = zero;
        if (mcparams.sweeps - lastSweepLogged <= mcparams.measureInterval) {
            //after the first sweep lastSweepLogged==1 and so on --> here the simulation is finished.
            //we can only calculate an error estimate if we have multiple jackknife blocks
            if (jkBlockCount > 1 and not mcparams.sweepsHasChanged) {
                uint32_t jkBlockSizeSamples = countValues / jkBlockCount;
                uint32_t jkTotalSamples = countValues - jkBlockSizeSamples;
//              std::cout << jkTotalSamples << std::endl;
                std::vector<ObsType> jkBlockAverages = par_jkBlockValues[cpi];   //copy
                for (uint32_t jb = 0; jb < jkBlockCount; ++jb) {
                    jkBlockAverages[jb] /= jkTotalSamples;
                }
                error = jackknife(jkBlockAverages, mean, zero);
            }
        }
        return std::make_tuple(mean, error);
    }
}


// Below here we explicitly use `double` instead of `num` because the
// MPI calls explicitly use MPI_DOUBLE.  This can easily be extended
// if some other floating precision type is ever to be used for num.



//specialized ObservableHandlerPT that uses num as a value type
// -- can store time series, can be output into common files "results*.values"
//    for all scalar observables [in subdirectories]
class ScalarObservableHandlerPT : public ObservableHandlerPTCommon<double> {
public:
    ScalarObservableHandlerPT(const ScalarObservable& localObservable,
                              const std::vector<int>& current_process_par,
                              const DetQMCParams& simulationParameters,
                              const DetQMCPTParams& ptParams,
                              const MetadataMap& metadataToStoreModel,
                              const MetadataMap& metadataToStoreMC,
                              const MetadataMap& metadataToStorePT
        );

    // Log newly measured observable value at each replica via the the
    // reference contained in this->obs, pass the number of the
    // current sweep.  Measurements from all replicas are gathered at
    // the root process.  Measurements do not need to be stored at
    // every sweep, but the number of skipped sweeps must be constant.
    void insertValue(uint32_t curSweep);    

    //If we don't have multiple jackknife blocks and the whole timeseries is stored
    //in memory, this can also give a naive variance estimate for the error
    std::tuple<double, double> evaluateJackknife(int control_parameter_index) const;

    //update timeseries files, discard batch of
    //data written to files from memory
    void outputTimeseries(); 

    friend void outputResults(
        const std::vector<std::unique_ptr<ScalarObservableHandlerPT>>& obsHandlers);
protected:
    //in addition to base class functionality supports adding to the timeseries buffers
    void handleValues(uint32_t curSweep);
    
    // time series entries added since last call to writeData(),
    // for each control parameter value
    std::vector<std::vector<double>> par_timeseriesBuffer;

    std::vector<std::unique_ptr<DoubleVectorWriterSuccessive>> par_storage;
    std::vector<char> par_storageFileStarted; // avoiding vector<bool>, but using it equivalently [true/false]
public:
    // only functions that can pass the key to this function have access
    // -- in this way access is granted only to DetQMC::serializeContents
    template<class Archive>
    void serializeContents(SerializeContentsKey const &sck, Archive &ar) {
        ObservableHandlerPTCommon<double>::serializeContents(sck, ar);
        ar & par_timeseriesBuffer;
        ar & par_storageFileStarted;
        //*storage should not need to be serialized.  It will always write to the end
        //of the timeseries file it finds at construction.
    }
};




//Vector valued observables.  We use Armadillo vectors as they support arithmetics.
//A fixed vector size must be specified at initialization. This indexes the vector from 0 to
//the vector size.

class VectorObservableHandlerPT : public ObservableHandlerPTCommon<arma::Col<double>> {
public:
    VectorObservableHandlerPT(const VectorObservable& localObservable,
                              const std::vector<int>& current_process_par,
                              const DetQMCParams& simulationParameters,
                              const DetQMCPTParams& ptParams,
                              const MetadataMap& metadataToStoreModel,
                              const MetadataMap& metadataToStoreMC,
                              const MetadataMap& metadataToStorePT);
    void insertValue(uint32_t curSweep);    //compare to ScalarObservableHandlerPT function    
    uint32_t getVectorSize() {
        return vsize;
    }
    friend void outputResults(
        const std::vector<std::unique_ptr<VectorObservableHandlerPT>>& obsHandlers);
protected:
    uint32_t vsize;
    arma::Col<double> indexes;
    std::string indexName;
    // at rank 0 this holds contiguous memory where the vector data
    // gathered from all replicas is stored
    std::vector<double> mpi_gather_buffer; 
};





//Vector indexed by arbitrary key
class KeyValueObservableHandlerPT : public VectorObservableHandlerPT {
public:
    KeyValueObservableHandlerPT(const KeyValueObservable& observable,
                                const std::vector<int>& current_process_par,                              
                                const DetQMCParams& simulationParameters,
                                const DetQMCPTParams& ptParams,
                                const MetadataMap& metadataToStoreModel,
                                const MetadataMap& metadataToStoreMC,
                                const MetadataMap& metadataToStorePT) :
        VectorObservableHandlerPT(observable, current_process_par,
                                  simulationParameters, ptParams,
                                  metadataToStoreModel, metadataToStoreMC,
                                  metadataToStorePT) {
        //this code is convenient, but sets the vector indexes twice upon construction
        indexes = observable.keys;
        indexName = observable.keyName;
    }
};


//Write expectation values and error bars for all observables to a file
//take metadata to store from the first entry in obsHandlers
//This is to be called by rank 0
void outputResults(const std::vector<std::unique_ptr<ScalarObservableHandlerPT>>& obsHandlers);

//write the results for each vector observable into a seperate file
void outputResults(const std::vector<std::unique_ptr<VectorObservableHandlerPT>>& obsHandlers);



    
#endif /* MPIOBSERVABLEHANDLERPT_H */
