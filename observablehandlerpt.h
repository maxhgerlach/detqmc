#ifndef OBSERVABLEHANDLERPT_H
#define OBSERVABLEHANDLERPT_H

// manage measurements of an observable, gather measurement values
// from various replicas, calculate expectation values and jackknife
// error bars; optionally store time series

#include <memory>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <armadillo>
#include "detqmcparams.h"
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
#pragma GCC diagnostic pop


template <typename ObsType>
class ObservableHandlerPTCommon {
public:
    ObservableHandlerPTCommon(
        const Observable<ObsType>& localObservable,
        const std::vector<uint32_t>& current_process_par;
        const DetQMCParams& simulationParameters,
        const DetQMCPTParams& ptParams,
        const MetadataMap& metadataToStoreModel,
        const MetadataMap& metadataToStoreMC,
        ObsType zeroValue = ObsType());
    virtual ~ObservableHandlerCommon() { }
    //To be called from rank 0:
    //return [mean value, error] at end of simulation
    //if jkBlockCount <= 1, only estimate an error (using variance()) if the whole
    //timeseries is in memory
    //
    //return [mean value, 0] if this is called earlier.
    //return [0, 0] if this is called from rank != 0.
    std::tuple<ObsType,ObsType> evaluateJackknife(uint32_t control_parameter_index) const;    
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
    uint32_t jkBlockCount;
    uint32_t jkBlockSizeSweeps;

    uint32_t lastSweepLogged;
    uint32_t countValues;

    //MPI specifics
    uint32_t numProcesses;      //total number of parallel processes
    uint32_t processIndex;      //MPI-rank of the current process
    
    // root process specifics:
    // -----------------------
    // This is indexed by process number, and contains the control
    // parameter value currently associated to each process. It is a reference
    // to the vector held and kept uptodate by DetQMCPT. It is only sensible
    // at the root process.
    const std::vector<uint32_t>& process_par;
    //This is the receive buffer of most recently measured observable
    //values for each replica, ordered by process.  The method
    //insertValue() of the derived classes below fill this, than
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
    }
};

template <typename ObsType>
ObservableHandlerPTCommon<ObsType>::ObservableHandlerPTCommon(
    const Observable<ObsType>& localObservable,
    const std::vector<uint32_t>& current_process_par;
    const DetQMCParams& simulationParameters,
    const DetQMCPTParams& ptParams,
    const MetadataMap& metadataToStoreModel,
    const MetadataMap& metadataToStoreMC,
    ObsType zeroValue = ObsType())
    : localObs(localObservable), name(localObs.name),
      zero(zeroValue),          //ObsType() may not be a valid choice!
      mcparams(simulationParameters),
      ptparams(ptParams),
      metaModel(metadataToStoreModel),
      metaMC(metadataToStoreMC),
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
      par_metaModel(),
{
    MPI_Comm_rank(MPI_COMM_WORLD, &processIndex);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    assert(ptparams.controlParameterValues.size() == numProcesses);
    if (processIndex == 0) {
        process_cur_value.resize(numProcesses, zero);
        par_jkBlockValues.resize(numProcesses, std::vector<ObsType>(jkBlockCount, zero));
        par_total.resize(numProcesses, zero);
        par_metaModel.resize(numProcesses, metaModel);
        for (uint32_t cpi = 0; cpi < numProcesses; ++cpi) {
            par_metaModel[cpi][ptparams.controlParameterName] =
                numToString(ptparams.controlParameterValues[cpi]);
        }
    }
}


template <typename ObsType>
void ObservableHandlerPTCommon<ObsType>::handleValues(uint32_t curSweep) {
    if (processIndex == 0) {
        uint32_t curJkBlock = curSweep / jkBlockSizeSweeps;
        for (uint32_t p_i = 0; p_i < numProcesses; ++p_i) {
            uint32_t controlParameterIndex = process_par[p_i];
            for (uint32_t jb = 0; jb < jkBlockCount; ++jb) {
                if (jb != curJkBlock) {
                    par_jkBlockValues[controlParameterIndex][jb] += process_cur_value[p_i];
                }
            }
            par_total[controlParameterIndex] += value;
        }
    }
    ++countValues;
    lastSweepLogged = curSweep;
}

template <typename ObsType>
std::tuple<ObsType,ObsType> ObservableHandlerPTCommon<ObsType>::evaluateJackknife(
    uint32_t cpi
    ) const {
    if (processIndex != 0) {
        return std::make_tuple(zero, zero);
    } else {
        ObsType mean = total[cpi] / countValues;
        ObsType error = zero;
        if (mcparams.sweeps - lastSweepLogged <= mcparams.measureInterval) {
            //after the first sweep lastSweepLogged==1 and so on --> here the simulation is finished.
            //we can only calculate an error estimate if we have multiple jackknife blocks
            if (jkBlockCount > 1 and not mcparams.sweepsHasChanged) {
                uint32_t jkBlockSizeSamples = countValues / jkBlockCount;
                uint32_t jkTotalSamples = countValues - jkBlockSizeSamples;
//              std::cout << jkTotalSamples << std::endl;
                std::vector<ObsType> jkBlockAverages = jkBlockValues[cpi];   //copy
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
                              const std::vector<uint32_t>& current_process_par;
                              const DetQMCParams& simulationParameters,
                              const DetQMCPTParams& ptParams,
                              const MetadataMap& metadataToStoreModel,
                              const MetadataMap& metadataToStoreMC);

    // Log newly measured observable value at each replica via the the
    // reference contained in this->obs, pass the number of the
    // current sweep.  Measurements from all replicas are gathered at
    // the root process.  Measurements do not need to be stored at
    // every sweep, but the number of skipped sweeps must be constant.
    void insertValue(uint32_t curSweep);    

    //If we don't have multiple jackknife blocks and the whole timeseries is stored
    //in memory, this can also give a naive variance estimate for the error
    std::tuple<double, double> evaluateJackknife(uint32_t control_parameter_index) const;

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
        ObservableHandlerCommon<double>::serializeContents(sck, ar);
        ar & par_timeseriesBuffer;
        ar & par_storageFileStarted;
        //*storage should not need to be serialized.  It will always write to the end
        //of the timeseries file it finds at construction.
    }
};


ScalarObservableHandlerPT::ScalarObservableHandlerPT(
    const ScalarObservable& localObservable,
    const std::vector<uint32_t>& current_process_par;
    const DetQMCParams& simulationParameters,
    const DetQMCPTParams& ptParams,
    const MetadataMap& metadataToStoreModel,
    const MetadataMap& metadataToStoreMC)
    : ObservableHandlerCommon<double>(localObservable, current_process_par,
                                      simulationParameters, ptParams,
                                      metadataToStoreModel, metadataToStoreMC),
      par_timeseriesBuffer(),       
      par_storage(),    
      par_storageFileStarted()
{
    if (processIndex == 0) {
        //by default one empty vector for each control parameter value
        //(as many as there are processes)
        par_timeseriesBuffer.resize(numProcesses, std::vector<double>());
        //initialize to a vector of something like a nullptr
        par_storage.resize(numProcesses);
        //boolean false stored as char
        par_storageFileStarted.resize(numProcesses, false);
    }
}


//to be called by insertvalue
void ScalarObservableHandlerPT::handleValues(uint32_t curSweep) {
    if (processIndex == 0) {
        if (mcparams.timeseries) {
            for (uint32_t p_i = 0; p_i < numProcesses; ++p_i) {
                uint32_t controlParameterIndex = process_par[p_i];
                timeseriesBuffer[controlParameterIndex].push_back(process_cur_value[p_i]);
            }

        }
    }
    ObservableHandlerCommon<double>::handleValues(curSweep);
}

void ScalarObservableHandlerPT::insertValue(uint32_t curSweep) {
    //MPI: gather the value of localObs from each replica in the
    //buffer at the root process: process_cur_value
    MPI_Gather( &(localObs.valRef.get()), // sendbuf: pass memory address of what localObs references
                1,                        // sendcount
                MPI_DOUBLE,               // sendtype
                process_cur_value.data(), // recvbuf: pointer to internal memory of vector [at root process]
                1,                        // recvcount
                MPI_DOUBLE,               // recvtype
                MPI_COMM_WORLD            // comm
        );
                
    this->handleValues(curSweep);
}



std::tuple<double, double> ScalarObservableHandlerPT::evaluateJackknife(
    uint32_t control_parameter_index) const {
    if (processIndex != 0) {
        return std::make_tuple(0.0, 0.0);
    } else {
        double mean;
        double error;

        std::tie(mean, error) = ObservableHandlerCommon<double>::evaluateJackknife(control_parameter_index);

        if (jkBlockCount <= 1 and timeseriesBuffer[control_parameter_index].size() == countValues) {
            error = std::sqrt(variance(timeseriesBuffer[control_parameter_index], mean));
        }
        return std::make_tuple(mean, error);
    }
}


void ScalarObservableHandlerPT::outputTimeseries() {
    //TODO: float precision
    if (processIndex == 0 and mcparams.timeseries) {
        for (uint32_t p_i = 0; p_i < numProcesses; ++p_i) {
            uint32_t cpi = process_par[p_i];

            std::string subdirectory = "p" + numToString(cpi) + "_" +
                mcparams.controlParameterName +
                numToString(mcparams.controlParameterValues[cpi]);
        
            if (not par_storage[cpi]) {
                std::string filename = subdirectory + "/" + name + ".series";
                if (par_storageFileStarted[cpi]) {
                    par_storage[cpi] =
                        std::unique_ptr<DoubleVectorWriterSuccessive>(
                            new DoubleVectorWriterSuccessive(filename,
                                                             false // create a new file
                                ));
                    par_storage->addHeaderText("Timeseries for observable " + name);
                    par_storage->addMetadataMap(par_metaModel[cpi]);
                    par_storage->addMetadataMap(metaMC);
                    par_storage->addMeta("observable", name);
                    par_storage->writeHeader();
                    par_storageFileStarted[cpi] = true;
                } else {
                    par_storage[cpi] = std::unique_ptr<DoubleVectorWriterSuccessive>(
                        new DoubleVectorWriterSuccessive(filename,
                                                         true// append to file
                            ));
                }
            }
            //append last batch of measurements
            par_storage[cpi]->writeData(par_timeseriesBuffer[cpi]);
            //no need to keep it in memory anymore
            par_timeseriesBuffer[cpi].resize(0);
        }
    }
}


//Vector valued observables.  We use Armadillo vectors as they support arithmetics.
//A fixed vector size must be specified at initialization. This indexes the vector from 0 to
//the vector size.

class VectorObservableHandlerPT : public ObservableHandlerPTCommon<arma::Col<double>> {
public:
    VectorObservableHandlerPT(const VectorObservable& localObservable,
                              const std::vector<uint32_t>& current_process_par;
                              const DetQMCParams& simulationParameters,
                              const DetQMCPTParams& ptParams,
                              const MetadataMap& metadataToStoreModel,
                              const MetadataMap& metadataToStoreMC);
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

VectorObservableHandlerPT::VectorObservableHandlerPT(const VectorObservable& localObservable,
                                                     const std::vector<uint32_t>& current_process_par;
                                                     const DetQMCParams& simulationParameters,
                                                     const DetQMCPTParams& ptParams,
                                                     const MetadataMap& metadataToStoreModel,
                                                     const MetadataMap& metadataToStoreMC)
    : ObservableHandlerPTCommon<arma::Col<double>>(
        localObservable, current_process_par,
        simulationParameters, ptParams,
        metadataToStoreModel, metadataToStoreMC,
        arma::zeros<arma::Col<double>>(observable.vectorSize)),
    vsize(observable.vectorSize), indexes(vsize), indexName("site"),
    mpi_gather_buffer()
{
    for (uint32_t counter = 0; counter < vsize; ++counter) {
        indexes[counter] = counter;
    }
    if (processIndex == 0) {
        mpi_gather_buffer.resize(numProcesses * vsize, 0.0);
    }
}

void VectorObservableHandlerPT::insertValue(uint32_t curSweep) {
    //MPI: gather the value of localObs from each replica in the
    //buffer at the root process: process_cur_value
    assert(vsize == localObs.valRef.n_elem);
    MPI_Gather( localObs.valRef.memptr(),  // sendbuf
                                           // pass arma vector data behind localObs reference
                vsize,                     // sendcount
                MPI_DOUBLE,                // sendtype
                mpi_gather_buffer.data(),  // recvbuf: pointer to internal memory of vector [at root process]
                vsize,                     // recvcount
                MPI_DOUBLE,                // recvtype
                MPI_COMM_WORLD             // comm
        );

    // use the contiguous memory of mpi_gather_buffer to hold the data
    // for the Armadillo vectors used for the individual replica measurements
    // at the root process
    if (processIndex == 0) {
        assert(mpi_gather_buffer.size() == numProcesses * vsize);
        for (uint32_t p_i = 0; p_i < numProcesses; ++p_i) {
            assert(process_cur_value[p_i].n_elem == vsize);
            process_cur_value[p_i] = arma::Col<double>(
                    mpi_gather_buffer.data() + p_i * vsize, // aux_mem*  [typed pointer: we do not need a sizeof(double) factor]
                    vsize,          // number_of_elements
                    false,          // copy_aux_mem [will continue to use the mpi_gather_buffer memory]
                    true            // strict [vector will remain bound to this memory for its lifetime]
                );
        }
    }
                
    this->handleValues(curSweep);
}




//Vector indexed by arbitrary key
class KeyValueObservableHandlerPT : public VectorObservableHandlerPT {
public:
    KeyValueObservableHandlerPT(const KeyValueObservable& observable,
                                const std::vector<uint32_t>& current_process_par,                              
                                const DetQMCParams& simulationParameters,
                                const DetQMCPTParams& ptParams,
                                const MetadataMap& metadataToStoreModel,
                                const MetadataMap& metadataToStoreMC) :
        VectorObservableHandlerPT(observable, current_process_par,
                                  simulationParameters, ptParams,
                                  metadataToStoreModel, metadataToStoreMC) {
        //this code is convenient but sets the vector indexes twice upon construction
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



    
#endif /* OBSERVABLEHANDLERPT_H */
