#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/filesystem.hpp"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "boost/mpi.hpp"
#pragma GCC diagnostic pop

#include "mpiobservablehandlerpt.h"

namespace fs = boost::filesystem;
namespace mpi = boost::mpi;


ScalarObservableHandlerPT::ScalarObservableHandlerPT(
        const ScalarObservable& localObservable,
        const std::vector<int>& current_process_par,
        const DetQMCParams& simulationParameters,
        const DetQMCPTParams& ptParams,
        const MetadataMap& metadataToStoreModel,
        const MetadataMap& metadataToStoreMC,
        const MetadataMap& metadataToStorePT)
    : ObservableHandlerPTCommon<double>(localObservable, current_process_par,
                                        simulationParameters, ptParams,
                                        metadataToStoreModel, metadataToStoreMC,
                                        metadataToStorePT),
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
            for (int p_i = 0; p_i < numProcesses; ++p_i) {
                int controlParameterIndex = process_par[p_i];
                par_timeseriesBuffer[controlParameterIndex].push_back(process_cur_value[p_i]);
            }

        }
    }
    ObservableHandlerPTCommon<double>::handleValues(curSweep);
}

void ScalarObservableHandlerPT::insertValue(uint32_t curSweep) {
    //MPI: gather the value of localObs from each replica in the
    //buffer at the root process: process_cur_value
    
    // MPI_Gather( const_cast<double*>(&(localObs.valRef.get())), // sendbuf :
    //             // pass memory address of what localObs references | need to cast away const for mpi < 3.0
    //             1,                        // sendcount
    //             MPI_DOUBLE,               // sendtype
    //             process_cur_value.data(), // recvbuf: pointer to internal memory of vector [at root process]
    //             1,                        // recvcount
    //             MPI_DOUBLE,               // recvtype
    //             0,                        // root
    //             MPI_COMM_WORLD            // comm
    //     );
    mpi::communicator world;
    mpi::gather(world,
                localObs.valRef.get(), // send: what localObs references
                process_cur_value,     // recv
                0);
                
    this->handleValues(curSweep);
}



std::tuple<double, double> ScalarObservableHandlerPT::evaluateJackknife(
    int control_parameter_index) const {
    if (processIndex != 0) {
        return std::make_tuple(0.0, 0.0);
    } else {
        double mean;
        double error;

        std::tie(mean, error) = ObservableHandlerPTCommon<double>::evaluateJackknife(control_parameter_index);

        if (jkBlockCount <= 1 and par_timeseriesBuffer[control_parameter_index].size() == countValues) {
            error = std::sqrt(variance(par_timeseriesBuffer[control_parameter_index], mean));
        }
        return std::make_tuple(mean, error);
    }
}


void ScalarObservableHandlerPT::outputTimeseries() {
    //TODO: float precision
    if (processIndex == 0 and mcparams.timeseries) {
        for (int p_i = 0; p_i < numProcesses; ++p_i) {
            int cpi = process_par[p_i];

            std::string subdirectory = "p" + numToString(cpi) + "_" +
                ptparams.controlParameterName +
                numToString(ptparams.controlParameterValues[cpi]);
            fs::create_directories(fs::path(subdirectory));            
        
            if (not par_storage[cpi]) {
                std::string filename = (fs::path(subdirectory) / fs::path(name + ".series")).string();
                if (not par_storageFileStarted[cpi]) {
                    par_storage[cpi] =
                        std::unique_ptr<DoubleVectorWriterSuccessive>(
                            new DoubleVectorWriterSuccessive(filename,
                                                             false // create a new file
                                ));
                    par_storage[cpi]->addHeaderText("Timeseries for observable " + name);
                    par_storage[cpi]->addMetadataMap(par_metaModel[cpi]);
                    par_storage[cpi]->addMetadataMap(metaMC);
                    par_storage[cpi]->addMetadataMap(metaPT);
                    par_storage[cpi]->addMeta("observable", name);
                    par_storage[cpi]->writeHeader();
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


VectorObservableHandlerPT::VectorObservableHandlerPT(const VectorObservable& localObservable,
                                                     const std::vector<int>& current_process_par,
                                                     const DetQMCParams& simulationParameters,
                                                     const DetQMCPTParams& ptParams,
                                                     const MetadataMap& metadataToStoreModel,
                                                     const MetadataMap& metadataToStoreMC,
                                                     const MetadataMap& metadataToStorePT)
: ObservableHandlerPTCommon<arma::Col<double>>(
    localObservable, current_process_par,
    simulationParameters, ptParams,
    metadataToStoreModel, metadataToStoreMC,
    metadataToStorePT,
    arma::zeros<arma::Col<double>>(localObservable.vectorSize)),
    vsize(localObservable.vectorSize), indexes(vsize), indexName("site"),
    mpi_gather_buffer()
{
    for (uint32_t counter = 0; counter < vsize; ++counter) {
        indexes[counter] = counter;
    }
    if (processIndex == 0) {
        mpi_gather_buffer.resize(numProcesses * vsize, 0.0);
    } else {
        mpi_gather_buffer.resize(1, 0.0);
    }
}

void VectorObservableHandlerPT::insertValue(uint32_t curSweep) {
    mpi::communicator world;
    //MPI: gather the value of localObs from each replica in the
    //buffer at the root process: process_cur_value
    assert(vsize == localObs.valRef.get().n_elem);
    // MPI_Gather( const_cast<double*>(localObs.valRef.get().memptr()),  // sendbuf
    //             // pass arma vector data behind localObs reference,
    //             // need to cast away const for MPI < 3.0
    //             vsize,                           // sendcount
    //             MPI_DOUBLE,                      // sendtype
    //             mpi_gather_buffer.data(),        // recvbuf: pointer to internal memory of vector [at root process]
    //             vsize,                           // recvcount
    //             MPI_DOUBLE,                      // recvtype
    //             0,                               // root
    //             MPI_COMM_WORLD                   // comm
    //     );
    mpi::gather(world,
                localObs.valRef.get().memptr(), // send: pass arma vector data behind localObs reference
                vsize,                          // sendcount
                mpi_gather_buffer,              // recv
                0                               // root
        );

    // use the contiguous memory of mpi_gather_buffer to hold the data
    // for the Armadillo vectors used for the individual replica measurements
    // at the root process
    if (processIndex == 0) {
        assert(mpi_gather_buffer.size() == numProcesses * vsize);
        for (int p_i = 0; p_i < numProcesses; ++p_i) {
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




void outputResults(const std::vector<std::unique_ptr<ScalarObservableHandlerPT>>& obsHandlers) {
    boost::mpi::communicator world;
    int processIndex = world.rank();
    int numProcesses = world.size();
    if (processIndex == 0 and obsHandlers.size() > 0) {
        typedef std::map<std::string, num> StringNumMap;
        typedef std::shared_ptr<StringNumMap> StringNumMapPtr;

        for (int cpi = 0; cpi < numProcesses; ++cpi) {
            StringNumMapPtr values(new StringNumMap);
            StringNumMapPtr errors(new StringNumMap);

            std::string subdirectory = "p" + numToString(cpi) + "_" +
                (*obsHandlers.begin())->ptparams.controlParameterName +
                numToString((*obsHandlers.begin())->ptparams.controlParameterValues[cpi]);
            fs::create_directories(fs::path(subdirectory));            

            for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
                num val, err;
                std::tie(val, err) = (*p)->evaluateJackknife(cpi);
                std::string obsname = (*p)->name;
                (*values)[obsname] = val;
                (*errors)[obsname] = err;
            }
            DataMapWriter<std::string, num> output;
            output.setData(values);
            output.setErrors(errors);
            output.addHeaderText("Monte Carlo results for observable expectation values");
            output.addMetadataMap((*obsHandlers.begin())->par_metaModel[cpi]);
            output.addMetadataMap((*obsHandlers.begin())->metaMC);
            output.addMetadataMap((*obsHandlers.begin())->metaPT);            
            output.addMeta("key", "observable");
            output.addHeaderText("observable\t value \t error");
            output.writeToFile((fs::path(subdirectory) / fs::path("results.values")).string());
        }
    }
}

void outputResults(const std::vector<std::unique_ptr<VectorObservableHandlerPT>>& obsHandlers) {
    boost::mpi::communicator world;
    int processIndex = world.rank();
    int numProcesses = world.size();
    if (processIndex == 0 and obsHandlers.size() > 0) {    
        typedef std::map<num, num> NumMap;
        typedef std::shared_ptr<NumMap> NumMapPtr;
        typedef DataMapWriter<num,num> NumMapWriter;
        
        for (int cpi = 0; cpi < numProcesses; ++cpi) {
            std::string subdirectory = "p" + numToString(cpi) + "_" +
                (*obsHandlers.begin())->ptparams.controlParameterName +
                numToString((*obsHandlers.begin())->ptparams.controlParameterValues[cpi]);

            for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
                const std::unique_ptr<VectorObservableHandlerPT>& obsptr = *p;
                uint32_t numberIndexes = obsptr->getVectorSize();
                arma::Col<num> values, errors;
                std::tie(values, errors) = obsptr->evaluateJackknife(cpi);
                NumMapPtr valmap(new NumMap);
                NumMapPtr errmap(new NumMap);
                for (uint32_t counter = 0; counter < numberIndexes; ++counter) {
                    num index = obsptr->indexes[counter];
                    valmap->insert(std::make_pair(index, values[counter]));
                    errmap->insert(std::make_pair(index, errors[counter]));
                }
                NumMapWriter output;
                output.setData(valmap);
                output.setErrors(errmap);
                output.addHeaderText("Monte Carlo results for vector observable " + obsptr->name +
                                     " expectation values");
                output.addMetadataMap(obsptr->par_metaModel[cpi]);
                output.addMetadataMap(obsptr->metaMC);
                output.addMetadataMap(obsptr->metaPT);                
                output.addMeta("key", obsptr->indexName);
                output.addMeta("observable", obsptr->name);
                output.addHeaderText("key\t value \t error");
                output.writeToFile(subdirectory + "/results-" + obsptr->name + ".values");
            }
        }
    }
}
