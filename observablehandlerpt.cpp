#include "observablehandlerpt.h"

void outputResults(const std::vector<std::unique_ptr<ScalarObservableHandlerPT>>& obsHandlers) {
    MPI_Comm_rank(MPI_COMM_WORLD, &processIndex);
    if (processIndex == 0) {
        typedef std::map<std::string, num> StringNumMap;
        typedef std::shared_ptr<StringNumMap> StringNumMapPtr;

        for (uint32_t cpi = 0; cpi < numProcesses; ++cpi) {
            StringNumMapPtr values(new StringNumMap);
            StringNumMapPtr errors(new StringNumMap);

            std::string subdirectory = "p" + numToString(cpi) + "_" +
                (*obsHandlers.begin())->mcparams.controlParameterName +
                numToString((*obsHandlers.begin())->mcparams.controlParameterValues[cpi]);

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
            output.addMeta("key", "observable");
            output.addHeaderText("observable\t value \t error");
            output.writeToFile(subdirectory + "/results.values");
        }
    }
}

void outputResults(const std::vector<std::unique_ptr<VectorObservableHandlerPT>>& obsHandlers) {
    MPI_Comm_rank(MPI_COMM_WORLD, &processIndex);
    if (processIndex == 0) {    
        typedef std::map<num, num> NumMap;
        typedef std::shared_ptr<NumMap> NumMapPtr;
        typedef DataMapWriter<num,num> NumMapWriter;

        for (uint32_t cpi = 0; cpi < numProcesses; ++cpi) {
            std::string subdirectory = "p" + numToString(cpi) + "_" +
                (*obsHandlers.begin())->mcparams.controlParameterName +
                numToString((*obsHandlers.begin())->mcparams.controlParameterValues[cpi]);

            for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
                const std::unique_ptr<VectorObservableHandler>& obsptr = *p;
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
                output.addMeta("key", obsptr->indexName);
                output.addMeta("observable", obsptr->name);
                output.addHeaderText("key\t value \t error");
                output.writeToFile(subdirectory + "/results-" + obsptr->name + ".values");
            }
        }
    }
}
