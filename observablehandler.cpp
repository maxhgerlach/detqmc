/*
 * observablehandler.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#include "observablehandler.h"
#include "datamapwriter.h"
#include "statistics.h"

ObservableHandler::ObservableHandler(const std::string& observableName,
		const MCParams& simulationParameters,
		const MetadataMap& metadataToStoreModel,
		const MetadataMap& metadataToStoreMC)
	: name(observableName), mcparams(simulationParameters),
	  metaModel(metadataToStoreModel), metaMC(metadataToStoreMC),
	  jkBlockCount(mcparams.jkBlocks),
	  jkBlockSizeSweeps(mcparams.sweeps / jkBlockCount),
	  lastSweepLogged(0), countValues(0),
	  jkBlockValues(jkBlockCount, 0),
	  total(0)
{
	if (mcparams.timeseries) {
		std::string filename = observableName + ".series";
		storage = std::unique_ptr<DoubleVectorWriterSuccessive>(
				new DoubleVectorWriterSuccessive(filename));
		storage->addHeaderText("Timeseries for observable " + observableName);
		storage->addMetadataMap(metaModel);
		storage->addMetadataMap(metaMC);
		storage->addMeta("observable", observableName);
		storage->writeHeader();
	}
}

ObservableHandler::~ObservableHandler() {
}

void ObservableHandler::insertValue(num value, unsigned curSweep) {
	if (mcparams.timeseries) {
		timeseriesBuffer.push_back(value);
	}

	unsigned curJkBlock = curSweep / jkBlockSizeSweeps;
	for (unsigned jb = 0; jb < jkBlockCount; ++jb) {
		if (jb != curJkBlock) {
			jkBlockValues[jb] += value;
		}
	}
	total += value;
	++countValues;
	lastSweepLogged = curSweep;
}

std::tuple<num, num> ObservableHandler::evaluateJackknife() const {
	num mean = total / countValues;
	num error = 0;
	if (lastSweepLogged == mcparams.sweeps) {
		//after the first sweep lastSweepLogged==1 and so on --> here the simulation is finished.
		//we can only calculate an error estimate if we have multiple jackknife blocks, or by
		//considering the naive variance estimate of the timeseries, the whole timeseries stored
		//in memory
		if (jkBlockCount <= 1) {
			if (timeseriesBuffer.size() == countValues) {
				error = variance(timeseriesBuffer, mean);
			} else {
				error = 0;
			}
		} else {
			unsigned jkBlockSizeSamples = countValues / jkBlockCount;
			unsigned jkTotalSamples = countValues - jkBlockSizeSamples;
			std::vector<num> jkBlockAverages = jkBlockValues;	//copy
			for (unsigned jb = 0; jb < jkBlockCount; ++jb) {
				jkBlockAverages[jb] /= jkTotalSamples;
			}
			error = jackknife(jkBlockAverages, mean);
		}
	}
	return std::make_tuple(mean, error);
}

void ObservableHandler::outputTimeseries() {
	//TODO: reserve reasonable amount of memory for data to be added afterwards
	//TODO: float precision
	if (mcparams.timeseries) {
		storage->writeData(timeseriesBuffer);	//append last batch of measurements
		timeseriesBuffer.resize(0);				//no need to keep it in memory anymore
	}
}

void outputResults(const std::vector<std::unique_ptr<ObservableHandler>>& obsHandlers) {
	typedef std::map<std::string, num> StringNumMap;
	typedef std::shared_ptr<StringNumMap> StringNumMapPtr;
	typedef DataMapWriter<std::string, num> StringNumWriter;

	StringNumMapPtr values(new StringNumMap);
	StringNumMapPtr errors(new StringNumMap);
	for (auto p = obsHandlers.cbegin(); p != obsHandlers.cend(); ++p) {
		num val, err;
		std::tie(val, err) = (*p)->evaluateJackknife();
		std::string obsname = (*p)->name;
		(*values)[obsname] = val;
		(*errors)[obsname] = err;
	}
	DataMapWriter<std::string, num> output;
	output.setData(values);
	output.setErrors(errors);
	output.addHeaderText("Monte Carlo results for observable expectation values");
	output.addMetadataMap((*obsHandlers.begin())->metaModel);
	output.addMetadataMap((*obsHandlers.begin())->metaMC);
	output.addMeta("key", "observable");
	output.addHeaderText("observable\t value \t error");
	output.writeToFile("results.values");
}




