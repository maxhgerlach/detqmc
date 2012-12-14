/*
 * observablehandler.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef OBSERVABLEHANDLER_H_
#define OBSERVABLEHANDLER_H_


// manage measurements of an observable
// calculate expectation values and jackknife error bars
// optionally store time series

#include <memory>
#include <string>
#include <vector>
#include <tuple>
#include "parameters.h"
#include "metadata.h"
#include "dataserieswritersucc.h"

class ObservableHandler {
public:
	ObservableHandler(const std::string& observableName,
			const MCParams& simulationParameters,
			const MetadataMap& metadataToStore);
	virtual ~ObservableHandler();

	// Log a newly measured observable value, pass the number of the current sweep.
	// Measurements do not need to be stored at every sweep, but the number of skipped
	// sweeps must be constant.
	void insertValue(num value, unsigned curSweep);

	//return [mean value, error] at end of simulation
	//return [mean value, 0] if this is called earlier
	std::tuple<num,num> evaluateJackknife() const;

	//update timeseries file, discard batch of
	//data written to file from memory
	void outputTimeseries();

	friend void ::outputResults(const std::vector<ObservableHandler>& obsHandlerso);
protected:
	std::string name;
	MCParams mcparams;
	MetadataMap meta;
	unsigned jkBlockCount;
	unsigned jkBlockSizeSweeps;

	unsigned lastSweepLogged;
	unsigned countValues;

	std::vector<num> jkBlockValues;			// running counts of jackknife block values
	num total;								// running accumulation regardless of jackknife block
	std::vector<num> timeseriesBuffer;		// time series entries added since last call to writeData()

	std::unique_ptr<DoubleVectorWriterSuccessive> storage;
};


//Write expectation values and error bars for all observables to a file
//take metadata to store from the first entry in obsHandlers
void outputResults(const std::vector<ObservableHandler>& obsHandlers);


#endif /* OBSERVABLEHANDLER_H_ */
