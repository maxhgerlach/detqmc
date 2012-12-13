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

#include <string>
#include <vector>
#include <tuple>
#include "parameters.h"

class ObservableHandler {
public:
	ObservableHandler(const std::string& observableName, const MCParams& simulationParameters);
	virtual ~ObservableHandler();

	// Log a newly measured observable value, pass the number of the current sweep.
	// Measurements do not need to be stored at every sweep, but the number of skipped
	// sweeps must be constant.
	void insertValue(num value, unsigned curSweep);

	//return [mean value, error] at end of simulation
	//return [mean value, 0] if this is called to early
	std::tuple<num,num> evaluateJackknife();

protected:
	std::string name;
	MCParams mcparams;
	unsigned jkBlockCount;
	unsigned jkBlockSizeSweeps;

	unsigned lastSweepLogged;
	unsigned countValues;

	std::vector<num> jkBlockValues;			// running counts of jackknife block values, after evaluateJackknife: block averages
	num total;								// running accumulation regardles of jackknife block
};

#endif /* OBSERVABLEHANDLER_H_ */
