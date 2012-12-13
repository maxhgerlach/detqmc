/*
 * observablehandler.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#include "observablehandler.h"
#include "statistics.h"

ObservableHandler::ObservableHandler(const std::string& observableName,
		const MCParams& simulationParameters)
	: name(observableName), mcparams(simulationParameters),
	  jkBlockCount(mcparams.jkBlocks),
	  jkBlockSizeSweeps(mcparams.sweeps / jkBlockCount),
	  lastSweepLogged(0), countValues(0),
	  jkBlockValues(jkBlockCount, 0),
	  total(0)
{ }

ObservableHandler::~ObservableHandler() {
}

void ObservableHandler::insertValue(num value, unsigned curSweep) {
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

std::tuple<num, num> ObservableHandler::evaluateJackknife() {
	num mean = total / countValues;
	num error = 0;
	if (lastSweepLogged == mcparams.sweeps) {
		//after the first sweep lastSweepLogged==1 and so on --> here the simulation is finished
		//perform jackknife analysis
		unsigned jkBlockSizeSamples = countValues / jkBlockCount;
		unsigned jkTotalSamples = countValues - jkBlockSizeSamples;
		if (jkBlockCount == 1) {
			jkTotalSamples = countValues;       //would be zero otherwise
		}
		for (unsigned jb = 0; jb < jkBlockCount; ++jb) {
			jkBlockValues[jb] /= jkTotalSamples;
		}
		error = jackknife(jkBlockValues, mean);
	}
	return std::make_tuple(mean, error);
}


