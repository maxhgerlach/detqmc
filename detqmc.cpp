/*
 * detqmc.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include <ctime>
#include <functional>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include "boost/assign/std/vector.hpp"
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"
#include "detqmc.h"
#include "dethubbard.h"
#include "detsdw.h"
#include "tools.h"
#include "git-revision.h"
#include "exceptions.h"

using std::cout;
using std::endl;

DetQMC::DetQMC(const ModelParams& parsmodel_, const MCParams& parsmc_) :
		parsmodel(parsmodel_), parsmc(parsmc_),
		//proper initialization of default initialized members done in method body
		greenUpdateType(), sweepFunc(), sweepThermalizationFunc(),
		modelMeta(), mcMeta(), rng(), replica(),
		obsHandlers(), vecObsHandlers(),
		sweepsDone(0)
{
	//check parameters
	if (parsmodel.specified.count("model") == 0) {
		throw ParameterMissing("model");
	}
	using namespace boost::assign;
	std::vector<std::string> neededMCPars;
	neededMCPars += "sweeps", "thermalization", "jkBlocks", "timeseries", "measureInterval";
	for (auto p = neededMCPars.cbegin(); p != neededMCPars.cend(); ++p) {
		if (parsmc.specified.count(*p) == 0) {
			throw ParameterMissing(*p);
		}
	}
	if (parsmc.specified.count("saveInterval") == 0) {
		parsmc.saveInterval = parsmc.sweeps;		//only save at end
	}

	if (parsmc.specified.count("rngSeed") == 0) {
		cout << "No rng seed specified, will use std::time(0)" << endl;
		parsmc.rngSeed = (long unsigned) std::time(0);
	}
	rng = RngWrapper(parsmc.rngSeed);

	if (parsmodel.model == "hubbard") {
		replica = createDetHubbard(rng, parsmodel);
	} else if (parsmodel.model == "sdw") {
		replica = createDetSDW(rng, parsmodel);
	}
	else {
		throw ParameterWrong("model", parsmodel.model);
	}

	if (parsmc.greenUpdateType == "simple") {
		greenUpdateType = GreenUpdateType::Simple;
	} else if (parsmc.greenUpdateType == "stabilized") {
		greenUpdateType = GreenUpdateType::Stabilized;
	} else {
		throw ParameterWrong("greenUpdateType", parsmc.greenUpdateType);
	}

	if (greenUpdateType == GreenUpdateType::Simple) {
		sweepFunc = [this]() {replica->sweepSimple();};
		sweepThermalizationFunc= [this]() {replica->sweepSimpleThermalization();};
	} else if (greenUpdateType == GreenUpdateType::Stabilized) {
		sweepFunc = [this]() {replica->sweep();};
		sweepThermalizationFunc = [this]() {replica->sweepThermalization();};
	}

	//some parameter consistency checking:
	if (parsmc.sweeps % parsmc.jkBlocks != 0) {
		throw ParameterWrong("Number of jackknife blocks " + numToString(parsmc.jkBlocks)
				+ " does not match number of sweeps " + numToString(parsmc.sweeps));
	}
	if ((parsmc.measureInterval > parsmc.sweeps) or
		(parsmc.sweeps % parsmc.measureInterval != 0)) {
		throw ParameterWrong("Measurement interval " + numToString(parsmc.measureInterval)
				+ " ill-chosen for number of sweeps " + numToString(parsmc.sweeps));
	}

	//prepare metadata
	modelMeta = replica->prepareModelMetadataMap();
	mcMeta = prepareMCMetadataMap();

	//prepare observable handlers
	auto scalarObs = replica->getScalarObservables();
	for (auto obsP = scalarObs.cbegin(); obsP != scalarObs.cend(); ++obsP) {
		obsHandlers.push_back(ObsPtr(
				new ScalarObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
	}
	auto vectorObs = replica->getVectorObservables();
	for (auto obsP = vectorObs.cbegin(); obsP != vectorObs.cend(); ++obsP) {
		vecObsHandlers.push_back(VecObsPtr(
				new VectorObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
	}
	auto keyValueObs = replica->getKeyValueObservables();
	for (auto obsP = keyValueObs.cbegin(); obsP != keyValueObs.cend(); ++obsP) {
		vecObsHandlers.push_back(VecObsPtr(
				new KeyValueObservableHandler(*obsP, parsmc, modelMeta, mcMeta)));
	}

	cout << "\nSimulation initialized, parameters: " << endl;
	cout << metadataToString(mcMeta, " ") << metadataToString(modelMeta, " ") << endl;
}


DetQMC::~DetQMC() {
}


void DetQMC::run() {
	thermalize(parsmc.thermalization);
	cout << "Starting measurements for " << parsmc.sweeps << " sweeps..." << endl;
	for (unsigned sw = 0; sw < parsmc.sweeps; sw += parsmc.saveInterval) {
		measure(parsmc.saveInterval, parsmc.measureInterval);
		cout << "  " << sw + parsmc.saveInterval << " ... saving results ...";
		saveResults();
		cout << endl;
	}
}

void DetQMC::thermalize(unsigned numSweeps) {
	cout << "Thermalization for " << numSweeps << " sweeps..." << endl;
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
		std::cout << sw << '\n';
		sweepFunc();
	}
	cout << endl;
}

void DetQMC::measure(unsigned numSweeps, unsigned measureInterval) {
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
		std::cout << sw << '\n';
		sweepFunc();
		++sweepsDone;
		if (sw % measureInterval == 0) {
			replica->measure();
			//TODO: can this for loop be expressed in a more elegant and terse way?
			for (auto ph = obsHandlers.begin(); ph != obsHandlers.end(); ++ph) {
				(*ph)->insertValue(sweepsDone);
			}
			for (auto ph = vecObsHandlers.begin(); ph != vecObsHandlers.end(); ++ph) {
				(*ph)->insertValue(sweepsDone);
			}
		}
	}
}


MetadataMap DetQMC::prepareMCMetadataMap() const {
	MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(parsmc.VAR)
	META_INSERT(greenUpdateType);
	META_INSERT(sweeps);
	META_INSERT(thermalization);
	META_INSERT(jkBlocks);
	META_INSERT(measureInterval);
	META_INSERT(saveInterval);
	META_INSERT(rngSeed);
#undef META_INSERT
	meta["timeseries"] = (parsmc.timeseries ? "true" : "false");
	return meta;
}


void DetQMC::saveResults() {
	outputResults(obsHandlers);
	for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
		(*p)->outputTimeseries();
	}
	outputResults(vecObsHandlers);
	std::string commonInfoFilename = "info.dat";
	writeOnlyMetaData(commonInfoFilename, collectVersionInfo(),
			"Collected innformation about this determinantal quantum Monte Carlo simulation",
			false);
	writeOnlyMetaData(commonInfoFilename, modelMeta,
			"Model parameters:",
			true);
	writeOnlyMetaData(commonInfoFilename, mcMeta,
			"Monte Carlo parameters:",
			true);
	MetadataMap currentState;
	currentState["sweepsDone"] = numToString(sweepsDone);
	writeOnlyMetaData(commonInfoFilename, currentState,
			"Current state of simulation:",
			true);
}
