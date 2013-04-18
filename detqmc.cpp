/*
 * detqmc.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include <ctime>
#include <functional>
#include <fstream>
#include <array>
#include <armadillo>
#include "boost/assign/std/vector.hpp"

#include "boost_serialize_armadillo.h"
#include "boost_serialize_array.h"


#include "boost/archive/binary_oarchive.hpp"
#include "boost/archive/binary_iarchive.hpp"


#include "detqmc.h"
#include "dethubbard.h"
#include "detsdw.h"
#include "tools.h"
#include "git-revision.h"
#include "exceptions.h"
#include "timing.h"

using std::cout;
using std::endl;

void DetQMC::initFromParameters(const ModelParams& parsmodel_, const MCParams& parsmc_) {
	parsmodel = parsmodel_;
	parsmc = parsmc_;
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

DetQMC::DetQMC(const ModelParams& parsmodel_, const MCParams& parsmc_) :
		parsmodel(), parsmc(),
		//proper initialization of default initialized members done in initFromParameters
		greenUpdateType(), sweepFunc(), sweepThermalizationFunc(),
		modelMeta(), mcMeta(), rng(), replica(),
		obsHandlers(), vecObsHandlers(),
		sweepsDone(0), sweepsDoneThermalization()
{
	initFromParameters(parsmodel_, parsmc_);
}

DetQMC::DetQMC(const std::string& stateFileName, unsigned newSweeps) :
		parsmodel(), parsmc(),
		//proper initialization of default initialized members done by loading from archive
		greenUpdateType(), sweepFunc(), sweepThermalizationFunc(),
		modelMeta(), mcMeta(), rng(), replica(),
		obsHandlers(), vecObsHandlers(),
		sweepsDone(), sweepsDoneThermalization()
{
	std::ifstream ifs;
	ifs.exceptions(std::ifstream::badbit | std::ifstream::failbit);
	ifs.open(stateFileName.c_str(), std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	ModelParams parsmodel_;
	MCParams parsmc_;
	ia >> parsmodel_ >> parsmc_;
	if (newSweeps > parsmc_.sweeps) {
		parsmc_.sweeps = newSweeps;
	}
	parsmc_.stateFileName = stateFileName;
	//TODO: changing the number of target sweeps makes the on the fly Jackknife
	//error-estimation invalid
	initFromParameters(parsmodel_, parsmc_);
	serializeContents(ia);
}

void DetQMC::saveState() {
	timing.start("saveState");
	std::ofstream ofs;
	ofs.exceptions(std::ofstream::badbit | std::ofstream::failbit);
	ofs.open(parsmc.stateFileName.c_str(), std::ios::binary);
	boost::archive::binary_oarchive oa(ofs);
	oa << parsmodel << parsmc;
	serializeContents(oa);
	timing.stop("saveState");
}

DetQMC::~DetQMC() {
}


void DetQMC::run() {
	if (sweepsDoneThermalization < parsmc.thermalization) {
		cout << "Thermalization for " << parsmc.thermalization << " sweeps..." << endl;
		for (unsigned sw = sweepsDoneThermalization;
				sw < parsmc.thermalization; sw += parsmc.saveInterval) {
			thermalize(parsmc.saveInterval);
			cout  << "  " << sw + parsmc.saveInterval << " ... saving state...";
			saveState();
			cout << endl;
		}
	}
	cout << "Thermalization finished\n" << endl;
	replica->thermalizationOver();
	if (sweepsDone < parsmc.sweeps) {
		cout << "Measurements for " << parsmc.sweeps << " sweeps..." << endl;
		for (unsigned sw = sweepsDone; sw < parsmc.sweeps; sw += parsmc.saveInterval) {
			measure(parsmc.saveInterval, parsmc.measureInterval);
			cout << "  " << sw + parsmc.saveInterval << " ... saving results and state ...";
			saveResults();
			saveState();
			cout << endl;
		}
	}
	cout << "Measurements finished\n" << endl;
}

void DetQMC::thermalize(unsigned numSweeps) {
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
		sweepThermalizationFunc();
		++sweepsDoneThermalization;
	}
}

void DetQMC::measure(unsigned numSweeps, unsigned measureInterval) {
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
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
	timing.start("saveResults");
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
	currentState["sweepsDoneThermalization"] = numToString(sweepsDoneThermalization);
	currentState["sweepsDone"] = numToString(sweepsDone);
	writeOnlyMetaData(commonInfoFilename, currentState,
			"Current state of simulation:",
			true);
	timing.stop("saveResults");
}
