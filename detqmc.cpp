/*
 * detqmc.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include <boost/assign/std/vector.hpp>
#include "detqmc.h"
#include "dethubbard.h"
#include "tools.h"
#include "git-revision.h"
#include "exceptions.h"



DetQMC::DetQMC(const ModelParams& parsmodel_, const MCParams& parsmc_) :
		parsmodel(parsmodel_), parsmc(parsmc_), sweepsDone(0)
{
	//TODO: RNG seed!
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

	if (parsmodel.model == "hubbard") {
		replica = createDetHubbard(parsmodel);
	}

	modelMeta = replica->prepareModelMetadataMap();
	mcMeta = prepareMCMetadataMap();
	for (unsigned obsIndex = 0; obsIndex < replica->getNumberOfObservables(); ++obsIndex) {
		obsHandlers.push_back(ObsPtr(
				new ObservableHandler(replica->getObservableName(obsIndex), parsmc,
						modelMeta, mcMeta)));
	}
}


DetQMC::~DetQMC() {
}


void DetQMC::run() {
	thermalize(parsmc.thermalization);
	for (unsigned sw = 0; sw < parsmc.sweeps; sw += parsmc.saveInterval) {
		measure(parsmc.saveInterval, parsmc.measureInterval);
		saveResults();
	}
}

void DetQMC::thermalize(unsigned numSweeps) {
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
		replica->sweepSimple();
	}
}

void DetQMC::measure(unsigned numSweeps, unsigned measureInterval) {
	for (unsigned sw = 0; sw < numSweeps; ++sw) {
		replica->sweepSimple();
		++sweepsDone;
		if (sw % measureInterval == 0) {
			//TODO: can this for loop be expressed in a more elegant and terse way?
			for (unsigned oi = 0; oi < replica->getNumberOfObservables(); ++oi) {
				obsHandlers[oi]->insertValue(replica->obsNormalized(oi), sweepsDone);
			}
		}
	}
}

MetadataMap DetQMC::prepareMCMetadataMap() {
	MetadataMap meta;
	meta["sweeps"] = numToString(parsmc.sweeps);
	meta["thermalization"] = numToString(parsmc.thermalization);
	meta["jkBlocks"] = numToString(parsmc.jkBlocks);
	meta["measureInterval"] = numToString(parsmc.measureInterval);
	meta["saveInterval"] = numToString(parsmc.saveInterval);
	meta["timeseries"] = (parsmc.timeseries ? "true" : "false");
	return meta;
}


void DetQMC::saveResults() {
	outputResults(obsHandlers);
	for (auto p = obsHandlers.begin(); p != obsHandlers.end(); ++p) {
		(*p)->outputTimeseries();
	}
	std::string commonInfoFilename = "info.dat";
	writeOnlyMetaData(commonInfoFilename, collectVersionInfo(),
			"Collected innformation about this determinantal quantum Monte Carlo simulation",
			true);
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
