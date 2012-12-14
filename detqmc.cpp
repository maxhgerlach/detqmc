/*
 * detqmc.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include "detqmc.h"
#include "dethubbard.h"



DetQMC::DetQMC(const ModelParams& parsmodel, const MCParams& parsmc) :
		parsmodel(parsmodel), parsmc(parsmc)
{
	if (parsmodel.model == "hubbard") {
		replica = createDetHubbard(parsmodel);
	}

	modelMeta = replica->prepareModelMetadataMap();
	for (unsigned obsIndex = 0; obsIndex < replica->getNumberOfObservables(); ++obsIndex) {
		obsHandlers.push_back(ObsPtr(
				new ObservableHandler(replica->getObservableName(obsIndex), parsmc, modelMeta)));
	}
}

DetQMC::~DetQMC() {
}

