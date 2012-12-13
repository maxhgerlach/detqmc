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
}

DetQMC::~DetQMC() {
}

