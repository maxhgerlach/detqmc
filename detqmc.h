/*
 * detqmc.h
 *
 * Handling of determinantal quantum Monte Carlo simulations
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */

#ifndef DETQMC_H_
#define DETQMC_H_

#include <vector>
#include <memory>
#include "parameters.h"
#include "dethubbard.h"
#include "observablehandler.h"



// Class handling the simulation
class DetQMC {
public:
	DetQMC(const ModelParams& parsmodel, const MCParams& parsmc);
	virtual ~DetQMC();
protected:
	ModelParams parsmodel;
	MCParams parsmc;
	std::unique_ptr<DetHubbard> replica;		//extend to allow for derived classes in the future
	std::vector<std::unique_ptr<ObservableHandler>> obsHandlers;
};





#endif /* DETQMC_H_ */
