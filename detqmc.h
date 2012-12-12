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

#include "dethubbard.h"

struct Params;       //definition given below in this file

// Class handling the simulation
class DetQMC {
public:
	//init from command line
	DetQMC(const Params& par);
	virtual ~DetQMC();
protected:
};

// Struct representing simulation parameters
struct Params {
	std::string model;
	num t;
	num U;
	num mu;
	unsigned L;
	unsigned d;
	num beta;
	unsigned m;
};



#endif /* DETQMC_H_ */
