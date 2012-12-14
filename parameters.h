/*
 * parameters.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>

//Collect various structs defining various parameters.
//Passing around these reduces code duplication somewhat, and reduces errors caused by passing
//values to the wrong (positional) constructor argument.

typedef double num;		//possibility to switch to single precision if ever desired

// Struct representing model parameters
struct ModelParams {
	std::string model;
	num t;
	num U;
	num mu;
	unsigned L;
	unsigned d;
	num beta;
	unsigned m;
};


// Struct representing Monte Carlo simulation parameters
struct MCParams {
	unsigned sweeps;			// number of sweeps used for measurements
	unsigned thermalization;	// number of warm-up sweeps allowed before equilibrium is assumed
	unsigned jkBlocks;			// number of jackknife blocks for error estimation
	bool timeseries;			// if true, write time series of individual measurements to disk
};



#endif /* PARAMETERS_H_ */
