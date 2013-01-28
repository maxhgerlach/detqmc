/*
 * parameters.h
 *
 *  Created on: Dec 13, 2012
 *      Author: gerlach
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <string>
#include <set>

// Collect various structs defining various parameters.
// Passing around these reduces code duplication somewhat, and reduces errors caused by passing
// values to the wrong (positional) constructor argument.

// The set specified included in each struct contains string representations
// of all parameters actually specified.  This allows throwing an exception
// at the appropriate point in the program if a parameter is missing.
//TODO: Find a more elegant solution.

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
	unsigned m;		//either specify number of timeslices 'm'
	num dtau;			//or timeslice separation 'dtau'
	std::set<std::string> specified;
};


// Struct representing Monte Carlo simulation parameters
struct MCParams {
	unsigned sweeps;			// number of sweeps used for measurements
	unsigned thermalization;	// number of warm-up sweeps allowed before equilibrium is assumed
	unsigned jkBlocks;			// number of jackknife blocks for error estimation
	bool timeseries;			// if true, write time series of individual measurements to disk
	unsigned measureInterval;	// take measurements every measureInterval sweeps
	unsigned saveInterval;		// write measurements to disk every saveInterval sweeps
	unsigned long rngSeed;		// seed for random number generator

	std::string greenUpdateType; 	//"simple" or "stabilized"

	std::set<std::string> specified;
};



#endif /* PARAMETERS_H_ */
