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

#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"

// Collect various structs defining various parameters.
// Passing around these reduces code duplication somewhat, and reduces errors caused by passing
// values to the wrong (positional) constructor argument.

// The set specified included in each struct contains string representations
// of all parameters actually specified.  This allows throwing an exception
// at the appropriate point in the program if a parameter is missing.
//TODO: Find a more elegant solution.

typedef double num;     //possibility to switch to single precision if ever desired

// Struct representing model specific parameters
struct ModelParams {
    //TODO: split this up so that model specific parameters are in substructs

    std::string model;
    bool timedisplaced;     //also evaluate timedisplaced Green functions and derived quantities
    bool checkerboard;      //SDW: use a checkerboard decomposition for computing the propagator
    std::string checkerboardMethod;		//SDW if checkerboard: "santos" or "assaad" or "assaad_berg"
    num t;      //Hubbard
    num U;      //Hubbard
    num r;      //SDW
    num txhor;  //SDW: hopping constants depending on direction and band
    num txver;  //SDW
    num tyhor;  //SDW
    num tyver;  //SDW
    num mu;
    uint32_t L;
    uint32_t d;
    num beta;
    uint32_t m;     //either specify number of timeslices 'm'
    num dtau;       //or timeslice separation 'dtau'
    uint32_t s;     //separation of timeslices where the Green function is calculated
                    //from scratch
    num accRatio;   //for SDW: target acceptance ratio for tuning spin update box size

    std::string bc; //boundary conditions: For SDW: "pbc", "apbc-x", "apbc-y" or "apbc-xy"

    bool rescale;	//perform global rescale move?
    uint32_t rescaleInterval;		//attempt global rescale move every # sweeps
    num rescaleGrowthFactor;		// factor by which to size up the fields
    num rescaleShrinkFactor;		// factor by which to size down the fields

    std::set<std::string> specified;

    ModelParams() :
            model(), timedisplaced(), checkerboard(), checkerboardMethod(), t(), U(), r(),
            txhor(), txver(), tyhor(), tyver(),
            mu(), L(), d(),
            beta(), m(), dtau(), s(), accRatio(), bc("pbc"),
            rescale(), rescaleInterval(), rescaleGrowthFactor(), rescaleShrinkFactor(),
            specified() {
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar & model & timedisplaced & checkerboard
           & t & U & r & txhor & txver & tyhor & tyver
           & mu & L & d & beta & m & dtau & s & accRatio & bc
           & rescale & rescaleInterval & rescaleGrowthFactor & rescaleShrinkFactor
           & specified;
    }
};


// Struct representing Monte Carlo simulation parameters
struct MCParams {
    uint32_t sweeps;            // number of sweeps used for measurements
    uint32_t thermalization;    // number of warm-up sweeps allowed before equilibrium is assumed
    uint32_t jkBlocks;          // number of jackknife blocks for error estimation
    bool timeseries;            // if true, write time series of individual measurements to disk
    uint32_t measureInterval;   // take measurements every measureInterval sweeps
    uint32_t saveInterval;      // write measurements to disk every saveInterval sweeps
    uint32_t rngSeed;       // seed for random number generator

    std::string greenUpdateType;    //"simple" or "stabilized"

    std::string stateFileName;      //for serialization dumps
    bool sweepsHasChanged;          //true, if the number of target sweeps has changed after resuming

    std::set<std::string> specified;

    MCParams() : sweeps(), thermalization(), jkBlocks(), timeseries(false), measureInterval(), saveInterval(),
            rngSeed(), greenUpdateType(), stateFileName(), sweepsHasChanged(false), specified()
    { }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar & sweeps & thermalization & jkBlocks & timeseries
           & measureInterval & saveInterval & rngSeed & greenUpdateType
           & stateFileName
           & sweepsHasChanged
           & specified;
    }
};



#endif /* PARAMETERS_H_ */
