//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * reweightingresult.h
 *
 *  Created on: Jun 29, 2011
 *      Author: gerlach
 */


// generalized for SDW DQMC (2015-02-06 - )


#ifndef REWEIGHTINGRESULT_H_
#define REWEIGHTINGRESULT_H_

#include "histogram.h"
#include "tools.h"

//this structure holds mean values, susceptibility, binder parameter and their errors
//determined at a certain beta by multi-reweighting
struct ReweightingResult {
    double energyAvg;           // for SDW: r/2 sum {phi^2}
    double energyError;
    double heatCapacity;        // will not be sensible for SDW
    double heatCapacityError;
    double obsAvg;
    double obsError;
    double obsSquared;          // systemSize * (<o^2>)
    doube  obsSquaredError;
    double obsSusc;             // systemSize * (<o^2> - <o>^2)
    double obsSuscError;
    double obsBinder;           // 1 - <o^4> / (3 <o^2>^2)
    double obsBinderError;
    double obsBinderRatio;      // <o^4> / <o^2>^2
    double obsBinderRatioError;


    //optional pointers to reweighted histograms
    HistogramDouble* energyHistogram;
    HistogramDouble* obsHistogram;

    //default constructor: all zero
    ReweightingResult() :
        energyAvg(0), energyError(0), heatCapacity(0), heatCapacityError(0),
        obsAvg(0), obsError(0), obsSquared(0), obsSquaredError(0),
        obsSusc(0), obsSuscError(0), obsBinder(0),
        obsBinderError(0), obsBinderRatio(0),
        obsBinderRatioError(0),
        energyHistogram(0), obsHistogram(0)
    { }

    //constructor setting only mean values, not error margins
    // ReweightingResult(double energy, double obs, double susc, double binder,
    //     ) :
    //     energyAvg(energy), energyError(-1), heatCapacity(-1),
    //             heatCapacityError(-1), obsAvg(obs), obsError(-1),
    //             obsSusc(susc), obsSuscError(-1), obsBinder(binder),
    //             obsBinderError(-1), energyHistogram(0), obsHistogram(0)
    // { }
    ReweightingResult(double energy, double heatCapacity, double obs,
                      double obsSquared, 
                      double susc, double binder, double binderRatio) :
        energyAvg(energy), energyError(-1), heatCapacity(heatCapacity),
        heatCapacityError(-1), obsAvg(obs), obsError(-1),
        obsSquared(obsSquared), obsSquaredError(-1),        
        obsSusc(susc), obsSuscError(-1), obsBinder(binder),
        obsBinderError(-1), obsBinderRatio(binderRatio),
        obsBinderRatioError(-1),
        energyHistogram(0), obsHistogram(0)
    { }

    //constructor setting values and errors
    // ReweightingResult(double energy, double energyErr, double obs,
    //         double obsErr, double susc, double suscErr, double binder,
    //         double binderErr) :
    //     energyAvg(energy), energyError(energyErr), heatCapacity(-1),
    //             heatCapacityError(-1), obsAvg(obs), obsError(obsErr),
    //             obsSusc(susc), obsSuscError(suscErr), obsBinder(binder),
    //             obsBinderError(binderErr), energyHistogram(0), obsHistogram(0)
    // { }
    ReweightingResult(double energy, double energyErr, double heatCapacity,
                      double heatCapacityError, double obs, double obsErr,
                      double obsSquared, double obsSquaredErr,
                      double susc, double suscErr,
                      double binder, double binderErr,
                      double binderRatio, double binderRatioErr) :
        energyAvg(energy), energyError(energyErr), heatCapacity(heatCapacity),
        heatCapacityError(heatCapacityError), obsAvg(obs),
        obsError(obsErr), obsSquared(obsSquared), obsSquaredError(obsSquaredErr),
        obsSusc(susc), obsSuscError(suscErr),
        obsBinder(binder), obsBinderError(binderErr),
        obsBinderRatio(binderRatio), obsBinderRatioError(binderRatioErr),        
        energyHistogram(0), obsHistogram(0)
    { }

    void freeMemory() {
        destroy(energyHistogram);
        destroy(obsHistogram);
    }

};


#endif /* REWEIGHTINGRESULT_H_ */
