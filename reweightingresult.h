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
    double energyAvg;
    double energyError;
    double heatCapacity;
    double heatCapacityError;
    double obsAvg;
    double obsError;
    double obsSusc;
    double obsSuscError;
    double obsBinder;
    double obsBinderError;

    //optional pointers to reweighted histograms
    HistogramDouble* energyHistogram;
    HistogramDouble* obsHistogram;

    //default constructor: all zero
    ReweightingResult() :
        energyAvg(0), energyError(0), heatCapacity(0), heatCapacityError(0),
        obsAvg(0), obsError(0), obsSusc(0), obsSuscError(0), obsBinder(0),
        obsBinderError(0), energyHistogram(0), obsHistogram(0)
    { }

    //constructor setting only mean values, not error margins
    ReweightingResult(double energy, double obs, double susc, double binder) :
        energyAvg(energy), energyError(-1), heatCapacity(-1),
                heatCapacityError(-1), obsAvg(obs), obsError(-1),
                obsSusc(susc), obsSuscError(-1), obsBinder(binder),
                obsBinderError(-1), energyHistogram(0), obsHistogram(0)
    { }
    ReweightingResult(double energy, double heatCapacity, double obs,
            double susc, double binder) :
        energyAvg(energy), energyError(-1), heatCapacity(heatCapacity),
                heatCapacityError(-1), obsAvg(obs), obsError(-1),
                obsSusc(susc), obsSuscError(-1), obsBinder(binder),
                obsBinderError(-1), energyHistogram(0), obsHistogram(0)
    { }

    //constructor setting values and errors
    ReweightingResult(double energy, double energyErr, double obs,
            double obsErr, double susc, double suscErr, double binder,
            double binderErr) :
        energyAvg(energy), energyError(energyErr), heatCapacity(-1),
                heatCapacityError(-1), obsAvg(obs), obsError(obsErr),
                obsSusc(susc), obsSuscError(suscErr), obsBinder(binder),
                obsBinderError(binderErr), energyHistogram(0), obsHistogram(0)
    { }
    ReweightingResult(double energy, double energyErr, double heatCapacity,
            double heatCapacityError, double obs, double obsErr, double susc,
            double suscErr, double binder, double binderErr) :
        energyAvg(energy), energyError(energyErr), heatCapacity(heatCapacity),
                heatCapacityError(heatCapacityError), obsAvg(obs),
                obsError(obsErr), obsSusc(susc), obsSuscError(suscErr),
                obsBinder(binder), obsBinderError(binderErr),
                energyHistogram(0), obsHistogram(0)
    { }

    void freeMemory() {
        destroy(energyHistogram);
        destroy(obsHistogram);
    }

};


#endif /* REWEIGHTINGRESULT_H_ */
