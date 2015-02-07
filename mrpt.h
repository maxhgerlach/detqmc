//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * multireweighthisto-pt.h
 *
 *  Created on: Jun 28, 2011
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )

#ifndef MULTIREWEIGHTHISTO_PT_H_
#define MULTIREWEIGHTHISTO_PT_H_


#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <boost/multi_array.hpp>
#include "exceptions.h"
#include "logval.h"
#include "histogram.h"
#include "reweightingresult.h"

//internally don't use the HistogramT class, which would use maps internally

//multiple histogram reweighting adapted for data from parallel tempering as in Chodera 2007
// [timeseries / per bin correlation handling optional]

// for classical (extended-) canonical MC simulations: uses energy and observable time series for different beta

// for SDW DQMC simulations: timeseries for different r: associated "energy" (0.5 sum phi^2), observable (e.g. |phi| magnetization)

class MultireweightHistosPT {

public:
    //Data kept public do make things simple. Handle with care!
    
    std::string controlParameterName; // name of control parameter, for classical (extended-) canonical MC simulations: inverse temperature beta

    std::string observable;                 //name of second observable (other than energy)
    int infoNumSamples;                     //number of samples, value taken from info.dat
    int infoSweepsBetweenMeasurements;      //taken from info.dat
    typedef std::vector<std::vector<double>*> DoubleSeriesCollection;       //map replica index -> pointer to time series of (floating point) measurements
    typedef std::vector<std::vector<int>*> IntSeriesCollection;             //map replica index -> pointer to series of integer values
    DoubleSeriesCollection energyTimeSeries;
    DoubleSeriesCollection observableTimeSeries;
    IntSeriesCollection cpiTimeSeries; //control parameter index time series for each replica, currently only used for the Fenwick-estimation of d.o.s. --> can be freed after histograms have been computed
    std::vector<double> controlParameterValues;
    unsigned numReplicas;
    unsigned systemN;       // should be systemL ** d
    unsigned systemL;
    unsigned systemSize;        // for classical MC == systemN; for QMC == systemN * beta / Delta_tau 

    //for energy histograms
    //energy-bins span: [U_m/N - binsize/2 ... U_m/N ... U_m/N + binsize/2]
    double minEnergyNormalized, maxEnergyNormalized;
    unsigned binCount;
    double binSize;     //Delta_U / N
    double deltaU;                      //not normalized
    double minEnergy, maxEnergy;        //not normalized
    LogVal lBinSize;
    std::vector<double> U_m;    //energy in the center of bin m (*not* normalized by system volume)
    typedef boost::multi_array<int, 2> Int2Array;
    Int2Array H_km;     //energy histogram: H_km[k][m] --> count of bin m from replica k
    typedef std::vector<int> Histo;
    Histo H_m;          //H_m = sum_k {H_km}
    typedef boost::multi_array<double, 2> Double2Array;
    Double2Array g_km;  //g_km[k][m] --> statistical inefficiency of bin m in replica k
    IntSeriesCollection m_kn;           //map replica index -> pointer to time series of occupied energy bins

    Int2Array H_lm;     //energy histogram: H_lm[l][m] --> count of bin m at inverse temperature l (accumulated from all replicas)

    Int2Array N_kl;             //counts of samples at temperature l in replica k

    std::vector<double> Heff_m;     // sum_k {(g_mk ^ -1) * H_mk
    std::vector<LogVal> lHeff_m;
    typedef boost::multi_array<LogVal, 2> LogVal2Array;
    Double2Array Neff_lm;           // sum_k {(g_mk ^ -1) * N_kl}
    LogVal2Array lNeff_lm;
    std::vector<double> Neff_l;     // sum_m {Neff_lm}
    LogVal2Array lPrecalc_lm;       // precalculated for updateDensityOfStates

    std::vector<LogVal> lZ_l;       //partition function at beta_l

    std::vector<LogVal> lOmega_m;   //estimate of density of states (up to multiplicative factor)

    //for observable histograms:
    double minObservableNormalized, maxObservableNormalized;

    std::ostream& out;


public:
    //Interface

    MultireweightHistosPT(std::ostream& outStream = std::cout);
    virtual ~MultireweightHistosPT();

    void addSimulationInfo(const std::string& filename);        //pass "info.dat"...
    unsigned getSystemN() const { return systemN; }
    unsigned getSystemL() const { return systemL; }
    std::string getObservableName() const { return observable; }
    std::string getControlParameterName() const { return controlParameterName; }

    // I) version for two-column timeseries [beta-index, obs-value],
    //    currently only supports control parameter beta,
    //    this is meant for the classical MC parallel tempering code.
    //distinguishes automatically by observable (denoted in meta data)
    //also extracts beta index timeseries
    //optionally subsamples the time series (subsample defines sample size)
    //optionally leaves out the first @discardEntries entries (extra thermalization)
    void addInputTimeSeries_twoColumn(const std::string& filename, unsigned subsample = 1,
                                      unsigned discardEntries = 0);

    // II) version for single-column timeseries, [obs-value] exclusively,
    //     which are already sorted by temperature
    //distinguishes automatically by observable (denoted in meta data)
    //optionally subsamples the time series (subsample defines sample size)
    //optionally leaves out the first @discardEntries entries (extra thermalization)
    void addInputTimeSeries_singleColumn(const std::string& filename, unsigned subsample = 1,
                                         unsigned discardEntries = 0);

    //compute averages directly at the original temperatures, only
    //from data at the original temperatures -- no reweighting
    //return a newly allocated map
    //    beta -> results
    typedef std::map<double, ReweightingResult> ResultsMap;
    virtual ResultsMap* directNoReweighting();

    //sort replica time series according to control parameter index, after this call
    //the timeseries for replica i only contain measurements at, e.g., temperature
    //beta_i.
    
    //This will probably hide correlations (the idea of the Chodera
    //article was to avoid this, we have found this not to matter in
    //the actual implementation).
    void sortTimeSeriesByControlParameter();

    //this should be done for the whole time series even when jackkniving
    //if saveAutocorr = true: store evaluated autocorr function points to (a) file(s)
    void setBinInefficienciesToUnity();     //most simple choice, results might be the worst -- actually little difference
    void measureBinInefficiencies(bool saveAutocorr = false);        //as in Chodera
    void measureGlobalInefficiencies(bool saveAutocorr = false);     //as in traditional Ferrenberg-Swendsen

    //only sensible if sorted by temperature:
    void writeOutEnergyTauInt(const std::string& filename);
    void writeOutObsTauInt(const std::string& filename);

    void updateEffectiveCounts();       //TODO: currently also called by findDensityOfStates

    virtual void createHistograms(int binCount_);
    virtual void createHistogramsIsing();      //special version for 2D Ising model (natural binning)

    //compute cross correlations of normalized histogram bin entries (H_km/N_k, H_kn/N_k)
    //write tables into sub-directories:
    void computeAndSaveHistogramCrossCorr();
    //use an alternative estimator (prob. less robust, but maybe less assumptions):
    void computeAndSaveHistogramCrossCorrAlt();

    void saveH_km(const std::string& filename);
    void saveg_km(const std::string& filename);
    void saveNeff_lm(const std::string& filename);
    void saveNeff_l(const std::string& filename);
    void savelNeff_lm(const std::string& filename);
    void saveU_m(const std::string& filename);

    virtual void findDensityOfStatesNonIteratively();       //using the method described in Fenwick, 2008 (uses histograms by temperature)
    virtual void findPartitionFunctionsAndDensityOfStates(double tolerance = 1E-7, int maxIterations = 1000);

    virtual void saveLogDensityOfStates(const std::string& filename);       //saves the logarithm. normalized to 0 for the central entry
    virtual void saveLogDensityOfStatesIsing(const std::string& filename);  //saves the logarithm. normalized to ln(2) for the lowest energy entry
    void savePartitionFunctions(const std::string& filename);
    void loadPartitionFunctions(const std::string& filename);       //afterwards call findPartitionFunctionsAndDensityOfStates again

    //FS-multi-reweighting without error estimation, all observables known
    //using time series:
    virtual ReweightingResult reweight(double targetControlParameter);

    //FS-multi-reweighting without error estimation, all observables known,
    // + energy and observable histograms at the target beta
    //the returned structure holds pointers to newly allocated memory
    // -- they should be deleted at some point
    //using time series:
    virtual ReweightingResult reweightWithHistograms(double targetControlParameter,
            unsigned obsBinCount);

    //FS-multi-reweighting without error estimation for one specific observable
    //using time series:
    virtual double reweightEnergy(double targetControlParameter);
    virtual double reweightSpecificHeat(double targetControlParameter);
    virtual double reweightObservable(double targetControlParameter);
    virtual double reweightObservableSusceptibility(double targetControlParameter);
    virtual double reweightObservableBinder(double targetControlParameter);

    //find the maximum of the susceptibility between control
    //parameters cpStart and cpEnd put its location into cpMax,
    //its value into cpMax and add all points evaluated to the map
    //pointsEvaluated (using time series)
    void findMaxObservableSusceptibility(double& cpMax, double& suscMax,
            std::map<double, double>& pointsEvaluated,
            double cpStart, double cpEnd);
    void findMaxSpecificHeatDiscrete(double& cpMax, double& suscMax,
            std::map<double, double>& pointsEvaluated,
            double cpStart, double cpEnd);
    void findMinBinder(double& cpMin, double& binderMin,
            std::map<double, double>& pointsEvaluated,
            double cpStart, double cpEnd);

    //Search between cpStart and cpEnd for a double-peak histogram
    //with as equal peak-heights as possible, obtain histograms from the
    //reweighting. .0<tolerance<1.0 should be chosen small enough to
    //distinguish the dip from the peaks and large enough to distinguish
    //the dip from noise...
    //Return the parameter location for this histogram in cpDouble
    //and the ratio max(max1,max2)/(min) in relativeDip. If no dip is found,
    //relativeDip is set to 1.0. A pointer to the final histogram is put
    //into histo.
    //The function for observable histograms has an additional parameter
    //specifying the number of bins to use for the reweighted histograms.
    void findEnergyEqualHeight(double& cpDouble, double& relativeDip,
            HistogramDouble*& histo,
            double cpStart, double cpEnd, double tolerance = 0.1);
    void findObsEqualHeight(double& cpDouble, double& relativeDip,
            HistogramDouble*& histo,
            double cpStart, double cpEnd, unsigned numBins,
            double tolerance = 0.1);
 
    //The following optimize for equal peak-*weight* (as in [JankeFirstOrder])
    //Additional input parameter: previously optimized equal *height*
    //histogram --> used to find cut-off
    void findEnergyEqualWeight(double& cpDouble, double& relativeDip,
            HistogramDouble*& histoResult,
            const HistogramDouble* equalHeightHisto,
            double cpStart, double cpEnd, double tolerance = 0.1);
    void findObsEqualWeight(double& cpDouble, double& relativeDip,
            HistogramDouble*& histoResult,
            const HistogramDouble* equalHeightHisto,
            double cpStart, double cpEnd, unsigned numBins,
            double tolerance = 0.1);

    //The following just reweight to obtain the histogram at the target cp,
    //return that and the relative dip
    void energyRelDip(double& relDip, HistogramDouble*& histoResult,
            double targetControlParameter, double tolerance = 0.1);
    void obsRelDip(double& relDip, HistogramDouble*& histoResult,
            double targetControlParameter, unsigned numBins, double tolerance = 0.1);


    //produces nice histograms for further processing / output at inverse
    //temperature @targetControlParameter. Returns a pointer to newly allocated memory
    //that has to be deleted at some point!
    HistogramDouble* reweightEnergyHistogramWithoutErrors(
            double targetControlParameter);
    HistogramDouble* reweightObservableHistogramWithoutErrors(
            double targetControlParameter, unsigned numBins);
    virtual HistogramDouble* reweightEnergyHistogram(double targetControlParameter) {
        return reweightEnergyHistogramWithoutErrors(targetControlParameter);
    }
    virtual HistogramDouble* reweightObservableHistogram(
            double targetControlParameter, unsigned numBins) {
        return reweightObservableHistogramWithoutErrors(targetControlParameter, numBins);
    }

    //obtain energy and specific heat estimates directly from the density of states
    //(do not re-consider the time series)
    virtual ReweightingResult reweightDiscrete(double beta);

protected:
    //interna

    bool basicConfig;

    unsigned addedEnergyTimeSeries;
    unsigned addedObservableTimeSeries;

    void createHistogramsHelper();
    void createHistogramsHelperDiscrete();
    void setUpHistograms(int binCount_);            //used in createHistograms (initializes data structures but does not set actual values)
    void setUpHistogramsIsing();

    void updateDensityOfStates();

    DoubleSeriesCollection computeWeights(double targetControlParameter);       //determine weights w_kn(cp) (occupies a lot of memory, free later)
    ReweightingResult reweightWithoutErrorsInternal(double targetControlParameter, const DoubleSeriesCollection& w_kn);

    //return the expectation value of the first moment of some observable
    //with given weights, put the result into firstMoment
    void reweight1stMomentInternalWithoutErrors(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn, double& firstMoment);
    //return the expectation value of the first and second moments of some
    //observable at a given temperature with given weights, put the result
    //into firstMoment, secondMoment
    void reweight1stMoment2ndMomentInternalWithoutErrors(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn,
        double& firstMoment, double& secondMoment);
    //likewise:
    void reweight2ndMoment4thMomentInternalWithoutErrors(
        const DoubleSeriesCollection& timeSeries,
        const DoubleSeriesCollection& w_kn,
        double& secondMoment, double& fourthMoment);


    HistogramDouble* reweightObservableHistogramUsingWeights(double targetControlParameter, unsigned numBins,
            const DoubleSeriesCollection& w_kn);

    //using pre-calculated weights @w_kn estimate a @numBins-histogram of the
    //observable at inverse temperature @targetControlParameter. Put the result into
    //output vector @histo
    void reweightObservableHistogramInternal(double targetControlParameter, unsigned numBins,
            const DoubleSeriesCollection& w_kn, std::vector<double>& histo);

};



#endif /* MULTIREWEIGHTHISTO_PT_H_ */
