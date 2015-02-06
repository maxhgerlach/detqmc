//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * multireweighthistosptjk.h
 *
 *  Created on: Jul 25, 2011
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )

#ifndef MULTIREWEIGHTHISTOSPTJK_H_
#define MULTIREWEIGHTHISTOSPTJK_H_

#include "mrpt.h"

//extends MultiReweightHistosPT with jackknifing for error estimation

class MultireweightHistosPTJK: public MultireweightHistosPT {
private:
    unsigned blockCount;

    //the data structures from the base class remain and still correspond to the estimates
    //done over the whole data set.

    //the following are related to data sets with one jackknife block left out

    typedef boost::multi_array<int, 3> Int3Array;
    Int3Array H_bkm;        //energy histogram: H_bkm[b][k][m] --> count of bin m from replica k, jackknife block b
    Int2Array H_bm;         //H_bm = sum_k {H_bkm}
    typedef boost::multi_array<double, 3> Double3Array;
    Double3Array g_bkm;     //g_bkm[b][k][m] --> statistical inefficiency of bin m in replica k, jackknife block b

    Int3Array H_blm;        //energy histogram: H_blm[b][l][m] --> count of bin m at inverse temperature l (accumulated from all replicas), jackknife block b

    Int3Array N_bkl;                //counts of samples at temperature l in replica k, jackknife block b

    Double2Array Heff_bm;       // sum_k {(g_mk ^ -1) * H_bmk
    LogVal2Array lHeff_bm;
    typedef boost::multi_array<LogVal, 3> LogVal3Array;
    Double3Array Neff_blm;          // sum_k {(g_mk ^ -1) * N_bkl}
    LogVal3Array lNeff_blm;
    LogVal3Array lPrecalc_blm;      // precalculated for updateDensityOfStates

    LogVal2Array lZ_bl;     //partition function at beta_l, jackknife block b

    LogVal2Array lOmega_bm; //estimate of density of states (up to multiplicative factor), jackknife block b

    void createHistogramsHelperJK();
    void createHistogramsHelperDiscreteJK();
    void setUpHistogramsHelperJK();
    void setUpHistogramsJK(int binCount_);
    void setUpHistogramsIsingJK();

    DoubleSeriesCollection computeWeightsJK(double targetBeta, unsigned jkBlock);
    void updateEffectiveCountsJK();                 //for all blocks
    void updateDensityOfStatesJK(unsigned block);   //for the indicated block

    ReweightingResult reweightJackknifeInternal(double targetBeta, unsigned jkBlock, const DoubleSeriesCollection& w_kn);

    //FS-multi-reweighting for one specific observable, one jackknife block
    void reweight1stMomentInternalJK(const DoubleSeriesCollection& timeSeries,
            const DoubleSeriesCollection& w_kn, unsigned jkBlock,
            double& firstMoment);
    void reweight1stMoment2ndMomentInternalJK(
            const DoubleSeriesCollection& timeSeries,
            const DoubleSeriesCollection& w_kn, unsigned jkBlock,
            double& firstMoment, double& secondMoment);
    void reweight2ndMoment4thMomentInternalJK(
            const DoubleSeriesCollection& timeSeries,
            const DoubleSeriesCollection& w_kn, unsigned jkBlock,
            double& secondMoment, double& fourthMoment);
    double reweightEnergyJK(double targetBeta, unsigned jkBlock);
    double reweightSpecificHeatJK(double targetBeta, unsigned jkBlock);
    double reweightSpecificHeatDiscreteJK(double targetBeta, unsigned jkBlock);
    double reweightObservableJK(double targetBeta, unsigned jkBlock);
    double reweightObservableSusceptibilityJK(double targetBeta,
            unsigned jkBlock);
    double reweightObservableBinderJK(double targetBeta, unsigned jkBlock);

    //reweight histograms taking only data from one of the jackknife blocks
    HistogramDouble* reweightEnergyHistogramJK(double targetBeta,
            unsigned jkBlock);
    HistogramDouble* reweightObservableHistogramJK(double targetBeta,
            unsigned numBins, unsigned jkBlock);

    friend class SuscMinCallableJK;
    friend class SpecificHeatDiscreteMinCallableJK;
    friend class EnergyHistogramPeakDiffMinCallableJK;
    friend class ObsHistogramPeakDiffMinCallableJK;
    friend class EnergyHistogramWeightDiffMinCallableJK;
    friend class ObsHistogramWeightDiffMinCallableJK;
    friend class BinderMinCallableJK;
    friend class BinderDiffMinCallableJK;
public:
    MultireweightHistosPTJK(unsigned jkTotalBlocks, std::ostream& outStream = std::cout);
    virtual ~MultireweightHistosPTJK();

    virtual void createHistograms(int binCount_);
    virtual void createHistogramsIsing();      //special version for 2D Ising model (natural binning)

    //compute averages directly at the original temperatures, only
    //from data at the original temperatures -- no reweighting
    //return a newly allocated map
    //    beta -> results
    //here: additional Jackknife error estimation
    virtual ResultsMap* directNoReweighting();

    //estimate jackknife errors of H_km, write them out as simple table
    //("companion" to saveH_km)
    void saveH_km_errors(const std::string& filename);

    virtual void findDensityOfStatesNonIteratively();
    virtual void findPartitionFunctionsAndDensityOfStates(double tolerance = 1E-7, int maxIterations = 1000);

    virtual void saveLogDensityOfStates(const std::string& filename);       //saves the logarithm. normalized to 0 for the central entry, this also writes out error estimates
    virtual void saveLogDensityOfStatesIsing(const std::string& filename);

    //jackknifes each time series into blockcount blocks (this allows for differing lengths, e.g. if they have been
    //sub-sampled before) --> additional error estimation
    virtual ReweightingResult reweight(double targetBeta);
    virtual HistogramDouble* reweightEnergyHistogram(double targetBeta);
    virtual HistogramDouble* reweightObservableHistogram(double targetBeta,
            unsigned numBins);
    virtual ReweightingResult reweightWithHistograms(double targetBeta,
            unsigned obsBinCount);

    //find the maximum of the susceptibility between betaStart and betaEnd
    //put its location into betaMax, its value into suscMax and add all points evaluated to the map pointsEvaluated
    //also use jackknifing to estimate errors on betaMax and suscMax and store the points evaluated in the passes
    //for each individual jackknife block into pointsEvaluatedJK
    void findMaxObservableSusceptibility(double& betaMax, double& betaMaxError, double& suscMax, double& suscMaxError,
            std::map<double, double>& pointsEvaluated, std::vector<std::map<double,double> >& pointsEvaluatedJK,
        double betaStart, double betaEnd);
    void findMaxSpecificHeatDiscrete(double& betaMax, double& betaMaxError, double& specHeatMax, double& specHeatMaxError,
            std::map<double, double>& pointsEvaluated, std::vector<std::map<double,double> >& pointsEvaluatedJK,
        double betaStart, double betaEnd);
    void findMinBinder(double& betaMin, double& betaMinError,
            double& binderMin, double& binderMinError,
            std::map<double, double>& pointsEvaluated,
            std::vector<std::map<double,double> >& pointsEvaluatedJK,
            double betaStart, double betaEnd);


    //Search between betaStart and betaEnd for a double-peak histogram
    //with as equal peak-heights as possible, obtain histograms from the
    //reweighting. .0<tolerance<1.0 should be chosen small enough to
    //distinguish the dip from the peaks and large enough to distinguish
    //the dip from noise...
    //Return the temperature location for this histogram in betaDouble
    //and the ratio max(max1,max2)/(2*min) in relativeDip. If no dip is found,
    //relativeDip is set to 1.0. A pointer to the final histogram is put
    //into histo.
    //The function for observable histograms has an additional parameter
    //specifying the number of bins to use for the reweighted histograms.
    //Additional jackknife error estimation
    void findEnergyEqualHeight(double& betaDouble, double& betaDoubleError,
            double& relativeDip, double& relDipError,
            HistogramDouble*& histoResult,
            double betaStart, double betaEnd, double tolerance = 0.1);
    void findObservableEqualHeight(double& betaDouble, double& betaDoubleError,
            double& relativeDip, double& relDipError,
            HistogramDouble*& histoResult,
            double betaStart, double betaEnd, unsigned numBins,
            double tolerance = 0.1);

    //The following additionally
    //optimize for equal peak-*weight* (as in [JankeFirstOrder]).
    //One function for both operations as this makes the jackknife
    //error estimation faster.
    void findEnergyEqualHeightWeight(
            double& betaDoubleEH, double& betaDoubleErrorEH,
            double& relativeDipEH, double& relDipErrorEH,
            HistogramDouble*& histoResultEH,
            double& betaDoubleEW, double& betaDoubleErrorEW,
            double& relativeDipEW, double& relDipErrorEW,
            HistogramDouble*& histoResultEW,
            double betaStart, double betaEnd, double tolerance = 0.1);
    void findObservableEqualHeightWeight(
            double& betaDoubleEH, double& betaDoubleErrorEH,
            double& relativeDipEH, double& relDipErrorEH,
            HistogramDouble*& histoResultEH,
            double& betaDoubleEW, double& betaDoubleErrorEW,
            double& relativeDipEW, double& relDipErrorEW,
            HistogramDouble*& histoResultEW,
            double betaStart, double betaEnd, unsigned numBins,
            double tolerance = 0.1);


    //The following just reweight to obtain the histogram at the target beta,
    //return that and the relative dip
    //additional JK error estimation
    void energyRelDip(double& relDip, double& relDipError,
            HistogramDouble*& histoResult,
            double targetBeta, double tolerance = 0.1);
    void obsRelDip(double& relDip, double& relDipError,
            HistogramDouble*& histoResult,
            double targetBeta, unsigned numBins, double tolerance = 0.1);

    //obtain energy and specific heat estimates directly from the density of states
    //(do not re-consider the time series), estimate Jackknife errors
    virtual ReweightingResult reweightDiscrete(double beta);


};

#endif /* MULTIREWEIGHTHISTOSPTJK_H_ */
