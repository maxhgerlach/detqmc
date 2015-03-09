
//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * mrpt-highlevel.h
 *
 *  Created on: Aug 5, 2011
 *      Author: max
 */

#ifndef MRPT_HIGHLEVEL_H_
#define MRPT_HIGHLEVEL_H_


// Common routines
//////////////////

//set options first, then call init
void setOutputDirectory(const char* dir);       //if set: write all files to this sub-directory
void setSubsample(unsigned samplesSize = 1);
void setBins(unsigned bins);
void setJackknife(bool useJackknife, unsigned blocks = 1);
void setQuiet(bool quiet);

//for finding density of states
void setMaxIterations(unsigned maxIterations);
void setTolerance(double tolerance);

void getOriginalControlParameterValues(int* outK, double** outArray1);



// Single boundary conditions, single instance of MRPT
//////////////////////////////////////////////////////

//functions wrapped around an instance of the multireweighting class -- this is
//used in the main source file, can be exported as a python module
//more functionality could easily be added later on (e.g. connection to numpy...)


//set options first, then call init
void setInfoFilename(const char* filename);


//pass command line options to be parsed -- this is the only way to add input
//time series
//if the relevant options are found, call the other functions too
void initFromCommandLine(int argc, char** argv);

void init();    //after this we can proceed estimating quantities using reweighting

void directResults();   //compute averages at original control parameters (cps), no reweighting

void reweight(double cp, bool createHistogramsToo = false);
void reweightRange(double cpMin, double cpMax, double cpStep, bool createHistogramsToo = false);
void reweightDiscreteRange(double cpMin, double cpMax, double cpStep, bool createHistogramsToo = false);

void reweightEnergyHistogram(double cp);
void reweightObservableHistogram(double cp);

void findMaxSusc(double cpStart, double cpEnd);
void findMaxSpecificHeat(double cpStart, double cpEnd);
void findMinBinder(double cpStart, double cpEnd);
void findEnergyDoublePeak(double cpStart, double cpEnd, double tolerance = 0.1);
void findObservableDoublePeak(double cpStart, double cpEnd, double tolerance = 0.1);

void findEnergyRelDip(double targetControlParameter, double tolerance = 0.1);
void findObservableRelDip(double targetControlParameter, double tolerance = 0.1);


//accessor functions to arrays (these give pointers directly to internal data!)
void getEnergyTimeSeries(unsigned k, int* outN_k, double** outArray1);
void getObservableTimeSeries(unsigned k, int* outN_k, double** outArray1);
void getBetaIndexTimeSeries(unsigned k, int* outN_k, int** outArray1);
void getU_m(int* outM, double** outArray1);
void getH_km(int* outK, int* outM, int** outArray2);
void getH_m(int* outM, int** outArray1);
void getg_km(int* outK, int* outM, double** outArray2);
void getH_lm(int* outL, int* outM, int** outArray2);
void getN_kl(int* outK, int* outL, int** outArray2);
void getHeff_m(int* outM, double** outArray1);
void getNeff_lm(int* outL, int* outM, double** outArray2);




// Multiple boundary conditions (to be averaged), 4 instances of MRPT
/////////////////////////////////////////////////////////////////////

//set options first, then call init
void setInfoFilenamesBC(const char* filename1, const char* filename2, const char* filename3, const char* filename4);

//pass command line options to be parsed -- this is the only way to add input
//time series
//if the relevant options are found, call the other functions too
void initFromCommandLineBC(int argc, char** argv);

void initBC();    //after this we can proceed estimating quantities using reweighting

void reweightBC(double cp);
void reweightRangeBC(double cpMin, double cpMax, double cpStep);
void directResultsBC();







#endif /* MRPT_HIGHLEVEL_H_ */
