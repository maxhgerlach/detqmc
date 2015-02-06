/*
 * main-binder-intersect.cpp
 *
 *  Created on: Aug 17, 2012
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )


// Using multiple histogram reweighting estimate the crossing temperature
// of the Binder cumulant of an observable at two different lattice sizes

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <exception>
#include <map>
#include <string>
#include <boost/filesystem.hpp>
#include "dlib/cmd_line_parser.h"
#include "tools.h"
#include "datamapwriter.h"
#include "metadata.h"
#include "mrpt-jk.h"
#include "mrpt.h"
#include "mrpt-binder-intersect.h"

using namespace std;

//variables in the following anonymous namespace are local to this file
namespace {
    MultireweightHistosPT* mr1 = 0;
    MultireweightHistosPT* mr2 = 0;

    unsigned binCount = 0;
    bool discreteIsingBins = false;

    bool use_jackknife = false;
    unsigned jackknifeBlocks = 0;

    bool non_iterative = true;
    unsigned maxIterations = 10000;
    double iterationTolerance = 1E-7;

    bool be_quiet = false;
    ofstream dev_null("/dev/null");
    string outputDirPrefix;

    string inputDir1, inputDir2;

    string infoFilename1, infoFilename2;

    typedef dlib::cmd_line_parser<char>::check_1a_c clp;
    clp parser;

    unsigned subsample = 1;
    unsigned discardSamples = 0;
    bool sortByBeta = false;
    bool saveTauInt = false;
    string zInFile1, zInFile2;
    string zOutFile1, zOutFile2;

    bool globalTau = false;
    bool noTau = false;
    bool autocorrPlots = false;

    string headerSuffix;

    double betaMin, betaMax;
}

void findBinderIntersection() {
	double beta;
	double betaError = 0;
	double binderDiff;
	double binderDiffError = 0;
	map<double,double> binderDiffPoints;

	string observable = mr1->observable;
	assert(observable == mr2->observable);

	cout << "Searching for intersection of Binder cumulants ";

	if (not use_jackknife) {
		cout << " without error bars." << endl;
		mrptFindBinderIntersectionBeta(beta, binderDiff, binderDiffPoints, mr1, mr2, betaMin, betaMax);
	} else {
		cout << " with jackknife error bars." << endl;
		mrptFindBinderIntersectionBetaJK(beta, betaError, binderDiff, binderDiffError, binderDiffPoints,
				(MultireweightHistosPTJK*)mr1, (MultireweightHistosPTJK*)mr2, betaMin, betaMax, jackknifeBlocks);
	}

    MetadataMap meta;
    string L1 = numToString(mr1->systemL);
    string N1 = numToString(mr1->systemN);
    string L2 = numToString(mr2->systemL);
    string N2 = numToString(mr2->systemN);
    meta["L1"] = L1;
    meta["N1"] = N1;
    meta["L2"] = L2;
    meta["N2"] = N2;
    meta["beta"] = numToString(beta, 16);
    meta["binderDiff"] = numToString(binderDiff, 16);
    if (use_jackknife) {
        meta["betaError"] = numToString(betaError, 16);
        meta["binderDiffError"] = numToString(binderDiffError, 16);
    }
    string comments = "Estimated intersection point of the Binder cumulants of " +
    		observable + " from lattice sizes L1 and L2, from MRPT in two instances\n";
    if (use_jackknife) {
        comments += "Jackknife error estimation, blockCount: "
                + numToString(jackknifeBlocks) + "\n";
    }
    writeOnlyMetaData(outputDirPrefix + "mrpt-binder-intersect-l"+L1+"l"+L2+".dat", meta,
            comments);

    if (not binderDiffPoints.empty()) {
        DoubleMapWriter binderDiffPointsOut;
        binderDiffPointsOut.setData(&binderDiffPoints);
        binderDiffPointsOut.addHeaderText("MRPT estimates of the difference of " + observable + " Binder cumulants");
        binderDiffPointsOut.addMeta("energyBins", binCount);
        binderDiffPointsOut.addMeta("observable", "binderdiff-" + observable);
        binderDiffPointsOut.addMeta("L1", mr1->getSystemL());
        binderDiffPointsOut.addMeta("N1", mr1->getSystemN());
        binderDiffPointsOut.addMeta("L2", mr2->getSystemL());
        binderDiffPointsOut.addMeta("N2", mr2->getSystemN());
        binderDiffPointsOut.addHeaderText("beta\t binderdiff-" + observable);
        binderDiffPointsOut.writeToFile(outputDirPrefix + "mrpt-binderdiff-"+observable+"l"+L1+"l"+L2+".values");
    }
}

void setJackknife(bool useJackknife, unsigned blocks) {
    use_jackknife = useJackknife;
    jackknifeBlocks = blocks;
    headerSuffix = (use_jackknife ? "\t error" : "");
}

void init() {
    destroy(mr1);
    destroy(mr2);

    mr1 = (use_jackknife ?
        new MultireweightHistosPTJK(jackknifeBlocks, be_quiet ? dev_null : cout) :
        new MultireweightHistosPT(be_quiet ? dev_null : cout));
    mr2 = (use_jackknife ?
        new MultireweightHistosPTJK(jackknifeBlocks, be_quiet ? dev_null : cout) :
        new MultireweightHistosPT(be_quiet ? dev_null : cout));

    mr1->addSimulationInfo(infoFilename1);
    mr2->addSimulationInfo(infoFilename2);

    int L1 = mr1->systemL;
    int L2 = mr2->systemL;
    assert(L1 != L2);

    string pre1 = "mrpt-l" + numToString(L1) + "-";
    string pre2 = "mrpt-l" + numToString(L2) + "-";

    for (unsigned arg = 0; arg < parser.number_of_arguments(); ++arg) {
    	string filename = parser[arg];
    	MetadataMap meta = readOnlyMetadata(filename);
    	int L = dlib::sa = meta["L"];

    	if (L == L1) {
    		mr1->addInputTimeSeries(filename, subsample, discardSamples);
    	} else if (L == L2) {
    		mr2->addInputTimeSeries(filename, subsample, discardSamples);
    	} else {
    		throw GeneralException("Time series " + filename + " with L=" +
    				numToString(L) + " matches neither L1=" + numToString(L1) +
    				"nor L2=" + numToString(L2));
    	}
    }

    if (sortByBeta) {
        mr1->sortTimeSeriesByBeta();
        mr2->sortTimeSeriesByBeta();
    }

    if (discreteIsingBins) {
        mr1->createHistogramsIsing();
        mr2->createHistogramsIsing();
    } else {
        mr1->createHistograms(binCount);
        mr2->createHistograms(binCount);
    }
    mr1->saveH_km(outputDirPrefix + pre1 + "Hkm.table");
    mr2->saveH_km(outputDirPrefix + pre2 + "Hkm.table");

    if (use_jackknife) {
        dynamic_cast<MultireweightHistosPTJK*>(mr1)->
                saveH_km_errors(outputDirPrefix + pre1 + "Hkm-errors.table");
        dynamic_cast<MultireweightHistosPTJK*>(mr2)->
                saveH_km_errors(outputDirPrefix + pre2 + "Hkm-errors.table");
    }
    mr1->saveU_m(outputDirPrefix + pre1 + "Um.table");
    mr2->saveU_m(outputDirPrefix + pre2 + "Um.table");

    if (non_iterative) {
        mr1->findDensityOfStatesNonIteratively();
        mr2->findDensityOfStatesNonIteratively();
    }

    if (noTau) {
        mr1->setBinInefficienciesToUnity();
        mr2->setBinInefficienciesToUnity();
    } else if (globalTau) {
        mr1->measureGlobalInefficiencies(autocorrPlots);
        mr2->measureGlobalInefficiencies(autocorrPlots);
    } else {
        mr1->measureBinInefficiencies(autocorrPlots);
        mr2->measureBinInefficiencies(autocorrPlots);
    }

    if (saveTauInt) {
        mr1->writeOutEnergyTauInt(pre1 + "tauint-energy.dat");
        mr1->writeOutObsTauInt(pre1 + "tauint-" + mr1->observable + ".dat");
        mr2->writeOutEnergyTauInt(pre2 + "tauint-energy.dat");
        mr2->writeOutObsTauInt(pre2 + "tauint-" + mr2->observable + ".dat");
    }

    mr1->saveg_km(outputDirPrefix + pre1 + "gkm.table");
    mr2->saveg_km(outputDirPrefix + pre2 + "gkm.table");

    mr1->updateEffectiveCounts();
    mr2->updateEffectiveCounts();

    //load partition functions as starting point for iteration
    if (zInFile1 != "") {
        mr1->loadPartitionFunctions(zInFile1);
    }
    if (zInFile2 != "") {
        mr2->loadPartitionFunctions(zInFile2);
    }

    if (maxIterations > 0) {
        mr1->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
        mr2->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
    }

    if (discreteIsingBins) {
        mr1->saveLogDensityOfStatesIsing(outputDirPrefix + pre1 + "dos.dat");
        mr2->saveLogDensityOfStatesIsing(outputDirPrefix + pre2 + "dos.dat");
    } else {
    	mr1->saveLogDensityOfStates(outputDirPrefix + pre1 + "dos.dat");
        mr2->saveLogDensityOfStates(outputDirPrefix + pre2 + "dos.dat");
    }

    if (zOutFile1 != "") {
        mr1->savePartitionFunctions(zOutFile1);
    }
    if (zOutFile2 != "") {
        mr2->savePartitionFunctions(zOutFile2);
    }
}


void initFromCommandLine(int argc, char** argv) {
    //Command line parsing
    parser.add_option("help", "display this help message");
    parser.add_option("q", "be less verbose");
    parser.add_option("dir1", "directory containing data for lattice size L1", 1);
    parser.add_option("dir2", "directory containing data for lattice size L2", 1);
    parser.add_option("info1", "info generated by simulation (\"info.dat\") for lattice size L1", 1);
    parser.add_option("info2", "info generated by simulation (\"info.dat\") for lattice size L2", 1);
    parser.add_option("b", "number of energy bins", 1);
    parser.add_option("discrete-ising-bins", "Pass this instead of -b to set up energy bins that match the natural discrete energies of the 2D Ising model of this system size.");
    parser.add_option("i", "max number of iterations to determine Z[beta]", 1);
    parser.add_option("t", "tolerance in iterative determination of Z[beta]", 1);
    parser.add_option("loadz1", "load partition functions for lattice size L1 from the indicated file", 1);
    parser.add_option("savez1", "save partition functions for lattice size L1 to the indicated file", 1);

    parser.add_option("loadz2", "load partition functions for lattice size L2 from the indicated file", 1);
    parser.add_option("savez2", "save partition functions for lattice size L2 to the indicated file", 1);

    parser.add_option("non-iterative", "first do a non-iterative estimation of the density of states as in Fenwick, 2008");

    parser.add_option("j", "use jack-knife error estimation, indicate number of blocks", 1);
    parser.add_option("h", "write reweighted histograms at each target beta");

    parser.add_option("sub-sample", "Sub samples the time series as they are read in (pass number of data points to be put into one sample", 1);
    parser.add_option("d", "Discard the first samples of the time series (pass number of samples to be left out) -- allows for further thermalization.", 1);
    parser.add_option("sort", "Sort replica timeseries by temperatures before doing any processing (simulate canonical data). This could hide correlations.");
    parser.add_option("save-tau-int", "Estimate and write out integrated autocorrelation times. Requires --sort");

    parser.add_option("autocorr-plots", "Write out the data points used to estimate (bin) statistical inefficiencies to check the form of the auto-correlation functions");

    parser.add_option("global-tau", "Do not estimate statistical inefficiencies for individual bins, but only globally for each temperature -- g_km = g_k");
    parser.add_option("no-tau", "Ignore all differences in statistical inefficiencies when reweighting -- g_km = 1");

    parser.add_option("beta-range", "Takes two arguments determining the range of inverse temperatures between which to search for the intersection of the Binder cumulants", 2);

    //further general arguments: file names of energy/observable time series
    parser.parse(argc, argv);

    //echo whole commandline:
    cout << "command line: ";
    for (int arg = 0; arg < argc; ++arg) {
        cout << argv[arg] << " ";
    }
    cout << endl;

    const char* one_time_opts[] = {"info1", "info2", "dir1", "dir2", "b", "loadz1", "savez1", "loadz2", "savez2", "j", "i",
    		"sub-sample", "sort"};
    parser.check_one_time_options(one_time_opts);
    const char* incompatible1[] = {"global-tau", "no-tau"};
    parser.check_incompatible_options(incompatible1);
    const char* incompatible2[] = {"autocorr-plots", "no-tau"};
    parser.check_incompatible_options(incompatible2);
    const char* incompatible3[] = {"b", "discrete-ising-bins"};
    parser.check_incompatible_options(incompatible3);

    if (parser.option("help")) {
        cout << "Multihistogram reweighting for time series originating from parallel tempering or canonical simulations" << endl;
        cout << "used to determine the crossing points of the Binder cumulants at two different lattice sizes." << endl;
        cout << "Command line options understood:" << endl;
        parser.print_options(cout);
        cout << endl;
        return;
    }

    if (const clp::option_type& jk = parser.option("j")) {
        setJackknife(true, fromString<unsigned>(jk.argument()));
    }

    if (not parser.option("dir1")) {
    	cerr << "specify directory with data for lattice size L1 (option --dir1)" << endl;
    	exit(1);
    } else {
    	inputDir1 = parser.option("dir1").argument();
    }
    if (not parser.option("dir2")) {
    	cerr << "specify directory with data for lattice size L2 (option --dir2)" << endl;
    	exit(1);
    } else {
    	inputDir2 = parser.option("dir2").argument();
    }

    be_quiet = parser.option("q");

    using namespace boost::filesystem;
    path p1;
    p1 /= inputDir1;
    p1 /= "info.dat";
    path p2;
    p2 /= inputDir2;
    p2 /= "info.dat";

    infoFilename1 = (parser.option("info1") ?
            parser.option("info1").argument() : p1.string());
    infoFilename2 = (parser.option("info2") ?
            parser.option("info2").argument() : p2.string());

    if (not parser.option("b")) {
        if (parser.option("discrete-ising-bins")) {
            discreteIsingBins = true;
        } else {
            cerr << "energy bin count not specified (option -b)!" << endl;
            exit(1);
        }
    } else {
        binCount = dlib::sa = parser.option("b").argument();
    }


    non_iterative = parser.option("non-iterative");
    maxIterations = (non_iterative ? 0 : 10000);            //if non-iterative estimation is attempted, by default don't do any iterations, else default to 1000
    if (parser.option("i")) {
        maxIterations = dlib::sa = parser.option("i").argument();
    }
    if (parser.option("t")) {
        iterationTolerance = dlib::sa = parser.option("t").argument();
    }

    if (const clp::option_type& ss = parser.option("sub-sample")) {
        subsample = dlib::sa = ss.argument();
    }

    if (const clp::option_type& dd = parser.option("d")) {
        discardSamples = dlib::sa = dd.argument();
    }

    sortByBeta = parser.option("sort");
    saveTauInt = sortByBeta and parser.option("save-tau-int");

    zInFile1  = (parser.option("loadz1") ? parser.option("loadz1").argument() : "");
    zOutFile1 = (parser.option("savez1") ? parser.option("savez1").argument() : "");
    zInFile2  = (parser.option("loadz2") ? parser.option("loadz2").argument() : "");
    zOutFile2 = (parser.option("savez2") ? parser.option("savez2").argument() : "");

    globalTau = parser.option("global-tau");
    noTau = parser.option("no-tau");
    autocorrPlots = parser.option("autocorr-plots");

    const clp::option_type& br = parser.option("beta-range");
    if (not br) {
    	cerr << "Specify inverse temperatures between which to search for intersection of Binder cumulants (option --beta-range)" << endl;
    	exit(2);
    }
    betaMin = dlib::sa = br.argument(0);
    betaMax = dlib::sa = br.argument(1);
    cout << betaMin << "\t" << betaMax << endl;

    init();

    findBinderIntersection();
}




int main(int argc, char **argv) {
	initFromCommandLine(argc, argv);
	return 0;
}

