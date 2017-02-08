/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

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
#include "exceptions.h"
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "dlib/cmd_line_parser.h"
#pragma GCC diagnostic pop
#include "tools.h"
#include "datamapwriter.h"
#include "metadata.h"
#include "mrpt-jk.h"
#include "mrpt.h"
#include "mrpt-binderratio-intersect.h"

using namespace std;

//variables in the following anonymous namespace are local to this file
namespace {
    MRPT_Pointer mr1;
    MRPT_Pointer mr2;

    unsigned binCount = 0;

    bool use_jackknife = false;
    unsigned jackknifeBlocks = 0;

    bool non_iterative = true;
    unsigned maxIterations = 10000;
    double iterationTolerance = 1E-7;

    bool be_quiet = false;
    ofstream dev_null("/dev/null");
    string outputDirPrefix;

    string infoFilename1, infoFilename2;

    enum { COL1, COL2 } timeSeriesFormat;
    
    typedef dlib::cmd_line_parser<char>::check_1a_c clp;
    clp parser;

    unsigned subsampleHowMuch = 1;
    unsigned discardSamples = 0;
    bool sortByCp = false;

    bool globalTau = false;
    bool noTau = false;

    string headerSuffix;

    double cpMin, cpMax;
}

void findBinderIntersection() {
    double cp;
    double cpError = 0;

    string observable = mr1->getObservableName();
    assert(observable == mr2->getObservableName());

    string cpName = mr1->getControlParameterName();
    assert(cpName == mr2->getControlParameterName());

    // limit cpMin and cpMax to range available in data
    auto iterator_pair1 = std::minmax_element(mr1->controlParameterValues.begin(),
                                              mr1->controlParameterValues.end());
    double cpMin_mr1 = *(iterator_pair1.first);
    double cpMax_mr1 = *(iterator_pair1.second);
    auto iterator_pair2 = std::minmax_element(mr2->controlParameterValues.begin(),
                                              mr2->controlParameterValues.end());
    double cpMin_mr2 = *(iterator_pair2.first);
    double cpMax_mr2 = *(iterator_pair2.second);

    //// this does not work with the Intel Compiler 13.1
    // cpMin = std::max( {cpMin, cpMin_mr1, cpMin_mr2} );
    // cpMax = std::min( {cpMax, cpMax_mr1, cpMax_mr2} );
    
    // find the maximum of {cpMin, cpMin_mr1, cpMin_mr2}
    double cpMin_to_max[] = {cpMin_mr1, cpMin_mr2};
    for (auto x : cpMin_to_max) {
        cpMin = std::max(cpMin, x);
    }
    // find the minimum of {cpMax, cpMax_mr1, cpMax_mr2}
    double cpMax_to_min[] = {cpMax_mr1, cpMax_mr2};
    for (auto x : cpMax_to_min) {
        cpMax = std::min(cpMax, x);
    }
    
    cout << "Searching for intersection of Binder cumulants between "
         << cpMin << " and " << cpMax;

    bool ok = true;
    if (not use_jackknife) {
        cout << " without error bars." << endl;
        findBinderRatioIntersect(cp, ok, mr1, mr2, cpMin, cpMax);
    } else {
        cout << " with jackknife error bars." << endl;
        findBinderRatioIntersectError(cp, cpError, ok,
                                      std::dynamic_pointer_cast<MultireweightHistosPTJK>(mr1),
                                      std::dynamic_pointer_cast<MultireweightHistosPTJK>(mr2),
                                      cpMin, cpMax, jackknifeBlocks);
    }

    if (not ok) {
        cout << "Apparently findRoot has not converged" << endl;
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
    if (ok) {
        meta[cpName] = numToString(cp, 16);
    }
    if (use_jackknife and ok) {
        meta[cpName + "Error"] = numToString(cpError, 16);
    }
    string comments = "Estimated intersection point of the Binder cumulants of " +
        observable + " from lattice sizes L1 and L2, from MRPT in two instances\n";
    if (use_jackknife) {
        comments += "Jackknife error estimation, blockCount: "
            + numToString(jackknifeBlocks) + "\n";
    }
    writeOnlyMetaData(outputDirPrefix + "mrpt-binder-intersect-l"+L1+"l"+L2+".dat", meta,
                      comments);
}

void setJackknife(bool useJackknife, unsigned blocks) {
    use_jackknife = useJackknife;
    jackknifeBlocks = blocks;
    headerSuffix = (use_jackknife ? "\t error" : "");
}

void setOutputDirectory(const char *dir) {
    string directory = dir;
    if (directory == "") directory = ".";
    outputDirPrefix = directory + "/";
}

void init() {
    mr1.reset();
    mr2.reset();

    mr1 = MRPT_Pointer(use_jackknife ?
           new MultireweightHistosPTJK(jackknifeBlocks, be_quiet ? dev_null : cout) :
           new MultireweightHistosPT(be_quiet ? dev_null : cout));
    mr2 = MRPT_Pointer(use_jackknife ?
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
            switch (timeSeriesFormat) {
            case COL1:
                mr1->addInputTimeSeries_singleColumn(filename, subsampleHowMuch, discardSamples);
                break;
            case COL2:
                mr1->addInputTimeSeries_twoColumn(filename, subsampleHowMuch, discardSamples);
                break;
            }
    	} else if (L == L2) {
            switch (timeSeriesFormat) {
            case COL1:
                mr2->addInputTimeSeries_singleColumn(filename, subsampleHowMuch, discardSamples);
                break;
            case COL2:
                mr2->addInputTimeSeries_twoColumn(filename, subsampleHowMuch, discardSamples);
                break;
            }
    	} else {
            throw GeneralError("Time series " + filename + " with L=" +
                               numToString(L) + " matches neither L1=" + numToString(L1) +
                               "nor L2=" + numToString(L2));
    	}
    }

    if (sortByCp) {
        mr1->sortTimeSeriesByControlParameter();
        mr2->sortTimeSeriesByControlParameter();
    }

    mr1->createHistograms(binCount);
    mr2->createHistograms(binCount);

    mr1->saveH_km(outputDirPrefix + pre1 + "Hkm.table");
    mr2->saveH_km(outputDirPrefix + pre2 + "Hkm.table");

    if (use_jackknife) {
        std::dynamic_pointer_cast<MultireweightHistosPTJK>(mr1)->
            saveH_km_errors(outputDirPrefix + pre1 + "Hkm-errors.table");
        std::dynamic_pointer_cast<MultireweightHistosPTJK>(mr2)->
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
        mr1->measureGlobalInefficiencies();
        mr2->measureGlobalInefficiencies();
    } else {
        mr1->measureBinInefficiencies();
        mr2->measureBinInefficiencies();
    }

    mr1->saveg_km(outputDirPrefix + pre1 + "gkm.table");
    mr2->saveg_km(outputDirPrefix + pre2 + "gkm.table");

    mr1->updateEffectiveCounts();
    mr2->updateEffectiveCounts();

    if (maxIterations > 0) {
        mr1->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
        mr2->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
    }

    mr1->saveLogDensityOfStates(outputDirPrefix + pre1 + "dos.dat");
    mr2->saveLogDensityOfStates(outputDirPrefix + pre2 + "dos.dat");
}


void initFromCommandLine(int argc, char** argv) {
    //Command line parsing
    parser.add_option("help", "display this help message");
    parser.add_option("q", "be less verbose");
    parser.add_option("info1", "info generated by simulation (\"info.dat\") for lattice size L1", 1);
    parser.add_option("info2", "info generated by simulation (\"info.dat\") for lattice size L2", 1);
    parser.add_option("b", "number of energy bins", 1);
    parser.add_option("i", "max number of iterations to determine Z[cp]", 1);
    parser.add_option("t", "tolerance in iterative determination of Z[cp]", 1);

    parser.add_option("non-iterative", "first do a non-iterative estimation of the density of states as in Fenwick, 2008");

    parser.add_option("j", "use jack-knife error estimation, indicate number of blocks", 1);

    parser.add_option("sub-sample", "Sub samples the time series as they are read in (pass number of data points to be put into one sample", 1);
    parser.add_option("d", "Discard the first samples of the time series (pass number of samples to be left out) -- allows for further thermalization.", 1);
    parser.add_option("sort", "Sort replica timeseries by temperatures before doing any processing (simulate canonical data). This could hide correlations.");

    parser.add_option("global-tau", "Do not estimate statistical inefficiencies for individual bins, but only globally for each temperature -- g_km = g_k");
    parser.add_option("no-tau", "Ignore all differences in statistical inefficiencies when reweighting -- g_km = 1");

    parser.add_option("cp-range", "Takes two arguments determining the range of inverse temperatures between which to search for the intersection of the Binder cumulants", 2);

    parser.add_option("time-series-format", "Set to number of columns: 1 (default) or 2; 1-column time series are already sorted by control parameter", 1);

    parser.add_option("outputDirectory", "directory to write intersection search results to", 1);
    
    //further general arguments: file names of energy/observable time series
    parser.parse(argc, argv);

    // //echo whole commandline:
    // cout << "command line: ";
    // for (int arg = 0; arg < argc; ++arg) {
    //     cout << argv[arg] << " ";
    // }
    // cout << endl;

    const char* one_time_opts[] = {"info1", "info2", "b", "j", "i", "sub-sample", "sort"};
    parser.check_one_time_options(one_time_opts);
    const char* incompatible1[] = {"global-tau", "no-tau"};
    parser.check_incompatible_options(incompatible1);

    if (parser.option("help")) {
        cout << "Multihistogram reweighting for time series originating from parallel tempering or canonical simulations" << endl;
        cout << "used to determine the crossing points of the Binder cumulants at two different lattice sizes." << endl;
        cout << "Command line options understood:" << endl;
        parser.print_options(cout);
        cout << endl;
        cout << "Remaining arguments: timeseries files" << endl;
        return;
    }

    if (const clp::option_type& jk = parser.option("j")) {
        setJackknife(true, fromString<unsigned>(jk.argument()));
    }

    be_quiet = parser.option("q");

    infoFilename1 = parser.option("info1").argument();
    std::cout << infoFilename1 << "\n";
    infoFilename2 = parser.option("info2").argument();
    std::cout << infoFilename2 << "\n";    

    if (not parser.option("b")) {
        cerr << "energy bin count not specified (option -b)!" << endl;
        exit(1);
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
        subsampleHowMuch = dlib::sa = ss.argument();
    }

    if (const clp::option_type& dd = parser.option("d")) {
        discardSamples = dlib::sa = dd.argument();
    }

    if (const clp::option_type& tsf_opt = parser.option("time-series-format")) {
        std::string tsf = tsf_opt.argument();
        if (tsf == "1") {
            timeSeriesFormat = COL1;
        } else if (tsf == "2") {
            timeSeriesFormat = COL2;
        } else {
            cerr << "Invalid time series format -- should be \"1\" or \"2\"" << endl;
        }
    }

    if (const clp::option_type& od = parser.option("outputDirectory")) {
        setOutputDirectory(od.argument().c_str());
    }
    
    sortByCp = parser.option("sort");

    globalTau = parser.option("global-tau");
    noTau = parser.option("no-tau");

    const clp::option_type& br = parser.option("cp-range");
    if (not br) {
    	cerr << "Specify control parameter values between which to search for intersection of Binder cumulants (option --cp-range)" << endl;
    	exit(2);
    }
    cpMin = dlib::sa = br.argument(0);
    cpMax = dlib::sa = br.argument(1);
    cout << cpMin << "\t" << cpMax << endl;

    init();

    findBinderIntersection();
}




int main(int argc, char **argv) {
    initFromCommandLine(argc, argv);
    return 0;
}

