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
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "boost/filesystem.hpp"
#include "dlib/cmd_line_parser.h"
#pragma GCC diagnostic pop
#include "exceptions.h"
#include "tools.h"
#include "datamapwriter.h"
#include "metadata.h"
#include "mrpt-jk.h"
#include "mrpt.h"
#include "mrpt-binderratio-intersect.h"

using namespace std;

//variables in the following anonymous namespace are local to this file
namespace {
//  enum BC { PBC=0, APBCX=1, APBCY=2, APBCXY=3, NONE };
    const std::array<BC, 4> all_BC = {{PBC, APBCX, APBCY, APBCXY}};
    std::array<MRPT_Pointer, 4> mrbc1, mrbc2;
    
    unsigned binCount = 0;

    bool use_jackknife = false;
    unsigned jackknifeBlocks = 0;

    bool non_iterative = true;
    unsigned maxIterations = 10000;
    double iterationTolerance = 1E-7;

    bool be_quiet = false;
    ofstream dev_null("/dev/null");
    string outputDirPrefix;

    std::array<std::string, 4> infoFilenamesBC1, infoFilenamesBC2;

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

    string observable = mrbc1[PBC]->getObservableName();
    assert(observable == mrbc2[PBC]->getObservableName());

    string cpName = mrbc1[PBC]->getControlParameterName();
    assert(cpName == mrbc2[PBC]->getControlParameterName());
    
    cout << "Searching for intersection of Binder cumulants ";

    bool ok = true;
    if (not use_jackknife) {
        cout << " without error bars." << endl;
        findBinderRatioIntersectBC(cp, ok, mrbc1, mrbc2, cpMin, cpMax);
    } else {
        cout << " with jackknife error bars." << endl;
        std::array<MRPTJK_Pointer, 4> mrjkbc1, mrjkbc2;
        for (BC bc: all_BC) {
            mrjkbc1[bc] = std::dynamic_pointer_cast<MultireweightHistosPTJK>(mrbc1[bc]);
            mrjkbc2[bc] = std::dynamic_pointer_cast<MultireweightHistosPTJK>(mrbc2[bc]);            
        }
        
        findBinderRatioIntersectBCError(cp, cpError, ok,
                                        mrjkbc1, mrjkbc2,
                                        cpMin, cpMax, jackknifeBlocks);
    }

    if (not ok) {
        cout << "Apparently findRoot has not converged" << endl;
    }

    MetadataMap meta;
    string L1 = numToString(mrbc1[PBC]->systemL);
    string N1 = numToString(mrbc1[PBC]->systemN);
    string L2 = numToString(mrbc2[PBC]->systemL);
    string N2 = numToString(mrbc2[PBC]->systemN);
    meta["L1"] = L1;
    meta["N1"] = N1;
    meta["L2"] = L2;
    meta["N2"] = N2;
    meta[cpName] = numToString(cp, 16);
    if (use_jackknife) {
        meta[cpName + "Error"] = numToString(cpError, 16);
    }
    string comments = "Estimated intersection point of the Binder cumulants of " +
        observable + " from lattice sizes L1 and L2, from MRPTBC in two instances\n";
    if (use_jackknife) {
        comments += "Jackknife error estimation, blockCount: "
            + numToString(jackknifeBlocks) + "\n";
    }
    writeOnlyMetaData(outputDirPrefix + "mrptbc-binder-intersect-l"+L1+"l"+L2+".dat", meta,
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
    for (BC bc: all_BC) {
        std::shared_ptr<MultireweightHistosPT>& mr_instance1 = mrbc1[bc];
        std::shared_ptr<MultireweightHistosPT>& mr_instance2 = mrbc2[bc];        
        
        if (use_jackknife) {
            mr_instance1 = MRPT_Pointer(new MultireweightHistosPTJK(
                                            jackknifeBlocks,
                                            be_quiet ? dev_null : cout));
            mr_instance2 = MRPT_Pointer(new MultireweightHistosPTJK(
                                            jackknifeBlocks,
                                            be_quiet ? dev_null : cout));
        } else {
            mr_instance1 = MRPT_Pointer(new MultireweightHistosPT(
                                            be_quiet ? dev_null : cout));
            mr_instance2 = MRPT_Pointer(new MultireweightHistosPT(
                                            be_quiet ? dev_null : cout));
        }

        mr_instance1->addSimulationInfo(infoFilenamesBC1[bc]);
        mr_instance2->addSimulationInfo(infoFilenamesBC2[bc]);        
    }

    int L1 = mrbc1[PBC]->systemL;
    int L2 = mrbc2[PBC]->systemL;
    assert(L1 != L2);

    namespace fs = boost::filesystem;
    for (unsigned arg = 0; arg < parser.number_of_arguments(); ++arg) {
    	string series_filename = parser[arg];
    	MetadataMap meta = readOnlyMetadata(series_filename);
    	int L = dlib::sa = meta["L"];

    	if (L == L1) {
            // find bc from infofilename path
            BC this_bc = NONE;
            for (BC bc: all_BC) {
                if (fs::canonical(fs::path(infoFilenamesBC1[bc]).parent_path()) ==
                    fs::canonical(fs::path(series_filename).parent_path().parent_path())) {
                    this_bc = bc;
                    break;
                }
            }
            switch (timeSeriesFormat) {
            case COL1:
                mrbc1[this_bc]->addInputTimeSeries_singleColumn(series_filename, subsampleHowMuch, discardSamples);
                break;
            case COL2:
                mrbc1[this_bc]->addInputTimeSeries_twoColumn(series_filename, subsampleHowMuch, discardSamples);
                break;
            }
    	} else if (L == L2) {
            // find bc from infofilename path
            BC this_bc = NONE;
            for (BC bc: all_BC) {
                if (fs::canonical(fs::path(infoFilenamesBC2[bc]).parent_path()) ==
                    fs::canonical(fs::path(series_filename).parent_path().parent_path())) {
                    this_bc = bc;
                    break;
                }
            }
            switch (timeSeriesFormat) {
            case COL1:
                mrbc2[this_bc]->addInputTimeSeries_singleColumn(series_filename, subsampleHowMuch, discardSamples);
                break;
            case COL2:
                mrbc2[this_bc]->addInputTimeSeries_twoColumn(series_filename, subsampleHowMuch, discardSamples);
                break;
            }
    	} else {
            throw GeneralError("Time series " + series_filename + " with L=" +
                               numToString(L) + " matches neither L1=" + numToString(L1) +
                               "nor L2=" + numToString(L2));
    	}
    }
    if (sortByCp) {
        for (BC bc: all_BC) {
            mrbc1[bc]->sortTimeSeriesByControlParameter();
	    mrbc2[bc]->sortTimeSeriesByControlParameter();
        }
    }
    
    for (BC bc: all_BC) {
        mrbc1[bc]->createHistograms(binCount);
        mrbc2[bc]->createHistograms(binCount);
    }

    if (non_iterative) {
        for (BC bc: all_BC) {
            mrbc1[bc]->findDensityOfStatesNonIteratively();
            mrbc2[bc]->findDensityOfStatesNonIteratively();
        }
    }

    if (noTau) {
        for (BC bc: all_BC) {
            mrbc1[bc]->setBinInefficienciesToUnity();
            mrbc2[bc]->setBinInefficienciesToUnity();
        }
    } else if (globalTau) {
        for (BC bc: all_BC) {
            mrbc1[bc]->measureGlobalInefficiencies();
            mrbc2[bc]->measureGlobalInefficiencies();
        }
    } else {
        for (BC bc: all_BC) {
            mrbc1[bc]->measureBinInefficiencies();
            mrbc2[bc]->measureBinInefficiencies();
        }
    }

    for (BC bc: all_BC) {
        mrbc1[bc]->updateEffectiveCounts();
        mrbc2[bc]->updateEffectiveCounts();
    }

    if (maxIterations > 0) {
        for (BC bc: all_BC) {
            mrbc1[bc]->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
            mrbc2[bc]->findPartitionFunctionsAndDensityOfStates(iterationTolerance, maxIterations);
        }
    }
}


void initFromCommandLine(int argc, char** argv) {
    //Command line parsing
    parser.add_option("help", "display this help message");
    parser.add_option("q", "be less verbose");
    parser.add_option("info1-pbc", "info generated by simulation (\"info.dat\") for lattice size L1, pbc", 1);
    parser.add_option("info2-pbc", "info generated by simulation (\"info.dat\") for lattice size L2, pbc", 1);
    parser.add_option("info1-apbc-x", "info generated by simulation (\"info.dat\") for lattice size L1, apbc-x", 1);
    parser.add_option("info2-apbc-x", "info generated by simulation (\"info.dat\") for lattice size L2, apbc-x", 1);
    parser.add_option("info1-apbc-y", "info generated by simulation (\"info.dat\") for lattice size L1, apbc-y", 1);
    parser.add_option("info2-apbc-y", "info generated by simulation (\"info.dat\") for lattice size L2, apbc-y", 1);
    parser.add_option("info1-apbc-xy", "info generated bxy simulation (\"info.dat\") for lattice size L1, apbc-xy", 1);
    parser.add_option("info2-apbc-xy", "info generated bxy simulation (\"info.dat\") for lattice size L2, apbc-xy", 1);
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

    //echo whole commandline:
    cout << "command line: ";
    for (int arg = 0; arg < argc; ++arg) {
        cout << argv[arg] << " ";
    }
    cout << endl;

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
        return;
    }

    if (const clp::option_type& jk = parser.option("j")) {
        setJackknife(true, fromString<unsigned>(jk.argument()));
    }

    be_quiet = parser.option("q");

    infoFilenamesBC1[PBC]    = parser.option("info1-pbc").argument();
    infoFilenamesBC1[APBCX]  = parser.option("info1-apbc-x").argument();    
    infoFilenamesBC1[APBCY]  = parser.option("info1-apbc-y").argument();    
    infoFilenamesBC1[APBCXY] = parser.option("info1-apbc-xy").argument();    
    infoFilenamesBC2[PBC]    = parser.option("info2-pbc").argument();
    infoFilenamesBC2[APBCX]  = parser.option("info2-apbc-x").argument();    
    infoFilenamesBC2[APBCY]  = parser.option("info2-apbc-y").argument();    
    infoFilenamesBC2[APBCXY] = parser.option("info2-apbc-xy").argument();    

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
    	cerr << "Specify inverse temperatures between which to search for intersection of Binder cumulants (option --cp-range)" << endl;
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

