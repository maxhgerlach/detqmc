#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include "boost/program_options.hpp"
#include "boost/version.hpp"
#include "boost/filesystem.hpp"
#include "boost/mpi.hpp"
#pragma GCC diagnostic pop
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include "git-revision.h"
#include "metadata.h"
#include "exceptions.h"
#include "timing.h"
#include "detqmcpt.h"
#include "detsdw.h"
#include "detsdwparams.h"

namespace mpi = boost::mpi;

//Parse command line and configuration file to configure the parameters of our simulation.
//In case of invocation with --help or --version, only print some info.
//return a tuple (runSimulation = true or false, simulationParameterStructs)
std::tuple<bool,bool,ModelParamsDetSDW,DetQMCParams,DetQMCPTParams> configureSimulation(int argc, char **argv) {
    bool runSimulation = true;
    bool resumeSimulation = false;
    ModelParamsDetSDW modelpar;
    DetQMCParams mcpar;
    DetQMCPTParams ptpar;

    //parse options
    namespace po = boost::program_options;
    using std::string;
    string confFileName;

    po::options_description genericOptions("Generic options, command line only");
    genericOptions.add_options()
        ("version", "print version information (git hash, build date) and exit")
        ("help", "print help on allowed options and exit")
        ("conf", po::value<string>(&confFileName)->default_value("simulation.conf"),
         "specify configuration file to be used; settings in there will be overridden by command line arguments")
        ;

    po::options_description modelOptions("SDW Model parameters, specify via command line or config file");
    modelOptions.add_options()
        ("model", po::value<string>(&modelpar.model)->default_value("sdw"), "only the sdw model is supported")
        ("checkerboard", po::value<bool>(&modelpar.checkerboard)->default_value(false), "use a checkerboard decomposition to compute the propagator for the SDW model")
        ("spinProposalMethod", po::value<std::string>(&modelpar.spinProposalMethod_string)->default_value("box"), "SDW model: method how new field values are proposed for local values: box, rotate_then_scale, or rotate_and_scale")
        ("adaptScaleVariance", po::value<bool>(&modelpar.adaptScaleVariance)->default_value(true), "valid unless spinProposalMethod=='box' -- this controls if the variance of the spin updates should be adapted during thermalization")
        ("updateMethod", po::value<std::string>(&modelpar.updateMethod_string)->default_value("iterative"), "How to do the local updates: iterative, woodbury or delayed")
        ("delaySteps", po::value<uint32_t>(&modelpar.delaySteps)->default_value(16), "parameter to use with delayedUpdates")
        ("lambda", po::value<num>(&modelpar.lambda)->default_value(1.0), "fermion-boson coupling")
        ("mu", po::value<num>(&modelpar.mu)->default_value(0.5), "chemical potential")
        ("L", po::value<uint32_t>(&modelpar.L), "linear spatial extent")
        ("d", po::value<uint32_t>(&modelpar.d), "spatial dimension")
        ("beta", po::value<num>(&modelpar.beta), "inverse temperature (in units of 1/t, kB=1)")
        ("temp", po::value<num>(), "temperature (in units of t, kB=1)")
        ("dtau", po::value<num>(&modelpar.dtau), "imaginary time discretization step size (beta = m*dtau). Pass either this or m. If dtau is specified, m is chosen to be compatible with s and beta. In turn the value of dtau actually used in the simulation may be smaller than this.")
        ("m", po::value<uint32_t>(&modelpar.m), "number of imaginary time discretization levels (beta = m*dtau). Pass either this or dtau.")
        ("s", po::value<uint32_t>(&modelpar.s)->default_value(1), "separation of timeslices where the Green-function is calculated from scratch with stabilized updates.")
        ("accRatio", po::value<num>(&modelpar.accRatio)->default_value(0.5), "target acceptance ratio for tuning spin update box size")
        ("bc", po::value<string>(&modelpar.bc_string)->default_value("pbc"), "boundary conditions to use: pbc (periodic), apbc-x, apbc-y or apbc-xy (anti-periodic in x- and/or y-direction)")
        ("txhor", po::value<num>(&modelpar.txhor)->default_value(-1.0), "hopping x left-right")
        ("txver", po::value<num>(&modelpar.txver)->default_value(-0.5), "hopping x up-down")
        ("tyhor", po::value<num>(&modelpar.tyhor)->default_value(0.5), "hopping y left-right")
        ("tyver", po::value<num>(&modelpar.tyver)->default_value(1.0), "hopping y up-down")
        ("cdwU", po::value<num>(&modelpar.cdwU)->default_value(0.0), "coupling to extra interaction that should give an instability to CDW order")
        ("globalUpdateInterval", po::value<uint32_t>(&modelpar.globalUpdateInterval)->default_value(100), "perform global update move every # sweeps [must be even]")
        ("globalShift", po::value<bool>(&modelpar.globalShift)->default_value(false), "perform global constant shift move")
        ("wolffClusterUpdate", po::value<bool>(&modelpar.wolffClusterUpdate)->default_value(false), "perform global Wolff-like single cluster update")
        ("wolffClusterShiftUpdate", po::value<bool>(&modelpar.wolffClusterShiftUpdate)->default_value(false), "perform global Wolff-like single cluster update combined with the global shift move")
        ("repeatUpdateInSlice", po::value<uint32_t>(&modelpar.repeatUpdateInSlice)->default_value(1), "how often to repeat updateInSlice for eacht timeslice per sweep, default: 1")
        ;

    po::options_description mcOptions("Parameters for Monte Carlo simulation, specify via command line or config file");
    mcpar.saveInterval = 0;
    mcOptions.add_options()
        ("greenUpdate", po::value<std::string>(&mcpar.greenUpdateType_string)->default_value("stabilized"), "method to use for updating the Green function: simple or stabilized")
        ("sweeps", po::value<uint32_t>(&mcpar.sweeps), "number of sweeps used for measurements, must be even for serialization consistency")
        ("thermalization", po::value<uint32_t>(&mcpar.thermalization), "number of warm-up sweeps, must be even for serialization consistency")
        ("jkBlocks", po::value<uint32_t>(&mcpar.jkBlocks)->default_value(1), "number of jackknife blocks for error estimation")
        ("timeseries", po::bool_switch(&mcpar.timeseries)->default_value(false), "if specified, write time series of individual measurements to disk")
        ("measureInterval", po::value<uint32_t>(&mcpar.measureInterval)->default_value(1), "take measurements every [arg] sweeps")
        ("saveInterval", po::value<uint32_t>(&mcpar.saveInterval), "write measurements to disk every [arg] sweeps; default: only at end of simulation, must be even for serialization consistency")
        ("rngSeed", po::value<uint32_t>(&mcpar.rngSeed), "seed for pseudo random number generator")
        ("saveConfigurationStreamText", po::bool_switch(&mcpar.saveConfigurationStreamText)->default_value(false),
         "when measuring, also save raw system configurations to disk, in text format")
        ("saveConfigurationStreamBinary", po::bool_switch(&mcpar.saveConfigurationStreamText)->default_value(false),
         "when measuring, also save raw system configurations to disk, in binary format")
        //Multiple processes -- only use standard state file names
        // ("state", po::value<string>(&mcpar.stateFileName)->default_value("simulation.state"),
        //  "file, the simulation state will be dumped to.  If it exists, resume the simulation from here.  If you now specify a value for sweeps that is larger than the original setting, an according number of extra-sweeps will be performed.  However, on-the-fly calculation of error bars will no longer work.  Also the headers of timeseries files will still show the wrong number of sweeps")
        ;

    po::options_description ptOptions("Parameters for SDW model replica exchange / parallel tempering, specify via command line or config file");
    ptOptions.add_options()
        ("rValues", po::value<std::vector<num>>(&ptpar.controlParameterValues)->multitoken(), "values for r, the parameter tuning SDW transition")
        ("exchangeInterval", po::value<uint32_t>(&ptpar.exchangeInterval)->default_value(0), "interval in sweeps between attempted replica exchanges; 0 means \"never\"")
        ;
    ptpar.controlParameterName = "r";

    po::variables_map vm;

    //parse command line
    // short options like -v are disabled to not confuse the multitoken parser
    // for option rValues in case of negative options
    po::options_description cmdlineOptions;
    cmdlineOptions.add(genericOptions).add(modelOptions).add(mcOptions).add(ptOptions);
    po::store(po::parse_command_line(argc, argv, cmdlineOptions,
                                     po::command_line_style::unix_style
                                     ^ po::command_line_style::allow_short), vm);
    po::notify(vm);

    using std::cout; using std::endl;

    //state file -- depends on MPI process rank
    mpi::communicator world;
    int processRank = world.rank();
    mcpar.stateFileName = "simulation." + numToString(processRank) + ".state";    
    
    if (boost::filesystem::exists(mcpar.stateFileName)) {
        cout << "p" << processRank << ": Found simulation state file " << mcpar.stateFileName << ", will resume simulation" << endl;
        resumeSimulation = true;
    }

    if (vm.count("help")) {
        if (processRank == 0) {
            cout << "Usage:" << endl << endl
                 << genericOptions << endl
                 << modelOptions << endl
                 << mcOptions << endl
                 << ptOptions << endl;
        }
        runSimulation = false;
    }
    if (vm.count("version")) {
        runSimulation = false;
    }

    if (runSimulation) {
        if (processRank == 0) {
            cout << "Assume config file " << confFileName << endl;
        }
    }

    //parse config file, options specified there have lower precedence
    po::options_description confFileOptions;
    confFileOptions.add(modelOptions).add(mcOptions).add(ptOptions);
    std::ifstream ifsConf(confFileName);
    po::store(po::parse_config_file(ifsConf, confFileOptions), vm);
    po::notify(vm);

    if (vm.count("temp")) {
        if (vm.count("beta")) {
            throw ConfigurationError("Specify either parameter temp or beta, not both");
        } else {
            //manually specify option beta, remove option temp
            num beta = 1.0 / vm["temp"].as<num>();
            vm.insert(std::make_pair("beta", po::variable_value(beta, false)));
            modelpar.beta = beta;
            vm.erase("temp");
            po::notify(vm);
        }
    }

    //record which options have been specified
    auto record = [vm](const po::options_description& optDesc,
                       std::set<std::string>& set) {
        auto opts = optDesc.options();
        for (auto p = opts.begin(); p != opts.end(); ++p) {
            const std::string& o = (*p)->long_name();
            if (vm.count(o)) {
                set.insert(o);
            }
        }
    };
    record(modelOptions, modelpar.specified);
    record(mcOptions, mcpar.specified);
    record(ptOptions, ptpar.specified);    
    if (vm.count("rValues")) {
        ptpar.specified.insert("controlParameterValues");
    }
    ptpar.specified.insert("controlParameterName");
    
    return std::make_tuple(runSimulation, resumeSimulation, modelpar, mcpar, ptpar);
}


int main(int argc, char **argv) {
    mpi::environment env(argc, argv);
    mpi::communicator world;    // default constructor corresponds to MPI_COMM_WORLD

    int processRank = world.rank();

    if (processRank == 0) {
        std::cout << "Build info:\n"
                  << metadataToString(collectVersionInfo())
                  << "\n";
    }
    
    ModelParamsDetSDW parmodel;
    DetQMCParams parmc;
    DetQMCPTParams parpt;    
    bool runSimulation;
    bool resumeSimulation;
    std::tie(runSimulation, resumeSimulation, parmodel, parmc, parpt) = configureSimulation(argc, argv);

    // if (parmodel.timedisplaced) {
    // 	std::cerr << "timedisplaced routines not implemented currently \n";
    // 	return 1;
    // }

    timing.start("total");
    if (runSimulation) {
        if (parmodel.checkerboard) { // As long as CheckerboardMethod remains a template parameter we need this branching here
            if (not resumeSimulation) {
                DetQMCPT<DetSDW<CB_ASSAAD_BERG>, ModelParamsDetSDW> simulation(parmodel, parmc, parpt);
                simulation.run();
            } else if (resumeSimulation) {
                // parmc.stateFileName -> is set process dependent
                DetQMCPT<DetSDW<CB_ASSAAD_BERG>, ModelParamsDetSDW> simulation(parmc.stateFileName, parmc);
                //only very select parameters given in parmc are updated for the resumed simulation
                simulation.run();
            }
        }
        else {
            if (not resumeSimulation) {
                // parmc.stateFileName -> is set process dependent
                DetQMCPT<DetSDW<CB_NONE>, ModelParamsDetSDW> simulation(parmodel, parmc, parpt);
                simulation.run();
            } else if (resumeSimulation) {
                DetQMCPT<DetSDW<CB_NONE>, ModelParamsDetSDW> simulation(parmc.stateFileName, parmc);
                //only very select parameters given in parmc are updated for the resumed simulation
                simulation.run();
            }
        }
    }
    timing.stop("total");

    return 0;
}
