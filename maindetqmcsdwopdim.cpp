#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#include "boost/program_options.hpp"
#include "boost/version.hpp"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/join.hpp"
#pragma GCC diagnostic pop
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include "git-revision.h"
#include "metadata.h"
#include "exceptions.h"
#include "timing.h"
#include "detqmc.h"
#include "detsdwopdim.h"
#include "detsdwparams.h"
#include "detmodelloggingparams.h"


// for versions of the program restricted to only one or two
// variations of O(1), O(2), O(3), define one or two of the macros
//   DETSDW_NO_O1, DETSDW_NO_O2, DETSDW_NO_O3



//Parse command line and configuration file to configure the parameters of our simulation.
//In case of invocation with --help or --version, only print some info.
//return a tuple (runSimulation = true or false, simulationParameterStructs)
std::tuple<bool,bool,DetModelLoggingParams,ModelParamsDetSDW,DetQMCParams> configureSimulation(int argc, char **argv) {
    bool runSimulation = true;
    bool resumeSimulation = false;
    DetModelLoggingParams loggingpar;
    ModelParamsDetSDW modelpar;
    DetQMCParams mcpar;

    //parse options
    namespace po = boost::program_options;
    using std::string;
    string confFileName;

    po::options_description genericOptions("Generic options, command line only");
    genericOptions.add_options()
        ("version,v", "print version information (git hash, build date) and exit")
        ("help", "print help on allowed options and exit")
        ("conf,c", po::value<string>(&confFileName)->default_value("simulation.conf"),
         "specify configuration file to be used; settings in there will be overridden by command line arguments")
        ;
    
#if defined(DETSDW_NO_O3) && defined(DETSDW_NO_O2)
    const uint32_t default_opdim = 1;
#elif defined(DETSDW_NO_O3) && defined(DETSDW_NO_O1)
    const uint32_t default_opdim = 2;
#elif defined(DETSDW_NO_O2) && defined(DETSDW_NO_O1)
    const uint32_t default_opdim = 3;
#else
    const uint32_t default_opdim = 3;
#endif

    po::options_description loggingOptions("DetModel logging parameters, specify via command line or config file");
    loggingOptions.add_options()
        ("logSV", po::value<bool>(&loggingpar.logSV)->default_value(false), "log Green's function singular value range, max, min")
        ("checkAndLogDetRatio", po::value<bool>(&loggingpar.checkAndLogDetRatio)->default_value(false), "verify correctness of local spin update transition probabilities")
        ("checkAndLogGreen", po::value<bool>(&loggingpar.checkAndLogGreen)->default_value(false), "verify correctness of updated Green's functions in local spin update")
        ("logGreenConsistency", po::value<bool>(&loggingpar.logGreenConsistency)->default_value(false), "check green consistency between wrapping / advancing, log differences")
        ;

    po::options_description modelOptions("SDW Model parameters, specify via command line or config file");
    modelOptions.add_options()
        ("model", po::value<string>(&modelpar.model)->default_value("sdw"), "only the sdw model is supported")
        ("turnoffFermions", po::value<bool>(&modelpar.turnoffFermions)->default_value(false), "normally false, if true: simulate a pure O(opdim) model, without considering fermion determinants")
        ("opdim", po::value<uint32_t>(&modelpar.opdim)->default_value(default_opdim), "Dimension of the antiferromagneic order parameter.  O(1), O(2) and O(3) models are supported.  If specified explicitly, must agree with the template instantiations included in the compiled executable")
        ("checkerboard", po::value<bool>(&modelpar.checkerboard)->default_value(false), "use a checkerboard decomposition to compute the propagator for the SDW model")
        ("spinProposalMethod", po::value<std::string>(&modelpar.spinProposalMethod_string)->default_value("box"), "SDW model: method how new field values are proposed for local values: box, rotate_then_scale, or rotate_and_scale")
        ("adaptScaleVariance", po::value<bool>(&modelpar.adaptScaleVariance)->default_value(true), "valid unless spinProposalMethod=='box' -- this controls if the variance of the spin updates should be adapted during thermalization")
        ("updateMethod", po::value<std::string>(&modelpar.updateMethod_string)->default_value("iterative"), "How to do the local updates: iterative, woodbury or delayed")
        ("delaySteps", po::value<uint32_t>(&modelpar.delaySteps)->default_value(16), "parameter to use with delayedUpdates")
        ("phi2bosons", po::value<bool>(&modelpar.phi2bosons)->default_value(false), "if this is true: run calculations with a simple theory (r/2)\\sum_i \\phi_i^2 -- ignores parameter u and spatial terms in the bosonic action")
        ("r", po::value<num>(&modelpar.r), "parameter tuning SDW transition")
        ("u", po::value<num>(&modelpar.u)->default_value(1.0), "non-linear self-coupling of phi")
        ("lambda", po::value<num>(&modelpar.lambda)->default_value(1.0), "fermion-boson coupling")
        ("mu", po::value<num>(&modelpar.mu)->default_value(0.5), "chemical potential (for both orbitals) -- this is superseded, if mux and muy are given")
        ("mux", po::value<num>(&modelpar.mux), "chemical potential, x-orbital separately")
        ("muy", po::value<num>(&modelpar.muy), "chemical potential, y-orbital separately")        
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

#if defined(DETSDW_NO_O1) && defined(DETSDW_NO_O2)
    modelpar.opdim = 3;
#endif //defined(DETSDW_NO_O1) && defined(DETSDW_NO_O2)
#if defined(DETSDW_NO_O1) && defined(DETSDW_NO_O3)
    modelpar.opdim = 2;
#endif //defined(DETSDW_NO_O1) && defined(DETSDW_NO_O2)
#if defined(DETSDW_NO_O2) && defined(DETSDW_NO_O3)
    modelpar.opdim = 1;
#endif //defined(DETSDW_NO_O1) && defined(DETSDW_NO_O2)

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
        ("state", po::value<string>(&mcpar.stateFileName)->default_value("simulation.state"),
         "file, the simulation state will be dumped to.  If it exists, resume the simulation from here.  If you now specify a value for sweeps that is larger than the original setting, an according number of extra-sweeps will be performed.  However, on-the-fly calculation of error bars will no longer work.  Also the headers of timeseries files will still show the wrong number of sweeps")
        ("saveConfigurationStreamText", po::bool_switch(&mcpar.saveConfigurationStreamText)->default_value(false),
         "when measuring, also save raw system configurations to disk, in text format")
        ("saveConfigurationStreamBinary", po::bool_switch(&mcpar.saveConfigurationStreamBinary)->default_value(false),
         "when measuring, also save raw system configurations to disk, in binary format")
        ;

    po::variables_map vm;

    //parse command line
    po::options_description cmdlineOptions;
    cmdlineOptions.add(genericOptions).add(loggingOptions).add(modelOptions).add(mcOptions);
    po::parsed_options parsed = po::command_line_parser(argc, argv)
        .options(cmdlineOptions)
        .allow_unregistered()   // do not throw an exception for unknown program options
        .run();
    po::store(parsed, vm);
    po::notify(vm);

    using std::cout; using std::endl;

    //inform about unknown program options
    std::vector<std::string> unrecognized_options = po::collect_unrecognized(parsed.options, po::exclude_positional);
    if (not unrecognized_options.empty()) {
        cout << "Ignored the following unrecognized options: "
             << boost::algorithm::join(unrecognized_options, ", ")
             << endl << endl;
    }

    if (boost::filesystem::exists(mcpar.stateFileName)) {
        cout << "Found simulation state file " << mcpar.stateFileName << ", will resume simulation" << endl;
        resumeSimulation = true;
    }

    if (vm.count("help")) {
        cout << "Usage:" << endl << endl
             << genericOptions << endl
             << loggingOptions << endl
             << modelOptions << endl
             << mcOptions << endl;
        runSimulation = false;
    }
    if (vm.count("version")) {
        runSimulation = false;
    }

    if (runSimulation) {
        cout << "Assume config file " << confFileName << endl;
    }

    //parse config file, options specified there have lower precedence
    po::options_description confFileOptions;
    confFileOptions.add(loggingOptions).add(modelOptions).add(mcOptions);
    std::ifstream ifsConf(confFileName);
    po::store(po::parse_config_file(ifsConf, confFileOptions), vm);
    po::notify(vm);

    if (vm.count("temp")) {
        if (vm.count("beta")) {
            throw_ConfigurationError("Specify either parameter temp or beta, not both");
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
    record(loggingOptions, loggingpar.specified);
    record(modelOptions, modelpar.specified);
    record(mcOptions, mcpar.specified);



    return std::make_tuple(runSimulation, resumeSimulation, loggingpar, modelpar, mcpar);
}


int main(int argc, char **argv) {
    std::cout << "Build info:\n"
              << metadataToString(collectVersionInfo())
              << "\n";

    DetModelLoggingParams parlogging;
    ModelParamsDetSDW parmodel;
    DetQMCParams parmc;
    bool runSimulation;
    bool resumeSimulation;
    std::tie(runSimulation, resumeSimulation, parlogging, parmodel, parmc) = configureSimulation(argc, argv);

    int return_code = 0;

    timing.start("total");
#define RUN_CASE(cb, opdim) case opdim: {                               \
                                DetQMC<DetSDW<cb, opdim>, ModelParamsDetSDW> simulation(parmodel, parmc, parlogging); \
                                simulation.run();                       \
                                break;                                  \
                            }
#define RESUME_CASE(cb, opdim) case opdim: {                            \
                                   DetQMC<DetSDW<cb, opdim>, ModelParamsDetSDW> simulation(parmc.stateFileName, parmc); \
                                   simulation.run();                    \
                                   break;                               \
                               }
#define DEFAULT_CASE default: \
    std::cerr << "Invalid opdim: " << opdim << "\n"; \
    return_code = 1;                                 \
    break;

    if (runSimulation) {
        uint32_t opdim = parmodel.opdim;
        if (parmodel.checkerboard) { // As long as CheckerboardMethod
                                     // and OPDIM remain template
                                     // parameters we need this
                                     // branching here
            if (not resumeSimulation) {
                switch (opdim) {
                    #ifndef DETSDW_NO_O1
                    RUN_CASE(CB_ASSAAD_BERG, 1)
                    #endif
                    #ifndef DETSDW_NO_O2
		    RUN_CASE(CB_ASSAAD_BERG, 2)
                    #endif
                    #ifndef DETSDW_NO_O3
		    RUN_CASE(CB_ASSAAD_BERG, 3)
                    #endif
                    DEFAULT_CASE
                }
            } else if (resumeSimulation) {
                //only very select parameters given in parmc are updated for the resumed simulation
                switch (opdim) {
                    #ifndef DETSDW_NO_O1
                    RESUME_CASE(CB_ASSAAD_BERG, 1)
                    #endif
                    #ifndef DETSDW_NO_O2
                    RESUME_CASE(CB_ASSAAD_BERG, 2)
                    #endif
                    #ifndef DETSDW_NO_O3
		    RESUME_CASE(CB_ASSAAD_BERG, 3)
                    #endif
                    DEFAULT_CASE
                }
            }
        }
        else {
            if (not resumeSimulation) {
                switch (opdim) {
                    #ifndef DETSDW_NO_O1
                    RUN_CASE(CB_NONE, 1)
                    #endif
                    #ifndef DETSDW_NO_O2
                    RUN_CASE(CB_NONE, 2)
                    #endif
                    #ifndef DETSDW_NO_O3
		    RUN_CASE(CB_NONE, 3)
                    #endif
                    DEFAULT_CASE
                }
            } else if (resumeSimulation) {
                //only very select parameters given in parmc are updated for the resumed simulation
                switch (opdim) {
                    #ifndef DETSDW_NO_O1
                    RESUME_CASE(CB_NONE, 1)
                    #endif
                    #ifndef DETSDW_NO_O2
                    RESUME_CASE(CB_NONE, 2)
                    #endif
                    #ifndef DETSDW_NO_O3
		    RESUME_CASE(CB_NONE, 3)
                    #endif
                    DEFAULT_CASE
                }
            }
        }
    }
#undef RUN_CASE
#undef RESUME_CASE
#undef DEFAULT_CASE
    timing.stop("total");

    return return_code;
}
