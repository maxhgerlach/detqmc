/*
 * maindet.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#include "boost/program_options.hpp"
#include "boost/version.hpp"
#include "boost/filesystem.hpp"
#pragma GCC diagnostic warning "-Wunused-parameter"
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include "git-revision.h"
#include "metadata.h"
#include "detqmc.h"
#include "exceptions.h"
#include "timing.h"


//Parse command line and configuration file to configure the parameters of our simulation.
//In case of invocation with --help or --version, only print some info.
//return a tuple (runSimulation = true or false, simulationParameterStruct)
std::tuple<bool,bool,ModelParams,MCParams> configureSimulation(int argc, char **argv) {
	bool runSimulation = true;
	bool resumeSimulation = false;
	ModelParams modelpar;
	MCParams mcpar;

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

	po::options_description modelOptions("Model parameters, specify via command line or config file");
	modelOptions.add_options()
			("model", po::value<string>(&modelpar.model)->default_value("hubbard"), "model to be simulated: hubbard or sdw")
			("timedisplaced", po::bool_switch(&modelpar.timedisplaced)->default_value(false), "also evaluated imaginary time-displaced Green functions and related observables")
			("checkerboard", po::bool_switch(&modelpar.checkerboard)->default_value(false), "use a checkerboard decomposition to compute the propagator of the Hubbard model")
			("r", po::value<num>(&modelpar.r), "parameter tuning SDW transition")
			("t", po::value<num>(&modelpar.t), "hopping energy")
			("U", po::value<num>(&modelpar.U), "potential energy")
			("mu", po::value<num>(&modelpar.mu), "chemical potential")
			("L", po::value<unsigned>(&modelpar.L), "linear spatial extent")
			("d", po::value<unsigned>(&modelpar.d), "spatial dimension")
			("beta", po::value<num>(&modelpar.beta), "inverse temperature (in units of 1/t, kB=1)")
			("temp", po::value<num>(), "temperature (in units of t, kB=1)")
			("dtau", po::value<num>(&modelpar.dtau), "imaginary time discretization step size (beta = m*dtau). Pass either this or m. If dtau is specified, m is chosen to be compatible with s and beta. In turn the value of dtau actually used in the simulation may be smaller than this.")
			("m", po::value<unsigned>(&modelpar.m), "number of imaginary time discretization levels (beta = m*dtau). Pass either this or dtau.")
			("s", po::value<unsigned>(&modelpar.s)->default_value(1), "separation of timeslices where the Green-function is calculated from scratch with stabilized updates.")
			("accRatio", po::value<num>(&modelpar.accRatio)->default_value(0.5), "for SDW: target acceptance ratio for tuning spin update box size")
			;

	po::options_description mcOptions("Parameters for Monte Carlo simulation, specify via command line or config file");
	mcpar.saveInterval = 0;
	mcOptions.add_options()
			("greenUpdate", po::value<std::string>(&mcpar.greenUpdateType)->default_value("stabilized"), "method to use for updating the Green function: simple or stabilized")
			("sweeps", po::value<unsigned>(&mcpar.sweeps), "number of sweeps used for measurements")
			("thermalization", po::value<unsigned>(&mcpar.thermalization), "number of warm-up sweeps")
			("jkBlocks", po::value<unsigned>(&mcpar.jkBlocks)->default_value(1), "number of jackknife blocks for error estimation")
			("timeseries", po::bool_switch(&mcpar.timeseries)->default_value(false), "if specified, write time series of individual measurements to disk")
			("measureInterval", po::value<unsigned>(&mcpar.measureInterval)->default_value(1), "take measurements every [arg] sweeps")
			("saveInterval", po::value<unsigned>(&mcpar.saveInterval), "write measurements to disk every [arg] sweeps; default: only at end of simulation")
			("rngSeed", po::value<unsigned long>(&mcpar.rngSeed), "seed for pseudo random number generator")
			("state", po::value<string>(&mcpar.stateFileName)->default_value("simulation.state"),
					"file, the simulation state will be dumped to.  If it exists, resume the simulation from here.  If you now specify a value for sweeps that is larger than the original setting, an according number of extra-sweeps will be performed.  However, on-the-fly calculation of error bars will no longer work")
			;

	po::variables_map vm;

	//parse command line
	po::options_description cmdlineOptions;
	cmdlineOptions.add(genericOptions).add(modelOptions).add(mcOptions);
	po::store(po::parse_command_line(argc, argv, cmdlineOptions), vm);
	po::notify(vm);

	//parse config file, options specified there have lower precedence
	po::options_description confFileOptions;
	confFileOptions.add(modelOptions).add(mcOptions);
	std::ifstream ifsConf(confFileName);
	po::store(po::parse_config_file(ifsConf, confFileOptions), vm);
	po::notify(vm);

	using std::cout; using std::endl;
	cout << "Assume config file " << confFileName << endl;

	if (boost::filesystem::exists(mcpar.stateFileName)) {
		cout << "Found simulation state file " << stateFileName << ", will resume simulation" << endl;
		resumeSimulation = true;
	}

	if (vm.count("help")) {
		cout << "Usage:" << endl << endl
			 << genericOptions << endl
			 << modelOptions << endl
			 << mcOptions << endl;
		runSimulation = false;
	}
	if (vm.count("version")) {
		runSimulation = false;
	}

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

	return std::make_tuple(runSimulation, resumeSimulation, modelpar, mcpar);
}


int main(int argc, char **argv) {
	std::cout << "Build info:\n"
		<< metadataToString(collectVersionInfo())
		<< "\n";

	ModelParams parmodel;
	MCParams parmc;
	bool runSimulation;
	bool resumeSimulation;
	std::tie(runSimulation, resumeSimulation, parmodel, parmc) = configureSimulation(argc, argv);

	timing.start("total");
	if (runSimulation) {
		if (not resumeSimulation) {
			DetQMC simulation(parmodel, parmc);
			simulation.run();
		} else if (resumeSimulation) {
			DetQMC simulation(parmc.stateFileName, parmc.sweeps);
			simulation.run();
		}
	}
	timing.stop("total");

	return 0;
}
