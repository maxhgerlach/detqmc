/*
 * maindet.cpp
 *
 *  Created on: Dec 10, 2012
 *      Author: gerlach
 */


#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include "detqmc.h"


//Parse command line and configuration file to configure the parameters of our simulation.
//In case of invocation with --help or --version, only print some info.
//return a tuple (runSimulation = true or false, simulationParameterStruct)
std::tuple<bool,Params> configureSimulation(int argc, char **argv) {
	bool runSimulation = true;
	Params par;

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
	po::options_description configOptions("Model and simulation parameters, specify via command line or config file");
	configOptions.add_options()
			("model", po::value<string>(&par.model)->default_value("hubbard"))
			("t", po::value<num>(&par.t), "hopping energy")
			("U", po::value<num>(&par.U), "potential energy")
			("mu", po::value<num>(&par.mu), "chemical potential")
			("L", po::value<unsigned>(&par.L), "linear spatial extent")
			("d", po::value<unsigned>(&par.d), "spatial dimension")
			("beta", po::value<num>(&par.beta), "inverse temperature (in units of 1/t, kB=1)")
			("m", po::value<unsigned>(&par.m), "number of imaginary time discretization levels (beta = m*dtau)")
			;
	po::variables_map vm;
	//parse command line
	po::store(po::command_line_parser(argc, argv).options(genericOptions).options(configOptions).run(), vm);
	//parse config file, options specified there have lower precedence
	std::ifstream ifsConf(confFileName);
	po::store(po::parse_config_file(ifsConf, configOptions), vm);

	using std::cout; using std::endl;
	if (vm.count("help")) {
		cout << "Usage:" << endl << endl << genericOptions << endl << configOptions << endl;
		runSimulation = false;
	}

	return std::make_tuple(runSimulation, par);
}


int main(int argc, char **argv) {
	Params par;
	bool runSimulation;
	std::tie(runSimulation, par) = configureSimulation(argc, argv);

	if (runSimulation) {
		DetQMC simulation(par);
	}

	return 0;
}
