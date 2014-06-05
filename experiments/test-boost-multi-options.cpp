#include <iostream>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include "boost/program_options.hpp"

namespace po = boost::program_options;

using namespace std;

int main(int argc, char* argv[])
{
    vector<vector<string> > Lists(2);
			
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("List0", po::value<vector<string> >(&Lists[0])->multitoken(), "List0.")
        ("List1", po::value<vector<string> >(&Lists[1])->multitoken(), "List1.")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm); //assign the variables (if they were specified)

    //conf file, lower precedence
    std::ifstream ifsConf("test.conf");
    po::store(po::parse_config_file(ifsConf, desc), vm);
    po::notify(vm);
    

    if (vm.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }

    for(unsigned int list = 0; list < Lists.size(); list++)
    {
        cout << "List: " << list << endl << "-----------" << endl;
        for(unsigned int item = 0; item < Lists[list].size(); item++)
        {
            cout << "List: " << list << " Item: " << item << " " << Lists[list][item] << endl;
        }
    }
	
    return 0;
}
