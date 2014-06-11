#include "boost/mpi.hpp"
#include <vector>
#include <iostream>

namespace mpi = boost::mpi;

int main(int argc, char* argv[])
{
    mpi::environment env(argc, argv);
    mpi::communicator world;

    int values[3];
    for (int i = 0; i < 3; ++i) {
        values[i] = 3*world.rank() + i;
    }

    int* allvalues = 0;
    if (world.rank() == 0) {
        allvalues = new int[3 * world.size()];
    }

    mpi::gather(world,
                values,         // send
                3,
                allvalues,
                0);
    
    if (world.rank() == 0) {
        for (int i = 0; i < 3*world.size(); ++i) {
            std::cout << allvalues[i] << std::endl;
        }
    }

    std::cout << std::endl;

    std::vector<int> allvalues2;
    allvalues2.resize(3 * world.size());

    mpi::gather(world,
                values,         // send
                3,
                allvalues2,     // recv
                0);
    
    if (world.rank() == 0) {
        for (int i = 0; i < 3*world.size(); ++i) {
            std::cout << allvalues2 [i] << std::endl;
        }
    }

    
    return 0;
}
