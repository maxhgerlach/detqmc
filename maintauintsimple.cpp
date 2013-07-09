/*
 * maintauintsimple.cpp
 *
 *  Created on: Apr 19, 2013
 *      Author: gerlach
 */

#include <string>
#include <iostream>
#include "statistics.h"
#include "dataseriesloader.h"

int main(int argc, char** argv) {
    using namespace std;
    std::cout << argc << endl;
    for (int c = 1; c < argc; ++c) {
        string filename{argv[c]};
        DoubleSeriesLoader loader;
        loader.readFromFile(filename);
        cout << filename << endl;
        cout << " samples: " << loader.getData(0)->size() << endl;
        cout << " tauint: " << tauint(*(loader.getData(0))) << endl;
        cout << " tauint-adaptive: " << tauint_adaptive(*(loader.getData(0))) << endl;
    }
    return 0;
}


