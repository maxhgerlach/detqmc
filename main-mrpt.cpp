//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

// generalized for SDW DQMC (2015-02-06 - )

/*
 * main-mrpt.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: gerlach
 */

#include <iostream>
#include <fstream>
#include <exception>
#include "mrpt-highlevel.h"

using namespace std;

int main(int argc, char **argv) {
    try {
        initFromCommandLine(argc, argv);
    } catch (exception& e) {
        cout << e.what() << endl;
        return 1;
    }
    catch (...) {
        cout << "Some error occurred" << endl;
        return 2;
    }
    return 0;
}
