/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */


/*
 * tools.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: gerlach
 */

#include "tools.h"
#include <glob.h>           //POSIX glob()
#include <iostream>

// Taken from:
// http://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system

std::vector<std::string> glob(const std::string& path){
    using namespace std;
    glob_t glob_result;
    glob(path.c_str(), GLOB_TILDE, NULL, &glob_result);
    vector<string> ret;
    for(uint32_t i = 0; i < glob_result.gl_pathc; ++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}



// these can be called from the debugger
void printMatrixReal(const arma::Mat<num>& mat) {
    mat.print(std::cout);
}
void printMatrixComplex(const arma::Mat<cpx>& mat) {
    mat.print(std::cout);
}







