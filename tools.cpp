/*
 * tools.cpp
 *
 *  Created on: Apr 29, 2013
 *      Author: gerlach
 */

#include "tools.h"
#include <glob.h>           //POSIX glob()

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






