//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * rngwrapper.cpp
 *
 *  Created on: Oct 12, 2011
 *      Author: gerlach
 */

#include "rngwrapper.h"

extern "C" {
#include "dsfmt/dSFMT.h"
}

#include <fstream>
#include <cstdlib>
#include "tools.h"

std::string RngWrapper::getName() const {
    return "DSFMT " + numToString(DSFMT_MEXP);
}

RngWrapper::RngWrapper(uint32_t seed_, uint32_t processIndex_)
        : seed(seed_), processIndex(processIndex_), dsfmt() {
    //dSFMT
    //TODO: use full seed

    //this is not a very rigorous approach to have independent random
    //number generators (taken from Katzgraber's introduction to PRNG)
    uint32_t mySeed = (uint32_t) std::abs(
            ((seed * 181) * ((processIndex - 83) * 359)) % 104729);

    dsfmt_init_gen_rand(&dsfmt, mySeed);
}

void RngWrapper::saveState() const {
    //dSFMT
    char *s2 = dsfmt_state_to_str(const_cast<dsfmt_t*>(&dsfmt), NULL);
    char filename[100];
    sprintf(filename, "state-dsfmt-rng-process-%d.dat", processIndex);
    copyFileBackUp(std::string(filename));
    FILE *F = fopen(filename, "w");
    fputs(s2, F);
    fclose(F);
    free(s2);
}

void RngWrapper::loadState() {
    //dSFMT
    char filename[100];
    sprintf(filename, "state-dsfmt-rng-process-%d.dat", processIndex);
    FILE *F = fopen(filename, "r");
    assert(F);
    char* err = dsfmt_file_to_state(&dsfmt, F, NULL);
    if (err) {
        throw(err);
    }
    fclose(F);
}

std::string RngWrapper::stateToString() const {
    char* cstring = dsfmt_state_to_str(const_cast<dsfmt_t*>(&dsfmt), NULL);
    std::string s(cstring);
    free(cstring);
    return s;
}

void RngWrapper::stringToState(const std::string& stateString) {
    char* cstring = const_cast<char*>(stateString.c_str());
    dsfmt_str_to_state(&dsfmt, cstring, NULL);
}

