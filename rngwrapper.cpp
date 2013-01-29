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

std::string RngWrapper::getName() {
	return "DSFMT " + numToString(DSFMT_MEXP);
}

RngWrapper::RngWrapper(unsigned long seed_, int processIndex_)
		: seed(seed_), processIndex(processIndex_), dsfmt() {
	//dSFMT
	//TODO: use full seed

	//this is not a very rigorous approach to have independent random
	//number generators (taken from Katzgraber's introduction to PRNG)
	unsigned long mySeed = abs(
			((seed * 181) * ((processIndex - 83) * 359)) % 104729);

	dsfmt_init_gen_rand(&dsfmt, mySeed);
}

void RngWrapper::saveState() {
	//dSFMT
	char *s2 = dsfmt_state_to_str(&dsfmt, NULL);
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
