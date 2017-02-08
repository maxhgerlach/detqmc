/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * normaldistribution.h
 *
 *  Created on: Feb 12, 2014
 *      Author: gerlach
 */

#ifndef NORMALDISTRIBUTION_H_
#define NORMALDISTRIBUTION_H_

#include <stack>
#include <cmath>
#include "rngwrapper.h"

typedef double num;

class NormalDistribution {
	RngWrapper& rng;
	std::stack<num> standard_normal_random_variables;			//stack is really overkill, but won't matter for now
	void generate_box_muller();
public:
	NormalDistribution(RngWrapper& rng);

	num get(num sigma, num mean);

	void reset();			//clear internal state (standard_normal_random_variables) -- no worrying about serialization
};

inline NormalDistribution::NormalDistribution(RngWrapper& rng_)
	: rng(rng_) {
}

inline void NormalDistribution::reset() {
	standard_normal_random_variables = std::stack<num>();
}


// cf. http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform#Polar_form
// or Numerical Recipes
inline void NormalDistribution::generate_box_muller() {
	using namespace std;

	num v1, v2, rsq;
	do {
		num u1 = rng.rand01();
		num u2 = rng.rand01();
		v1 = 2.0*u1 - 1.0;
		v2 = 2.0*u2 - 1.0;
		rsq = v1*v1 + v2*v2;
	} while (rsq >= 1.0 or rsq == 0.0);		//pick v1, v2 on unit circle

	// [polar] Box-Muller transform
	num fac = sqrt(-2.0 * log(rsq) / rsq);
	standard_normal_random_variables.push(v1 * fac);
	standard_normal_random_variables.push(v2 * fac);
}


inline num NormalDistribution::get(num sigma, num mean) {
	if (standard_normal_random_variables.empty()) {
		generate_box_muller();
	}
	num var = standard_normal_random_variables.top();
	standard_normal_random_variables.pop();

	return mean + sigma * var;
}


#endif /* NORMALDISTRIBUTION_H_ */
