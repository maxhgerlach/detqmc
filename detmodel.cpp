/*
 * detmodel.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include "detmodel.h"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#include "exceptions.h"


MatNum computePropagator(num scalar, const MatNum& matrix) {
    using namespace arma;

    VecNum eigval;
    MatNum eigvec;
    eig_sym(eigval, eigvec, matrix);

    return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}

ModelParams updateTemperatureParameters(ModelParams pars) {
    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "beta", "s";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (pars.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

//check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x)   {                                           \
                                IF_NOT_POSITIVE(x) {                    \
                                    throw ParameterWrong(#x, pars.x);   \
                                }                                       \
                            }
    CHECK_POSITIVE(beta);
    CHECK_POSITIVE(m);
    CHECK_POSITIVE(s);
    CHECK_POSITIVE(dtau);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

    if (pars.specified.count("dtau") != 0) {
        if (pars.specified.count("m")) {
            throw ParameterWrong("Only specify one of the parameters m and dtau");
        }
        //require: at least m >= 2*s --> green function calculated from scratch at least for
        //two actual 2 time slices
        if (pars.s * pars.dtau >= pars.beta) {
        	// given value of s too high, find pair of m and s that works
        	// for this beta and dtau (or slightly smaller dtau)
        	uint32_t mI = uint32_t(std::ceil(pars.beta/pars.dtau));  // -> division rounded up --> m not too small
        	pars.s = uint32_t(std::ceil(mI / 2.0));                  // -> we want m >= 2*s !
        	pars.m = 2 * s;
        	pars.dtau = pars.beta / pars.m;
        } else {
        	//this may give rather low values of dtau -- but not for integer inverse temperatures
        	uint32_t n = uint32_t(std::ceil(pars.beta / (pars.s * pars.dtau)));
        	pars.m = pars.s * n;
        	pars.dtau = pars.beta / pars.m;
        }
    } else if (pars.specified.count("m") == 0) {
        throw ParameterMissing("m");
    }
    if (pars.m % pars.s != 0 or pars.m / pars.s < 2) {
        throw ParameterWrong("Parameters m=" + numToString(pars.m) + " and s=" + numToString(pars.s)
                + " do not agree.");
    }

    return pars;
}



