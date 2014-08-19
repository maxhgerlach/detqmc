#ifndef DETMODELPARAMS_H
#define DETMODELPARAMS_H

#include <vector>
#include <string>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop            


#include "metadata.h"
#include "exceptions.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
typedef double num;     //possibility to switch to single precision if ever desired
#pragma GCC diagnostic pop





// Template for struct representing model specific parameters
// -- needs to have a proper specialization for each model considered, which actually
//    implements the functions and provides data members
// for a class derived of DetModelGC this should at least be beta, m, s, dtau

// The set specified contains string representations of all parameters
// actually specified.  This allows throwing an exception at the
// appropriate point in the program if a parameter is missing.
template<class Model>
struct ModelParams {
    void check() { }
    MetadataMap prepareMetadataMap() { return MetadataMap(); }

    //If the model supports replica exchange these functions should
    //be implemented for its control parameter
    //void set_exchange_parameter_value(num val) { ... }
    //num  get_exchange_parameter_value() { return ... }    
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version; (void)ar;
    }    
};


//Special handling to allow passing either 'm' or 'beta', but not both.
//'dtau' must always be given.
//Also check that 's' is set matching.
template<class ModelParams>
ModelParams updateTemperatureParameters(ModelParams pars) {
    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "dtau", "s";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (pars.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

//check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x)   {                   \
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

    //we need exactly one of the parameters 'm' and 'beta'
    if (pars.specified.count("beta") != 0 and pars.specified.count("m") != 0) {
    	throw ParameterWrong("Only specify one of the parameters beta and m");
    }
    if (pars.specified.count("m") == 0 and pars.specified.count("beta") == 0) {
    	throw ParameterWrong("Specify either parameter m or beta");
    }


    if (pars.specified.count("m")) {
    	pars.beta = pars.m * pars.dtau;
    } else if (pars.specified.count("beta")) {
    	//this may result in a slightly lower inverse temperature beta
    	//if dtau is not chosen to match well
    	pars.m = uint32_t(pars.beta / pars.dtau);
    	pars.beta = pars.m * pars.dtau;
    }

    // m needs to be larger than s
    while (pars.m <= pars.s) {
    	// throw ParameterWrong("Parameters m=" + numToString(pars.m) + " and s=" + numToString(pars.s)
        //                      + " do not agree. m must be larger than s.");
        pars.s -= 1;
    }
    if (pars.s < 1) {
        throw ParameterWrong("Cannot choose parameter s obeying 0 < s < m =" + numToString(pars.m));
    }

    return pars;
}



// specializations in separate header files: DetHubbardParams, DetSDWParams



#endif /* DETMODELPARAMS_H */

