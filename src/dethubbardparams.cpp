#include "dethubbardparams.h"

#include <vector>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop


void ModelParams<DetHubbard>::check() {
    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "t", "U", "mu", "L", "d", "checkerboard";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw_ParameterMissing(*p);
        }
    }

    if (model != "hubbard") {
        throw_ParameterWrong_message("Parameters specify model: " + model + " instead of hubbard");
    }

    if (checkerboard and L % 2 != 0) {
        throw_ParameterWrong_message("Checker board decomposition only supported for even linear lattice sizes");
    }
    if (checkerboard and d != 2) {
        throw_ParameterWrong_message("Checker board decomposition only supported for 2d lattices");
    }
    if (bc != "pbc") {
        throw_ParameterWrong_message("Boundary conditions " + bc + " not supported for Hubbard model (only pbc)");
    }

    //check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (specified.count(#x) > 0 and x <= 0)
#define CHECK_POSITIVE(x)   {                   \
        IF_NOT_POSITIVE(x) {                    \
            throw_ParameterWrong(#x, x);        \
        }                                       \
    }
    CHECK_POSITIVE(L);
    CHECK_POSITIVE(d);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE
}


MetadataMap ModelParams<DetHubbard>::prepareMetadataMap() const {
    MetadataMap meta;
    meta["model"] = "hubbard";
    meta["checkerboard"] = (checkerboard ? "true" : "false");
//    meta["timedisplaced"] = (timedisplaced ? "true" : "false");
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
    META_INSERT(t);
    META_INSERT(U);
    META_INSERT(mu);
    META_INSERT(L);
    META_INSERT(d);
    META_INSERT(beta);
    META_INSERT(m);
    META_INSERT(dtau);
    META_INSERT(s);
#undef META_INSERT
    return meta;
}

