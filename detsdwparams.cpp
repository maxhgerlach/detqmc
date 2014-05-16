#include "detsdwparams.h"

#include <vector>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop

#include "exceptions.h"

void ModelParamsDetSDW::check() {
    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "mu", "L", "r", "lambda", "accRatio", "bc", "txhor", "txver", "tyhor",
        "tyver", "updateMethod", "spinProposalMethod", "repeatUpdateInSlice",
        "globalShift", "wolffClusterUpdate", "wolffClusterShiftUpdate";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }
    //Check parameters chosen correctly
    if (model != "sdw") {
        throw ParameterWrong("Parameters specify model: " + model + " instead of sdw");
    }
    if (d != 2) {
        throw ParameterWrong("d", d);
    }
    std::string possibleBC[] = {"pbc", "apbc-x", "apbc-y", "apbc-xy"};
    bool bc_is_one_of_the_possible = false;
    for (const std::string& test_bc : possibleBC) {
        if (test_bc == bc_string) bc_is_one_of_the_possible = true;
    }
    if (not bc_is_one_of_the_possible) {
        throw ParameterWrong("bc", bc_string);
    }
    std::string possibleUpdateMethods[] = {"iterative", "woodbury", "delayed"};
    bool updateMethod_is_one_of_the_possible = false;
    for (const std::string& test_updateMethod : possibleUpdateMethods) {
        if (test_updateMethod == updateMethod_string) updateMethod_is_one_of_the_possible = true;
    }
    if (not updateMethod_is_one_of_the_possible) {
        throw ParameterWrong("updateMethod", updateMethod_string);
    }
    if (specified.count("updateMethod") and updateMethod_string == "delayed") {
        if (not specified.count("delaySteps")) {
            throw ParameterMissing("delaySteps");
        }
        uint32_t N = static_cast<uint32_t>(std::pow(L, 2));
        if (delaySteps <= 0 or delaySteps > N) {
            throw ParameterWrong("delaySteps", delaySteps);
        }
    }
    std::string possibleSpinProposalMethods[] = {"box", "rotate_then_scale", "rotate_and_scale"};
    bool spinProposalMethod_is_one_of_the_possible = false;
    for (const std::string& test_spinProposalMethod: possibleSpinProposalMethods) {
        if (test_spinProposalMethod == spinProposalMethod_string) spinProposalMethod_is_one_of_the_possible = true;
    }
    if (not spinProposalMethod_is_one_of_the_possible) {
        throw ParameterWrong("spinProposalMethod", spinProposalMethod_string);
    }

    if ((globalShift or wolffClusterUpdate or wolffClusterShiftUpdate)
        and globalUpdateInterval == 0) {
        throw ParameterWrong("globalUpdateInterval", globalUpdateInterval);
    }

    if (wolffClusterShiftUpdate and (globalShift or wolffClusterUpdate)) {
        throw ParameterWrong("Either use combined wolffClusterShiftUpdate or individual global updates");
    }

    if (checkerboard and L % 2 != 0) {
        throw ParameterWrong("Checker board decomposition only supported for even linear lattice sizes");
    }
#define IF_NOT_POSITIVE(x) if (specified.count(#x) > 0 and x <= 0)
#define CHECK_POSITIVE(x)   {                   \
        IF_NOT_POSITIVE(x) {                    \
            throw ParameterWrong(#x, x);   \
        }                                       \
    }
    CHECK_POSITIVE(L);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

    // convert strings to enums
    if (bc_string == "pbc") {
        bc = PBC;
    } else if (bc_string == "apbc-x") {
        bc = APBC_X;
    } else if (bc_string == "apbc-y") {
        bc = APBC_Y;
    } else if (bc_string == "apbc-xy") {
        bc = APBC_XY;
    } else {
        // "safe default"
        bc = PBC;
    }
    if (updateMethod_string == "iterative") {
        updateMethod = ITERATIVE;
    } else if (updateMethod_string == "woodbury") {
        updateMethod = WOODBURY;
    } else if (updateMethod_string == "delayed") {
        updateMethod = DELAYED;
    } else {
        // "safe default"
        updateMethod = ITERATIVE;
    }
    if (spinProposalMethod_string == "box") {
        spinProposalMethod = BOX;
    } else if (spinProposalMethod_string == "rotate_then_scale") {
        spinProposalMethod = ROTATE_THEN_SCALE;
    } else if (spinProposalMethod_string == "rotate_and_scale") {
        spinProposalMethod = ROTATE_AND_SCALE;
    } else {
        // "safe default"
        spinProposalMethod = BOX;
    }

    // computed parameters
    N = L*L;
}


MetadataMap ModelParamsDetSDW::prepareMetadataMap() const {
    MetadataMap meta;
    #define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
    meta["model"] = "sdw";
    meta["checkerboard"] = (checkerboard ? "true" : "false");
    meta["updateMethod"] = updateMethodstr(updateMethod);
    meta["spinProposalMethod"] = spinProposalMethodstr(spinProposalMethod);
    if (spinProposalMethod != BOX) {
        META_INSERT(adaptScaleVariance);
    }
    if (updateMethod == DELAYED) {
        META_INSERT(delaySteps);
    }
    if (bc == PBC) {
        meta["bc"] = "pbc";
    } else if (bc == APBC_X) {
        meta["bc"] = "apbc-x";
    } else if (bc == APBC_Y) {
        meta["bc"] = "apbc-y";
    } else if (bc == APBC_XY) {
        meta["bc"] = "apbc-xy";
    }
    META_INSERT(accRatio);
    META_INSERT(r);
    META_INSERT(lambda);
    META_INSERT(txhor);
    META_INSERT(txver);
    META_INSERT(tyhor);
    META_INSERT(tyver);
    META_INSERT(cdwU);
    META_INSERT(mu);
    META_INSERT(L);
    META_INSERT(d);
    META_INSERT(N);
    META_INSERT(beta);
    META_INSERT(m);
    META_INSERT(dtau);
    META_INSERT(s);
    META_INSERT(globalShift);
    META_INSERT(wolffClusterUpdate);
    META_INSERT(wolffClusterShiftUpdate);
    if (globalShift or wolffClusterUpdate or wolffClusterShiftUpdate) {
        META_INSERT(globalUpdateInterval);
    }
    META_INSERT(repeatUpdateInSlice);
#undef META_INSERT
    return meta;
}




template struct ModelParams<DetSDW<CB_NONE>>;
template struct ModelParams<DetSDW<CB_ASSAAD_BERG>>;