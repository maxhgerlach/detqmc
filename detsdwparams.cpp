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
    neededModelPars += "L", "r", "lambda", "accRatio", "bc", "txhor", "txver", "tyhor",
        "tyver", "updateMethod", "spinProposalMethod", "repeatUpdateInSlice",
        "globalShift", "wolffClusterUpdate", "wolffClusterShiftUpdate";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw_ParameterMissing(*p);
        }
    }
    // need mu or (mux and muy)
    if (not (specified.count("mu") or (specified.count("mux") and specified.count("muy")))) {
        throw_ParameterMissing("mu or (mux and muy)");
    }
    
    //Check parameters chosen correctly
    if (model != "sdw") {
        throw_ParameterWrong_message("Parameters specify model: " + model + " instead of sdw");
    }
    if (d != 2) {
        throw_ParameterWrong("d", d);
    }
    if (not (opdim == 1 or opdim == 2 or opdim == 3)) {
        throw_ParameterWrong("opdim", opdim);
    }
    std::string possibleBC[] = {"pbc", "apbc-x", "apbc-y", "apbc-xy"};
    bool bc_is_one_of_the_possible = false;
    for (const std::string& test_bc : possibleBC) {
        if (test_bc == bc_string) bc_is_one_of_the_possible = true;
    }
    if (not bc_is_one_of_the_possible) {
        throw_ParameterWrong("bc", bc_string);
    }

    if (weakZflux and opdim !=2) {
        throw_ParameterWrong_message("Magnetic field specified for opdim=" + numToString(opdim) +
                                     ", but currently only supported for opdim=2");
    }

    std::string possibleUpdateMethods[] = {"iterative", "woodbury", "delayed"};
    bool updateMethod_is_one_of_the_possible = false;
    for (const std::string& test_updateMethod : possibleUpdateMethods) {
        if (test_updateMethod == updateMethod_string) updateMethod_is_one_of_the_possible = true;
    }
    if (not updateMethod_is_one_of_the_possible) {
        throw_ParameterWrong("updateMethod", updateMethod_string);
    }
    if (specified.count("updateMethod") and updateMethod_string == "delayed") {
        if (not specified.count("delaySteps")) {
            throw_ParameterMissing("delaySteps");
        }
        uint32_t N = static_cast<uint32_t>(std::pow(L, 2));
        if (delaySteps <= 0 or delaySteps > N) {
            throw_ParameterWrong("delaySteps", delaySteps);
        }
    }    
    
    std::string possibleSpinProposalMethods[] = {"box", "rotate_then_scale", "rotate_and_scale"};
    bool spinProposalMethod_is_one_of_the_possible = false;
    for (const std::string& test_spinProposalMethod: possibleSpinProposalMethods) {
        if (test_spinProposalMethod == spinProposalMethod_string) spinProposalMethod_is_one_of_the_possible = true;
    }
    if (not spinProposalMethod_is_one_of_the_possible) {
        throw_ParameterWrong("spinProposalMethod", spinProposalMethod_string);
    }

    if ((globalShift or wolffClusterUpdate or wolffClusterShiftUpdate)
        and globalUpdateInterval == 0) {
        throw_ParameterWrong("globalUpdateInterval", globalUpdateInterval);
    }

    if (wolffClusterShiftUpdate and (globalShift or wolffClusterUpdate)) {
        throw_ParameterWrong_message("Either use combined wolffClusterShiftUpdate or individual global updates");
    }

    if (checkerboard and L % 2 != 0) {
        throw_ParameterWrong_message("Checker board decomposition only supported for even linear lattice sizes");
    }

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

    // if fermions are turned off: no delayed updates
    if (turnoffFermions and updateMethod == DELAYED) {
        throw_ParameterWrong_message("Cannot turn off fermions and have delayed updates at the same time");
    }

    // overrelaxation moves: only if fermions are turned off
    if (overRelaxation and not turnoffFermions) {
        throw_ParameterWrong_message("Cannot have overRelaxation moves if the fermions are turned on");
    }

    // computed parameters
    N = L*L;

    // helper to parse repeat* settings
    auto handle_repeatSetting = [&](const std::string& parameterName,
                                    const std::string& parameterStringValue) -> uint32_t {
        uint32_t repeat = 0;
        if (parameterStringValue == "systemSize") {
            repeat = N * m;
        } else if (parameterStringValue == "systemL") {
            repeat = L;
        } else if (parameterStringValue == "systemm") {
            repeat = m;
        } else if (parameterStringValue == "sqrtSystemLm") {
            repeat = static_cast<uint32_t>(std::sqrt( L * m ));
        } else {
            repeat = fromString<uint32_t>(parameterStringValue);
            if (repeat == 0) {
                throw_ParameterWrong(parameterName, parameterStringValue);
            } 
        }
        return repeat;
    };

    // handle repeatWolffPerSweep and repeatOverRelaxation
    updateTemperatureParameters(*this); // it's ok to do this call repeatedly
    if (specified.count("repeatWolffPerSweep")) {
        repeatWolffPerSweep = handle_repeatSetting("repeatWolffPerSweep",
                                                   repeatWolffPerSweep_string);
    }
    if (specified.count("repeatOverRelaxation")) {
        repeatOverRelaxation = handle_repeatSetting("repeatOverRelaxation",
                                                    repeatOverRelaxation_string);
    }

#define IF_NOT_POSITIVE(x) if (specified.count(#x) > 0 and x <= 0)
#define CHECK_POSITIVE(x)   {                   \
        IF_NOT_POSITIVE(x) {                    \
            throw_ParameterWrong(#x, x);   \
        }                                       \
    }
    CHECK_POSITIVE(L);
    CHECK_POSITIVE(repeatWolffPerSweep);
    CHECK_POSITIVE(repeatOverRelaxation);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE
}


MetadataMap ModelParamsDetSDW::prepareMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
#define META_INSERT_TRUE_FALSE(VAR) {meta[#VAR] = (VAR ? "true" : "false");}
    meta["model"] = "sdw";
    meta["opdim"] = numToString(opdim);
    META_INSERT_TRUE_FALSE(checkerboard);
    META_INSERT_TRUE_FALSE(phi2bosons);
    META_INSERT_TRUE_FALSE(phiFixed);
    META_INSERT_TRUE_FALSE(dumpGreensFunction);
    META_INSERT_TRUE_FALSE(turnoffFermions);
    META_INSERT_TRUE_FALSE(turnoffFermionMeasurements);
    META_INSERT_TRUE_FALSE(overRelaxation);
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
    META_INSERT(u);
    META_INSERT(lambda);
    META_INSERT(txhor);
    META_INSERT(txver);
    META_INSERT(tyhor);
    META_INSERT(tyver);
    META_INSERT(cdwU);
    if (specified.count("mux") and specified.count("muy")) {
        META_INSERT(mux);
        META_INSERT(muy);    
    }
    else {
        META_INSERT(mu);
    }
    META_INSERT_TRUE_FALSE(weakZflux);
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
    if (wolffClusterUpdate or wolffClusterShiftUpdate) {
        META_INSERT(repeatWolffPerSweep);
    }
    if (overRelaxation) {
        META_INSERT(repeatOverRelaxation);
    }
    META_INSERT(repeatUpdateInSlice);
#undef META_INSERT
#undef META_INSERT_TRUE_FALSE
    return meta;
}




// template struct ModelParams<DetSDW<CB_NONE>>;
// template struct ModelParams<DetSDW<CB_ASSAAD_BERG>>;
