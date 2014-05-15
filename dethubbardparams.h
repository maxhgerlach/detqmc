#ifndef DETHUBBARDPARAMS_H
#define DETHUBBARDPARAMS_H

#include <vector>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop

#include "exceptions.h"

#include "detmodelparams.h"


template <bool TimeDisplaced, bool CheckerBoard>
class DetHubbard;


template<>
template<bool TimeDisplaced, bool CheckerBoard>
struct ModelParams<DetHubbard<TimeDisplaced,CheckerBoard>> {
    std::string model;          // should be hubbard

    bool checkerboard;

    num t;
    num U;
    num mu;

    uint32_t L;
    uint32_t d;
    num beta;
    uint32_t m;     //either specify number of timeslices 'm'
    num dtau;       //or timeslice separation 'dtau'
    uint32_t s;     //separation of timeslices where the Green function is calculated
                    //from scratch
    std::string bc;

    std::set<std::string> specified;

    ModelParams() :
        model("hubbard"), checkerboard(),
        t(), U(), mu(), L(), d(), beta(), m(), dtau(), s(), bc(),
        specified()
    { }
    
    void check();
    MetadataMap prepareMetaDataMap() const;
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar  & model & Checkerboard
            t & U & mu & L & d & beta & m
            & dtau & s & bc
            & specified;            
    }    
};


template<>
template<bool TimeDisplaced, bool CheckerBoard>
void ModelParams<DetHubbard<TimeDisplaced,CheckerBoard>>::check() {
    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "t", "U", "mu", "L", "d", "checkerboard";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (pars.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

    if (pars.checkerboard and pars.L % 2 != 0) {
        throw ParameterWrong("Checker board decomposition only supported for even linear lattice sizes");
    }
    if (pars.checkerboard and pars.d != 2) {
        throw ParameterWrong("Checker board decomposition only supported for 2d lattices");
    }
    if (pars.bc != "pbc") {
        throw ParameterWrong("Boundary conditions " + pars.bc + " not supported for Hubbard model (only pbc)");
    }

    //check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x)   {                   \
        IF_NOT_POSITIVE(x) {                    \
            throw ParameterWrong(#x, pars.x);   \
        }                                       \
    }
    CHECK_POSITIVE(L);
    CHECK_POSITIVE(d);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE
}

template<>
template<bool TimeDisplaced, bool CheckerBoard>
MetadataMap ModelParams<DetHubbard<TimeDisplaced,CheckerBoard>>::prepareMetadataMap() const {
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
    


#endif /* DETHUBBARDPARAMS_H */
