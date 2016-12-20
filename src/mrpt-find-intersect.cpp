#include <vector>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#include "boost/math/tools/roots.hpp"
#pragma GCC diagnostic pop
#include "mrpt-binderratio-intersect.h"
#include "statistics.h"


// Callables that we are going to minimize

// <|m|^4> / <|m|^2>^2

class BinderRatioDifference {
    MRPT_Pointer mr1, mr2;
public:
    BinderRatioDifference(MRPT_Pointer mr1_, MRPT_Pointer mr2_) :
        mr1(mr1_), mr2(mr2_) { }
    double operator()(double controlParameter) {
        double br1 = mr1->reweightObservableBinderRatio(controlParameter);
        double br2 = mr2->reweightObservableBinderRatio(controlParameter);
        double diff = br1 - br2;
        std::cout << diff << std::endl;
        return diff;
    }
};

class BinderRatioDifferenceJK {
    MRPTJK_Pointer mr1, mr2;
    unsigned jkBlock;
public:
    BinderRatioDifferenceJK(MRPTJK_Pointer mr1_, MRPTJK_Pointer mr2_, unsigned jkBlock_) :
        mr1(mr1_), mr2(mr2_), jkBlock(jkBlock_) { }
    double operator()(double controlParameter) {
        double br1 = mr1->reweightObservableBinderRatioJK(controlParameter, jkBlock);
        double br2 = mr2->reweightObservableBinderRatioJK(controlParameter, jkBlock);
        double diff = br1 - br2;
        return diff;
    }
};

// systemSize * <|m|^2> / L^(2-eta) with eta = 0.25

class ScaledKTSusceptibilityDifference {
    MRPT_Pointer mr1, mr2;
    unsigned systemL1, systemL2;
    double divisor1, divisor2;
public:
    ScaledKTSusceptibilityDifference(MRPT_Pointer mr1_, MRPT_Pointer mr2_) :
        mr1(mr1_), mr2(mr2_),
        systemL1(mr1->getSystemL()), systemL2(mr2->getSystemL()),
        divisor1(std::pow(double(systemL1), 2.0 - 0.25)),
        divisor2(std::pow(double(systemL2), 2.0 - 0.25)) { }
    double operator()(double controlParameter) {
        double s1 = mr1->reweightObservableSusceptibilityPart(controlParameter) / divisor1;
        double s2 = mr2->reweightObservableSusceptibilityPart(controlParameter) / divisor2;
        double diff = s1 - s2;
        std::cout << diff << std::endl;
        return diff;
    }
};

class ScaledKTSusceptibilityDifferenceJK {
    MRPTJK_Pointer mr1, mr2;
    unsigned jkBlock;
    unsigned systemL1, systemL2;
    double divisor1, divisor2;
public:
    ScaledKTSusceptibilityDifferenceJK(MRPTJK_Pointer mr1_, MRPTJK_Pointer mr2_, unsigned jkBlock_) :
        mr1(mr1_), mr2(mr2_), jkBlock(jkBlock_),
        systemL1(mr1->getSystemL()), systemL2(mr2->getSystemL()),
        divisor1(std::pow(double(systemL1), 2.0 - 0.25)),
        divisor2(std::pow(double(systemL2), 2.0 - 0.25)) { }
    double operator()(double controlParameter) {
        double s1 = mr1->reweightObservableSusceptibilityPartJK(controlParameter, jkBlock) / divisor1;
        double s2 = mr2->reweightObservableSusceptibilityPartJK(controlParameter, jkBlock) / divisor2;
        double diff = s1 - s2;
        return diff;
    }
};




// find root between a and b
template<class Callable>
double findRoot(Callable f, double a, double b, bool& ok) {
    using namespace boost::math::tools;
    
    const boost::uintmax_t max_iter_target=500;
    boost::uintmax_t max_iter=max_iter_target;

    try {
        auto r1 = toms748_solve( f , a, b,
                                 eps_tolerance<double>(7), max_iter );

        ok = (max_iter <= max_iter_target);

        return (r1.second + r1.first) / 2.0;
    } catch (const std::exception& e) {
        std::cerr << "findRoot failed. what(): " << e.what() << std::endl;
        ok = false;
        return 0;
    }
}


template<class DifferenceCallable>
void findIntersectImplementation(double& cpOut, bool& ok,
                                 MRPT_Pointer mr1, MRPT_Pointer mr2,
                                 double cpMin, double cpMax) {
    DifferenceCallable f(mr1, mr2);
    cpOut = findRoot(f, cpMin, cpMax, ok);
}

template<class DifferenceJKCallable>
void findIntersectErrorImplementation(double& cpOut, double& cpErrOut, bool& ok,
                                      MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                      double cpMin, double cpMax, unsigned jkBlockCount) {
    //jackknife blocking
    std::vector<double> jkCp_b(jkBlockCount);
    ok = true;
    for (signed b = 0; b < (signed)jkBlockCount; ++b) {
        DifferenceJKCallable f(mr1, mr2, b);
        bool ok_b = true;
        jkCp_b[b] = findRoot(f, cpMin, cpMax, ok_b);
        ok = ok && ok_b;
    }

    jackknife(cpOut, cpErrOut, jkCp_b);
}



void findBinderRatioIntersect(double& cpOut, bool& ok,
                              MRPT_Pointer mr1, MRPT_Pointer mr2,
                              double cpMin, double cpMax) {
    findIntersectImplementation<BinderRatioDifference>(
        cpOut, ok, mr1, mr2, cpMin, cpMax);
}


void findBinderRatioIntersectError(double& cpOut, double& cpErrOut, bool& ok,
                                   MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                   double cpMin, double cpMax, unsigned jkBlockCount) {
    findIntersectErrorImplementation<BinderRatioDifferenceJK>(
        cpOut, cpErrOut, ok, mr1, mr2, cpMin, cpMax, jkBlockCount);
}


void findScaledKTSusceptibilityIntersect(double& cpOut, bool& ok,
                              MRPT_Pointer mr1, MRPT_Pointer mr2,
                              double cpMin, double cpMax) {
    findIntersectImplementation<ScaledKTSusceptibilityDifference>(
        cpOut, ok, mr1, mr2, cpMin, cpMax);
}


void findScaledKTSusceptibilityIntersectError(double& cpOut, double& cpErrOut, bool& ok,
                                   MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                   double cpMin, double cpMax, unsigned jkBlockCount) {
    findIntersectErrorImplementation<ScaledKTSusceptibilityDifferenceJK>(
        cpOut, cpErrOut, ok, mr1, mr2, cpMin, cpMax, jkBlockCount);
}





