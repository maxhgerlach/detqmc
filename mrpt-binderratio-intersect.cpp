#include <vector>
#include "boost/math/tools/roots.hpp"
#include "mrpt-binderratio-intersect.h"
#include "statistics.h"


// Callables that we are going to minimize

class BinderRatioDifference {
    MRPT_Pointer mr1, mr2;
public:
    BinderRatioDifference(MRPT_Pointer mr1_, MRPT_Pointer mr2_) :
        mr1(mr1_), mr2(mr2_) { }
    double operator()(double controlParameter) {
        double br1 = mr1->reweightObservableBinderRatio(controlParameter);
        double br2 = mr1->reweightObservableBinderRatio(controlParameter);
        double diff = br1 - br2;
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
        double br2 = mr1->reweightObservableBinderRatioJK(controlParameter, jkBlock);
        double diff = br1 - br2;
        return diff;
    }
};





const std::array<BC, 4> all_BC = {{PBC, APBCX, APBCY, APBCXY}};


class BinderRatioDifferenceBC {
    std::array<MRPT_Pointer, 4> mrbc1, mrbc2;
public:
    BinderRatioDifferenceBC(std::array<MRPT_Pointer, 4> mrbc1_, std::array<MRPT_Pointer, 4> mrbc2_) {
        for (BC bc: all_BC) {
            mrbc1[bc] = mrbc1_[bc];
            mrbc2[bc] = mrbc2_[bc];            
        }
    }
    double operator()(double cp) {
        std::array<ReweightedMomentsJK, 4> moments_bc1, moments_bc2;        

        for (BC bc: all_BC) {
            moments_bc1[bc] = mrbc1[bc]->reweightObservableMoments(cp);
            moments_bc2[bc] = mrbc2[bc]->reweightObservableMoments(cp);            
        }

        ReweightedMomentsJK moments_averaged1, moments_averaged2;
        for (BC bc: all_BC) {
            moments_averaged1.o2 += moments_bc1[bc].o2;
            moments_averaged1.o4 += moments_bc1[bc].o4;
            moments_averaged2.o2 += moments_bc2[bc].o2;
            moments_averaged2.o4 += moments_bc2[bc].o4;
        }
        moments_averaged1.o2 /= 4.0;
        moments_averaged1.o4 /= 4.0;
        moments_averaged2.o2 /= 4.0;
        moments_averaged2.o4 /= 4.0;

        double br1 = moments_averaged1.o4 / std::pow(moments_averaged1.o2, 2);
        double br2 = moments_averaged2.o4 / std::pow(moments_averaged2.o2, 2);        

        double diff = br1 - br2;
        return diff;
    }
};



class BinderRatioDifferenceBCJK {
    std::array<MRPTJK_Pointer, 4> mrbc1, mrbc2;
    unsigned jkBlock;
public:
    BinderRatioDifferenceBCJK(std::array<MRPTJK_Pointer, 4> mrbc1_, std::array<MRPTJK_Pointer, 4> mrbc2_,
                              unsigned jkBlock_) {
        jkBlock = jkBlock_;
        for (BC bc: all_BC) {
            mrbc1[bc] = mrbc1_[bc];
            mrbc2[bc] = mrbc2_[bc];            
        }
    }
    double operator()(double controlParameter) {
        std::array<ReweightedMomentsJK, 4> moments_bc1, moments_bc2;        

        for (BC bc: all_BC) {
            mrbc1[bc]->reweightObservableSecondFourthMomentJK(moments_bc1[bc].o2,
                                                              moments_bc1[bc].o4,
                                                              controlParameter,
                                                              jkBlock);
            mrbc2[bc]->reweightObservableSecondFourthMomentJK(moments_bc2[bc].o2,
                                                              moments_bc2[bc].o4,
                                                              controlParameter,
                                                              jkBlock);
        }

        ReweightedMomentsJK moments_averaged1, moments_averaged2;
        for (BC bc: all_BC) {
            moments_averaged1.o2 += moments_bc1[bc].o2;
            moments_averaged1.o4 += moments_bc1[bc].o4;
            moments_averaged2.o2 += moments_bc2[bc].o2;
            moments_averaged2.o4 += moments_bc2[bc].o4;
        }
        moments_averaged1.o2 /= 4.0;
        moments_averaged1.o4 /= 4.0;
        moments_averaged2.o2 /= 4.0;
        moments_averaged2.o4 /= 4.0;

        double br1 = moments_averaged1.o4 / std::pow(moments_averaged1.o2, 2);
        double br2 = moments_averaged2.o4 / std::pow(moments_averaged2.o2, 2);        

        double diff = br1 - br2;
        return diff;
    }
};



// find root between a and b
template<class Callable>
double findRoot(Callable f, double a, double b, bool& ok) {
    using namespace boost::math::tools;
    
    const boost::uintmax_t max_iter_target=500;
    boost::uintmax_t max_iter=max_iter_target;

    auto r1 = toms748_solve( f , a, b,
                             eps_tolerance<double>(7), max_iter );

    ok = (max_iter <= max_iter_target);

    return (r1.second + r1.first) / 2.0;
}




void findBinderRatioIntersect(double& cpOut, bool& ok,
                              MRPT_Pointer mr1, MRPT_Pointer mr2,
                              double cpMin, double cpMax) {
    BinderRatioDifference f(mr1, mr2);
    cpOut = findRoot(f, cpMin, cpMax, ok);
}    


void findBinderRatioIntersectError(double& cpOut, double& cpErrOut, bool& ok,
                                   MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                   double cpMin, double cpMax, unsigned jkBlockCount) {
    //jackknife blocking
    std::vector<double> jkCp_b(jkBlockCount);
    ok = true;
    for (signed b = 0; b < (signed)jkBlockCount; ++b) {
        BinderRatioDifferenceJK f(mr1, mr2, b);
        bool ok_b = true;
        jkCp_b[b] = findRoot(f, cpMin, cpMax, ok_b);
        ok = ok && ok_b;
    }

    jackknife(cpOut, cpErrOut, jkCp_b);
}


void findBinderRatioIntersectBC(double& cpOut, bool& ok, 
                                std::array<MRPT_Pointer, 4> mrbc1,
                                std::array<MRPT_Pointer, 4> mrbc2,
                                double cpMin, double cpMax) {
    BinderRatioDifferenceBC f(mrbc1, mrbc2);
    cpOut = findRoot(f, cpMin, cpMax, ok);    
}

void findBinderRatioIntersectBCError(double& cpOut, double& cpErrOut, bool& ok,
                                     std::array<MRPTJK_Pointer, 4> mrbc1,
                                     std::array<MRPTJK_Pointer, 4> mrbc2,
                                     double cpMin, double cpMax, unsigned jkBlockCount) {
    //jackknife blocking
    std::vector<double> jkCp_b(jkBlockCount);
    ok = true;
    for (signed b = 0; b < (signed)jkBlockCount; ++b) {
        BinderRatioDifferenceBCJK f(mrbc1, mrbc2, b);
        bool ok_b = true;
        jkCp_b[b] = findRoot(f, cpMin, cpMax, ok_b);
        ok = ok && ok_b;
    }

    jackknife(cpOut, cpErrOut, jkCp_b);
    
}


