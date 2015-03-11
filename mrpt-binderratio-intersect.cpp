#include "mrpt-binderratio-intersect.h"


// Callables that we are going to minimize

typedef std::shared_ptr<MultireweightHistosPT> MRPT_Pointer;

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

typedef std::shared_ptr<MultireweightHistosPTJK> MRPTJK_Pointer;

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




enum BC { PBC=0, APBCX=1, APBCY=2, APBCXY=3, NONE };
const std::array<BC, 4> all_BC = {{PBC, APBCX, APBCY, APBCXY}};


class BinderRatioDifferenceBC {
    std::array<MRPT_Pointer>, 4> mrbc1, mrbc2;
public:
    BinderRatioDifference(std::array<MRPT_Pointer>, 4> mrbc1_, std::array<MRPT_Pointer>, 4> mrbc2_) {
        for (BC bc: all_BC) {
            mrbc1[bc] = mrbc1_[bc];
            mrbc2[bc] = mrbc2_[bc];            
        }
    }
    double operator()(double controlParameter) {
        std::array<ReweightedMomentsJK, 4> moments_bc1, moments_bc2;        

        for (BC bc: all_BC) {
            moments_bc1[bc] = mrbc1[bc].reweightObservableMoments(cp);
            moments_bc2[bc] = mrbc2[bc].reweightObservableMoments(cp);            
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
    std::array<MRPTJK_Pointer>, 4> mrbc1, mrbc2;
    unsigned jkBlock;
public:
    BinderRatioDifference(std::array<MRPTJK_Pointer>, 4> mrbc1_, std::array<MRPTJK_Pointer>, 4> mrbc2_,
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
            mrbc1[bc].reweightObservableSecondFourthMomentJK(moments_bc1.o2,
                                                             moments_bc1.o4,
                                                             controlParameter,
                                                             jkBlock);
            mrbc2[bc].reweightObservableSecondFourthMomentJK(moments_bc2.o2,
                                                             moments_bc2.o4,
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
