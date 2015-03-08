#ifndef REWEIGHTEDMOMENTSJK_H_
#define REWEIGHTEDMOMENTSJK_H_

#include <vector>

struct ReweightedMomentsJK {
    // expectation values of the moments for some observable o:
    double o;                   // <o>
    double o2;                  // <o**2>
    double o4;                  // <o**4>
    // jackknife-block estimatats of the above
    std::vector<double> jkBlocks_o;
    std::vector<double> jkBlocks_o2;
    std::vector<double> jkBlocks_o4; 

    // default constructor: no jackknife-block estimates
    ReweightedMomentsJK(double o_ = 0.0,
                        double o2_ = 0.0,
                        double o4_ = 0.0) :
        o(o_), o2(o2_), o4(o4_)
    { }

    ReweightedMomentsJK(double o_, double o2_, double o4_,
                        const std::vector<double>& jkBlocks_o_,
                        const std::vector<double>& jkBlocks_o2_,
                        const std::vector<double>& jkBlocks_o4_) :
        o(o_), o2(o2_), o4(o4_), jkBlocks_o(jkBlocks_o_), jkBlocks_o2(jkBlocks_o2_),
        jkBlocks_o4(jkBlocks_o4_)
    { }    
};


#endif //REWEIGHTEDMOMENTSJK_H_
