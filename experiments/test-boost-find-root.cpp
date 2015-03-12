#include "boost/math/tools/roots.hpp"
#include <cmath>
#include <iostream>

using namespace boost::math::tools;

class Test {
public:
    double operator()(const double x) {
        return std::sin(x) - std::cos(x);
    }
};


int main() {
    boost::uintmax_t max_iter=500;
    
    
    auto r1 = toms748_solve( [](double x)->double { return std::sin(x) - std::cos(x); } , 0.0, 2.0,
                             eps_tolerance<double>(7), max_iter );

    std::cout << "root bracketed: [ " << r1.first << " , " << r1.second <<  " ]" << std::endl;
    // std::cout << "f("<< r1.first << ")=" << t(r1.first) << std::endl;
    // std::cout << "f("<< r1.second << ")=" << t(r1.second) << std::endl;
    std::cout << "max_iter=" << max_iter << std::endl;
    
    return 0;         
}
