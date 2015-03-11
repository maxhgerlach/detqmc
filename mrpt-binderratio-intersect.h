#include <memory>
#include <array>
#include "mrpt.h"
#include "mrpt-jk.h"


typedef std::shared_ptr<MultireweightHistosPT> MRPT_Pointer;
typedef std::shared_ptr<MultireweightHistosPTJK> MRPTJK_Pointer;

void findBinderRatioIntersect(double& cpOut, bool& ok, MRPT_Pointer mr1, MRPT_Pointer mr2,
                              double cpMin, double cpMax);

void findBinderRatioIntersectError(double& cpOut, double& cpErrOut, bool& ok, 
                                   MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                   double cpMin, double cpMax, unsigned jkBlockCount);

enum BC { PBC=0, APBCX=1, APBCY=2, APBCXY=3, NONE };

void findBinderRatioIntersectBC(double& cpOut, bool& ok, 
                                std::array<MRPT_Pointer, 4> mrbc1,
                                std::array<MRPT_Pointer, 4> mrbc2,
                                double cpMin, double cpMax);

void findBinderRatioIntersectBCError(double& cpOut, double& cpErrOut, bool& ok,
                                     std::array<MRPTJK_Pointer, 4> mrbc1,
                                     std::array<MRPTJK_Pointer, 4> mrbc2,
                                     double cpMin, double cpMax, unsigned jkBlockCount);

