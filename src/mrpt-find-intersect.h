/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

#include <memory>
#include <array>
#include "mrpt.h"
#include "mrpt-jk.h"


// functions to find intersection of certain observables
//
// here we do not add support for averaging over boundary conditions
// (see mrpt-binderratio-intersect.{h,cpp} for that)


typedef std::shared_ptr<MultireweightHistosPT> MRPT_Pointer;
typedef std::shared_ptr<MultireweightHistosPTJK> MRPTJK_Pointer;

// <|m|^4> / <|m|^2>^2
void findBinderRatioIntersect(double& cpOut, bool& ok, MRPT_Pointer mr1, MRPT_Pointer mr2,
                              double cpMin, double cpMax);

void findBinderRatioIntersectError(double& cpOut, double& cpErrOut, bool& ok, 
                                   MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                   double cpMin, double cpMax, unsigned jkBlockCount);

// systemSize * <|m|^2> / L^(2-eta) with eta = 0.25
void findScaledKTSusceptibilityIntersect(double& cpOut, bool& ok, MRPT_Pointer mr1, MRPT_Pointer mr2,
                                         double cpMin, double cpMax);

void findScaledKTSusceptibilityIntersectError(double& cpOut, double& cpErrOut, bool& ok, 
                                              MRPTJK_Pointer mr1, MRPTJK_Pointer mr2,
                                              double cpMin, double cpMax, unsigned jkBlockCount);

