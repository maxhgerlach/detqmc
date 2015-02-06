/*
 * mrpt-binder-intersect.h
 *
 *  Created on: Aug 17, 2012
 *      Author: gerlach
 */


// generalized for SDW DQMC (2015-02-06 - ) .... we should use a
// solver for the zero of the difference instead.  Maybe do that in
// Python after all.

#ifndef MRPT_BINDER_INTERSECT_H_
#define MRPT_BINDER_INTERSECT_H_

#include <map>
#include "mrpt.h"
#include "mrpt-jk.h"

//minimize difference of binder cumulant between two instances of MRPT as a function of
//beta, give error estimates too
void mrptFindBinderIntersectionBeta(double& betaOut, double& binderDiffOut,
		std::map<double, double>& pointsEvaluated, MultireweightHistosPT* mr1,
		MultireweightHistosPT* mr2, double betaMin, double betaMax);

//minimize difference of binder cumulant between two instances of MRPT as a function of
//beta, give error estimates too
void mrptFindBinderIntersectionBetaJK(double& betaOut, double& betaErrorOut,
		double& binderDiffOut, double& binderDiffErrorOut,
		std::map<double, double>& pointsEvaluated, MultireweightHistosPTJK* mr1,
		MultireweightHistosPTJK* mr2, double betaMin, double betaMax,
		unsigned jkBlockCount);

#endif /* MRPT_BINDER_INTERSECT_H_ */
