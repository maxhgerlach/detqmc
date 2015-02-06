/*
 * mrpt-binder-intersect.cpp
 *
 *  Created on: Aug 17, 2012
 *      Author: gerlach
 */

// generalized for SDW DQMC (2015-02-06 - )


#include <map>
#include <vector>
#include <cmath>
#include "mrpt.h"
#include "mrpt-jk.h"
#include "numerics.h"
#include "statistics.h"
#include "mrpt-binder-intersect.h"

using namespace std;

//callable: computes difference of Binder cumulant estimates from two different
//MRPT instances at the same inverse temperature beta
class BinderDiffMinCallable{
	MultireweightHistosPT* mr1;
	MultireweightHistosPT* mr2;
	map<double, double>& pointsEvaluated;
	unsigned jkBlock;
public:
	BinderDiffMinCallable(MultireweightHistosPT* mr1, MultireweightHistosPT* mr2,
			map<double, double>& pointsEvaluated) :
				mr1(mr1), mr2(mr2), pointsEvaluated(pointsEvaluated) { }
	double operator()(double beta) {
		double binder1 = mr1->reweightObservableBinder(beta);
		double binder2 = mr2->reweightObservableBinder(beta);
		double diff = fabs(binder2 - binder1);
		pointsEvaluated[beta] = diff;
		return diff;
	}
};

//callable: computes difference of Binder cumulant estimates from two different
//MRPT instances at the same inverse temperature beta for one jackknife block
class BinderDiffMinCallableJK {
	MultireweightHistosPTJK* mr1;
	MultireweightHistosPTJK* mr2;
	map<double, double>& pointsEvaluated;
	unsigned jkBlock;
public:
	BinderDiffMinCallableJK(MultireweightHistosPTJK* mr1, MultireweightHistosPTJK* mr2,
			map<double, double>& pointsEvaluated, unsigned jkBlock) :
				mr1(mr1), mr2(mr2), pointsEvaluated(pointsEvaluated), jkBlock(jkBlock) { }
	double operator()(double beta) {
		double binder1 = mr1->reweightObservableBinderJK(beta, jkBlock);
		double binder2 = mr2->reweightObservableBinderJK(beta, jkBlock);
		double diff = fabs(binder2 - binder1);
		pointsEvaluated[beta] = diff;
		return diff;
	}
};


void mrptFindBinderIntersectionBeta(double& betaOut,
		double& binderDiffOut,
		map<double,double>& pointsEvaluated,
		MultireweightHistosPT* mr1, MultireweightHistosPT* mr2,
		double betaMin, double betaMax) {
	BinderDiffMinCallable f(mr1, mr2, pointsEvaluated);
    brentMinimize(betaOut, binderDiffOut, f, betaMin, betaMax);
    cout << "Found minimum of Binder intersection at beta = " << betaOut
    	 << " with binderDiff = " << binderDiffOut << endl;
}


void mrptFindBinderIntersectionBetaJK(double& betaOut, double& betaErrorOut,
		double& binderDiffOut, double& binderDiffErrorOut,
		map<double,double>& pointsEvaluated,
		MultireweightHistosPTJK* mr1, MultireweightHistosPTJK* mr2,
		double betaMin, double betaMax, unsigned jkBlockCount) {
    //jackknife blocking
    vector<double> jkBinderDiff_b(jkBlockCount);
    vector<double> jkBeta_b(jkBlockCount);
    #pragma omp parallel for
    for (signed b = 0; b < (signed)jkBlockCount; ++b) {
        BinderDiffMinCallableJK f(mr1, mr2, pointsEvaluated, b);
        brentMinimize(jkBeta_b[b], jkBinderDiff_b[b], f, betaMin, betaMax);
        cout << "  block " << b << ": beta = " << jkBeta_b[b]
             << ", binderDiff = " << jkBinderDiff_b[b] << endl;
    }

    jackknife(binderDiffOut, binderDiffErrorOut, jkBinderDiff_b);
    jackknife(betaOut, betaErrorOut, jkBeta_b);

    cout << "Found minimum of Binder intersection at beta = " << betaOut
         << " +/- " << betaErrorOut
    	 << " with binderDiff = " << binderDiffOut
    	 << " +/- " << binderDiffErrorOut << endl;
}

