/*
 * dethubbard.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: gerlach
 */

#include <vector>
#include <cmath>
#include <cassert>
#include <boost/assign/std/vector.hpp>    // 'operator+=()' for vectors
#include "tools.h"
#include "exceptions.h"
#include "rngwrapper.h"
#include "dethubbard.h"
#include "parameters.h"

//using std::acosh;    //Intel compiler chokes with std::acosh


std::unique_ptr<DetHubbard> createDetHubbard(const ModelParams& pars) {
	//check parameters
	using namespace boost::assign;
	std::vector<std::string> neededModelPars;
	neededModelPars += "t", "U", "mu", "L", "d", "beta", "m";
	for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
		if (pars.specified.count(*p) == 0) {
			throw ParameterMissing(*p);
		}
	}

	return std::unique_ptr<DetHubbard>(new DetHubbard(
			pars.t, pars.U, pars.mu, pars.L, pars.d, pars.beta, pars.m)	);
}



DetHubbard::DetHubbard(num t, num U, num mu, unsigned L, unsigned d, num beta, unsigned m) :
		t(t), U(U), mu(mu), L(L), d(d),
		latticeCoordination(2*d), N(static_cast<unsigned>(uint_pow(L,d))),
		beta(beta), m(m), dtau(beta/m),
		alpha(acosh(std::exp(dtau * U * 0.5))),
		nearestNeigbors(2*d, N),			//coordination number: 2*d
		tmat(N, N), proptmat(N,N),
		auxfield(N, m),              //m columns of N rows
		gUp(N,N,m), gDn(N,N,m),      //m slices of N columns x N rows
		obsNames(), obsShorts(), obsValPointers(), obsCount(0)
{
	createNeighborTable();
	setupRandomAuxfield();
	setupTmat();
	using namespace boost::assign;         // bring operator+=() into scope
	obsNames += "occupation spin up", "occupation spin down", "total occupation",
			"kinetic energy", "potential energy", "total energy";
	obsShorts += "nUp", "nDown", "n", "e_t", "e_U", "e";
	obsValPointers += &occUp, &occDn, &occTotal, &eKinetic, &ePotential, &eTotal;
	assert(obsNames.size() == obsShorts.size());
	assert(obsNames.size() == obsValPointers.size());
	obsCount = obsNames.size();
}

DetHubbard::~DetHubbard() {
}


MetadataMap DetHubbard::prepareModelMetadataMap() {
	MetadataMap meta;
	meta["model"] = "hubbard";
	//TODO: simplify the following with a macro...
	meta["t"] = numToString(t);
	meta["U"] = numToString(U);
	meta["mu"] = numToString(mu);
	meta["L"] = numToString(L);
	meta["d"] = numToString(d);
	meta["N"] = numToString(N);
	meta["beta"] = numToString(beta);
	meta["m"] = numToString(m);
	meta["dtau"] = numToString(dtau);
	return meta;
}

void DetHubbard::sweepSimple() {
	for (unsigned timeslice = 0; timeslice < m; ++timeslice) {
		gUp.slice(timeslice) = computeGreenFunction(timeslice, Spin::Up);
		gDn.slice(timeslice) = computeGreenFunction(timeslice, Spin::Down);
		for (unsigned site = 0; site < N; ++site) {
			num ratio =  weightRatioSingleFlip(site, timeslice);
			//Metropolis
			if (ratio > 1 or rng.rand01() < ratio) {
				updateGreenFunctionsAfterFlip(site, timeslice);
			}
		}
	}
}

void DetHubbard::measure() {
	//used to measure occupation:
	num sum_GiiUp = 0;
	num sum_GiiDn = 0;
	//used to measure kinetic energy:
	num sum_GneighUp = 0;
	num sum_GneighDn = 0;
	//used to measure potential energy:
	num sum_GiiUpDn = 0;
	for (unsigned timeslice = 0; timeslice < m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			//use diagonal elements of Green functions:
			sum_GiiUp += gUp(site, site, timeslice);
			sum_GiiDn += gDn(site, site, timeslice);
			sum_GiiUpDn += gUp(site, site, timeslice) * gDn(site, site, timeslice);
			//use nearest neighbor elements of Green functions:
			for (unsigned neighIndex = 0; neighIndex < latticeCoordination; ++neighIndex) {
				unsigned neigh = nearestNeigbors(neighIndex, site);
				sum_GneighUp += gUp(site, neigh, timeslice);
				sum_GneighDn += gDn(site, neigh, timeslice);
			}
		}
	}
	occUp = 1.0 - (1.0 / (N*m)) * sum_GiiUp;
	occDn = 1.0 - (1.0 / (N*m)) * sum_GiiDn;
	occTotal = occUp + occDn;
	ePotential = (U / (N*m)) * (sum_GiiUpDn + 0.5 * sum_GiiUp + 0.5 * sum_GiiDn);
	//Note: chemical potential term included in kinetic energy:
	eKinetic   = (t / (N*m)) * (sum_GneighUp + sum_GneighDn) + mu * occTotal;
	eTotal = eKinetic + ePotential;
}


unsigned DetHubbard::getNumberOfObservables() const {
	return obsCount;
}

num DetHubbard::obsNormalized(unsigned obsIndex) const {
	if (obsIndex < obsCount) {
		return *(obsValPointers[obsIndex]);
	} else {
		throw WrongObsIndex(obsIndex);
	}
}

std::string DetHubbard::getObservableName(unsigned obsIndex) const {
	if (obsIndex < obsCount) {
		return obsNames[obsIndex];
	} else {
		throw WrongObsIndex(obsIndex);
	}
}

std::string DetHubbard::getObservableShort(unsigned obsIndex) const {
	if (obsIndex < obsCount) {
		return obsShorts[obsIndex];
	} else {
		throw WrongObsIndex(obsIndex);
	}
}



inline unsigned DetHubbard::coordsToSite(const std::vector<unsigned>& coords) const {
    const unsigned dimensions = coords.size();
    int site = 0;
    for (unsigned dim = 0; dim < dimensions; ++dim) {
        site += coords[dim] * uint_pow(L, dim);
    }
    return site;
}

void DetHubbard::createNeighborTable() {
	using std::floor;
	nearestNeigbors.resize(2*d, N);
    std::vector<unsigned> curCoords(d);     //holds the x, y, z coordinate components of the current site
    std::vector<unsigned> newCoords(d);     //newly calculated coords of the neighbor
    for (unsigned site = 0; site < N; ++site) {
        int reducedSite = site;
        for (int dim = d - 1; dim >= 0; --dim) {
            curCoords[dim] = floor(reducedSite / uint_pow(L, dim));
            reducedSite -= curCoords[dim] * uint_pow(L, dim);
        }
        assert(reducedSite == 0);
        for (unsigned dim = 0; dim < d; ++dim) {
        	//neighbor in + direction, periodic
        	newCoords = curCoords;
        	newCoords[dim] = (newCoords[dim] + 1) % L;
        	nearestNeigbors(dim * 2, site) = coordsToSite(newCoords);
        	//neighbor in - direction, periodic
        	newCoords = curCoords;
        	newCoords[dim] = (newCoords[dim] - 1 + L) % L;
        	nearestNeigbors(dim * 2 + 1, site) = coordsToSite(newCoords);
        }
    }
}


void DetHubbard::setupRandomAuxfield() {
	for (unsigned timeslice = 0; timeslice < m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			if (rng.rand01() <= 0.5) {
				auxfield(site, timeslice) = +1;
			} else {
				auxfield(site, timeslice) = -1;
			}
		}
	}
}

void DetHubbard::setupTmat() {
	tmat = -mu * arma::eye(N, N);

	for (unsigned site = 0; site < N; ++site) {
		for (auto p = nearestNeigbors.begin_col(site);
				 p != nearestNeigbors.end_col(site); ++p) {
			tmat(*p, site) -= t;
		}
	}

	proptmat = computePropagator(dtau, tmat);
}


nummat DetHubbard::computePropagator(num scalar, nummat matrix) {
	using namespace arma;

	numvec eigval;
	nummat eigvec;
	eig_sym(eigval, eigvec, matrix);

	return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}


inline nummat DetHubbard::computeBmat(unsigned n2, unsigned n1, Spin spinz,
		const intmat& arbitraryAuxfield) {
	using namespace arma;

	assert(n2 > n1);
	assert(n2 <= m);
	//assert(n1 >= 0);

	int sign = static_cast<int>(spinz);

	nummat B = diagmat(exp(sign * alpha * arbitraryAuxfield.col(n2 - 1))) * proptmat;

	for (unsigned n = n2 - 1; n >= n1 + 1; --n) {
		B *= diagmat(exp(sign * alpha * arbitraryAuxfield.col(n - 1))) * proptmat;
	}

	return B;
}

inline nummat DetHubbard::computeBmat(unsigned n2, unsigned n1, Spin spinz) {
	return computeBmat(n2, n1, spinz, auxfield);
}

inline nummat DetHubbard::computeGreenFunction(
		const nummat& bTau0, const nummat& bBetaTau) {
	return arma::inv(arma::eye(N,N) + bTau0 * bBetaTau);
}

inline nummat DetHubbard::computeGreenFunction(unsigned timeslice,
		Spin spinz) {
	//TODO: should use stored B-matrices, for the timeslices that have not changed
	return computeGreenFunction(computeBmat(timeslice, 0, spinz),
			computeBmat(m, timeslice, spinz));
}


void DetHubbard::computeAllGreenFunctions() {
	auto compute = [this](Spin spinz, numcube& green) {
		for (unsigned timeslice = 0; timeslice < m; ++timeslice) {
			green.slice(timeslice) = computeGreenFunction(timeslice, spinz);
		}
	};
	compute(Spin::Up,   gUp);
	compute(Spin::Down, gDn);
}

num DetHubbard::weightRatioGeneric(const intmat& auxfieldBefore,
		const intmat& auxfieldAfter) {
	using namespace arma;
	return det(eye(N,N) + computeBmat(m, 0, Spin::Up, auxfieldAfter)) *
		   det(eye(N,N) + computeBmat(m, 0, Spin::Down, auxfieldAfter)) /
			(det(eye(N,N) + computeBmat(m, 0, Spin::Up, auxfieldBefore)) *
			 det(eye(N,N) + computeBmat(m, 0, Spin::Down, auxfieldBefore)));
}

inline num DetHubbard::weightRatioSingleFlip(unsigned site, unsigned timeslice) {
	//TODO: possibly precompute the exponential factors (auxfield is either +/- 1), would require an if though.
	return (1 + std::exp(-2 * alpha * auxfield(timeslice, site)) *
			(1 - gUp(site,site,timeslice)))
			*
		   (1 + std::exp(+2 * alpha * auxfield(timeslice, site)) *
			(1 - gDn(site,site,timeslice)));
}

inline void DetHubbard::updateGreenFunctionsAfterFlip(unsigned site, unsigned timeslice) {
	auto update = [this, site](nummat& green, num expfactor) {
		const nummat& greenOld = green;		//reference
		nummat greenNew = green;			//copy
		const nummat oneMinusGreenOld = arma::eye(N,N) - greenOld;
		num divisor = 1 + expfactor * (1 - greenOld(site, site));
		for (unsigned y = 0; y < N; ++y) {
			for (unsigned x = 0; x < N; ++x) {
				greenNew(x, y) -=  greenOld(x, site) * expfactor *
						oneMinusGreenOld(site, y) / divisor;
			}
		}
	};

	update(gUp.slice(timeslice), std::exp(-2 * alpha * auxfield(timeslice, site)));
	update(gDn.slice(timeslice), std::exp(+2 * alpha * auxfield(timeslice, site)));
}












