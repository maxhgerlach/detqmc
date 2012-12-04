/*
 * dethubbard.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: gerlach
 */

#include <vector>
#include <cmath>
#include <cassert>
#include "dethubbard.h"

using namespace std;

DetHubbard::DetHubbard(num t, num U, num mu, unsigned L, unsigned d, num beta, unsigned m) :
		t(t), U(U), mu(mu), L(L), d(d), N(static_cast<unsigned>(pow(L,d))),
		beta(beta), m(m), dtau(beta/m),
		alpha(acosh(exp(dtau * U * 0.5))),
		nearestNeigbors(2*d, N),			//coordination number: 2*d
		tmat(N, N), proptmat(N,N),
		auxfield(N, m)    //m columns of N rows
{
	createNeighborTable();
	setupTmat();
}

DetHubbard::~DetHubbard() {
	// TODO Auto-generated destructor stub
}

inline unsigned DetHubbard::coordsToSite(const std::vector<unsigned>& coords) {
    const unsigned dimensions = coords.size();
    int site = 0;
    for (unsigned dim = 0; dim < dimensions; ++dim) {
        site += coords[dim] * (unsigned)pow(L, dim);
    }
    return site;
}

void DetHubbard::createNeighborTable() {
	nearestNeigbors.resize(2*d, N);
    std::vector<unsigned> curCoords(d);     //holds the x, y, z coordinate components of the current site
    std::vector<unsigned> newCoords(d);     //newly calculated coords of the neighbor
    for (unsigned site = 0; site < N; ++site) {
        int reducedSite = site;
        for (int dim = d - 1; dim >= 0; --dim) {
            curCoords[dim] = floor(reducedSite / (unsigned)pow(L, dim));
            reducedSite -= curCoords[dim] * (unsigned)pow(L, dim);
        }
        assert(reducedSite == 0);
        for (unsigned dim = 0; dim < d; ++dim) {
        	//neighbor in + direction, periodic
        	newCoords = curCoords;
        	newCoords[dim] = (newCoords[dim] + 1) % L;
        	nearestNeigbors(dim * 2, N) = coordsToSite(newCoords);
        	//neighbor in - direction, periodic
        	newCoords = curCoords;
        	newCoords[dim] = (newCoords[dim] - 1 + L) % L;
        	nearestNeigbors(dim * 2 + 1, N) = coordsToSite(newCoords);
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
		const nummat& arbitraryAuxfield) {
	using namespace arma;

	assert(n2 > n1);
	assert(n2 <= m);
	assert(n1 >= 0);

	int sign = static_cast<int>(spinz);

	nummat B = diagmat(exp(sign * alpha * arbitraryAuxfield.col(n2 - 1))) * proptmat;

	for (unsigned n = n2 - 1; n >= n1 + 1; --n) {
		B *= diagmat(exp(sign * alpha * arbitraryAuxfield.col(n - 1))) * proptmat;
	}

	return B;
}

inline nummat DetHubbard::computeBmat(unsigned n2, unsigned n1, Spin spinz) {
	return computeBmat(n2, n1, spinz);
}

inline nummat DetHubbard::computeGreenFunction(
		const nummat& bTau0, const nummat& bBetaTau) {
	return arma::inv(arma::eye(N,N) + bTau0 * bBetaTau);
}

inline nummat DetHubbard::computeGreenFunction(unsigned timeslice,
		Spin spinz) {
	return computeGreenFunction(computeBmat(timeslice, 0, spinz),
			computeBmat(m, timeslice, spinz));
}


void DetHubbard::computeAllGreenFunctions() {
	auto compute = [this](Spin spinz, numcube& green) {
		for (unsigned timeslice = 0; timeslice < m; ++timeslice) {
			green.slice(timeslice) = computeGreenFunction(timeslice, spinz);
		}
	};
	compute(Spin::Up,   greenfctUp);
	compute(Spin::Down, greenfctDown);
}

num DetHubbard::weightRatioGeneric(const nummat& auxfieldBefore,
		const nummat& auxfieldAfter) {
	using namespace arma;
	return det(eye(N,N) + computeBmat(m, 0, Spin::Up, auxfieldAfter)) *
		   det(eye(N,N) + computeBmat(m, 0, Spin::Down, auxfieldAfter)) /
			(det(eye(N,N) + computeBmat(m, 0, Spin::Up, auxfieldBefore)) *
			 det(eye(N,N) + computeBmat(m, 0, Spin::Down, auxfieldBefore)));
}

inline num DetHubbard::weightRatioSingleFlip(unsigned site, unsigned timeslice) {
	//TODO: possibly precompute the exponential factors (auxfield is either +/- 1), would require an if though.
	return (1 + exp(-2 * alpha * auxfield(timeslice, site)) *
			(1 - greenfctUp(site,site,timeslice)))
			*
		   (1 + exp(+2 * alpha * auxfield(timeslice, site)) *
			(1 - greenfctDown(site,site,timeslice)));
}

inline void DetHubbard::updateGreenFunctionsAfterFlip(unsigned site, unsigned timeslice) {
	auto update = [N, site](nummat& green, num expfactor) {
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

	update(greenfctUp.slice(timeslice),   exp(-2 * alpha * auxfield(timeslice, site)));
	update(greenfctDown.slice(timeslice), exp(+2 * alpha * auxfield(timeslice, site)));
}










