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



nummat DetHubbard::computeBmat(unsigned n2, unsigned n1, DetHubbard::Spin spinComponent) {
	using namespace arma;

	assert(n2 > n1);
	assert(n2 <= m);
	assert(n1 >= 0);

	int sign = static_cast<int>(spinComponent);

	nummat B = diagmat(exp(sign * alpha * auxfield.col(n2 - 1))) * proptmat;

	for (unsigned n = n2 - 1; n >= n1 + 1; --n) {
		B *= diagmat(exp(sign * alpha * auxfield.col(n - 1))) * proptmat;
	}

	return B;
}



