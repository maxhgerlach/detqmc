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


std::unique_ptr<DetHubbard> createDetHubbard(RngWrapper& rng, const ModelParams& pars) {
	//check parameters
	using namespace boost::assign;
	std::vector<std::string> neededModelPars;
	neededModelPars += "t", "U", "mu", "L", "d", "beta", "m";
	for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
		if (pars.specified.count(*p) == 0) {
			throw ParameterMissing(*p);
		}
	}

	return std::unique_ptr<DetHubbard>(new DetHubbard(rng,
			pars.t, pars.U, pars.mu, pars.L, pars.d, pars.beta, pars.m)	);
}



DetHubbard::DetHubbard(RngWrapper& rng_,
		num t, num U, num mu, unsigned L, unsigned d, num beta, unsigned m) :
		rng(rng_),
		t(t), U(U), mu(mu), L(L), d(d),
		z(2*d), //coordination number: 2*d
		N(static_cast<unsigned>(uint_pow(L,d))),
		beta(beta), m(m), dtau(beta/m),
		alpha(acosh(std::exp(dtau * U * 0.5))),
		nearestNeigbors(z, N),
		proptmat(N,N),
		auxfield(N, m),              //m columns of N rows
		gUp(N,N,m), gDn(N,N,m),      //m slices of N columns x N rows
		obsNames(), obsShorts(), obsValPointers(), obsCount(0)
{
	createNeighborTable();
	setupRandomAuxfield();
//	auxfield.print(std::cout);
	setupPropTmat();
//	tmat.print(std::cout);
	using namespace boost::assign;         // bring operator+=() into scope
	obsNames += "occupationUp", "occupationDown", "totalOccupation",
			"doubleOccupation", "localMoment",
			"kineticEnergy", "potentialEnergy", "totalEnergy";
	obsShorts += "nUp", "nDown", "n", "n2", "m^2", "e_t", "e_U", "e";
	obsValPointers += &occUp, &occDn, &occTotal, &occDouble, &localMoment,
			&eKinetic, &ePotential, &eTotal;
	assert(obsNames.size() == obsShorts.size());
	assert(obsNames.size() == obsValPointers.size());
	obsCount = obsNames.size();
}

DetHubbard::~DetHubbard() {
}


MetadataMap DetHubbard::prepareModelMetadataMap() {
	MetadataMap meta;
	meta["model"] = "hubbard";
#define META_INSERT(VAR) meta[#VAR] = numToString(VAR)
	META_INSERT(t);
	META_INSERT(U);
	META_INSERT(mu);
	META_INSERT(L);
	META_INSERT(d);
	META_INSERT(N);
	META_INSERT(beta);
	META_INSERT(m);
	META_INSERT(dtau);
#undef META_INSERT
	return meta;
}

void DetHubbard::sweepSimple() {
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		gUp.slice(timeslice-1) = computeGreenFunctionNaive(timeslice, Spin::Up);
		gDn.slice(timeslice-1) = computeGreenFunctionNaive(timeslice, Spin::Down);
		// picking sites linearly: system seemed to alternate between two configurations
		// sweep after sweep
//		for (unsigned site = 0; site < N; ++site) {
		for (unsigned count = 0; count < N; ++count) {
			unsigned site = rng.randInt(0, N-1);

			num ratio =  weightRatioSingleFlip(site, timeslice);

			//DEBUG: comparison of weight ratio calculation
//			intmat newAuxfield = auxfield;
//			newAuxfield(site, timeslice - 1) *= -1;
//			num refRatio = weightRatioGeneric(auxfield, newAuxfield);
//			std::cout << ratio << " vs. " << refRatio << std::endl;

			assert(ratio > 0);

			//Metropolis
			if (ratio > 1 or rng.rand01() < ratio) {
//			if (refRatio > 1 or rng.rand01() < refRatio) {
//				std::cout << "acc" << '\n';
//				auxfield = newAuxfield;
				auxfield(site, timeslice-1) *= -1;
				updateGreenFunctionsAfterFlip(site, timeslice);
			}
//			} else {
//				std::cout << "no acc" << '\n';
//			}
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
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			//use diagonal elements of Green functions:
			sum_GiiUp += gUp(site, site, timeslice - 1);
			sum_GiiDn += gDn(site, site, timeslice - 1);
			sum_GiiUpDn += gUp(site, site, timeslice - 1) * gDn(site, site, timeslice - 1);
			//use nearest neighbor elements of Green functions:
			for (unsigned neighIndex = 0; neighIndex < z; ++neighIndex) {
				unsigned neigh = nearestNeigbors(neighIndex, site);
				sum_GneighUp += gUp(site, neigh, timeslice - 1);
				sum_GneighDn += gDn(site, neigh, timeslice - 1);
			}
		}
	}
	occUp = 1.0 - (1.0 / (N*m)) * sum_GiiUp;
	occDn = 1.0 - (1.0 / (N*m)) * sum_GiiDn;
	occTotal = occUp + occDn;

	occDouble = 1.0 + (1.0 / (N*m)) * (sum_GiiUpDn - sum_GiiUp - sum_GiiDn);
	localMoment = occTotal - 2*occDouble;

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
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			if (rng.rand01() <= 0.5) {
				auxfield(site, timeslice - 1) = +1;
			} else {
				auxfield(site, timeslice - 1) = -1;
			}
		}
	}
}

void DetHubbard::setupPropTmat() {
	nummat tmat = -mu * arma::eye(N, N);

	for (unsigned site = 0; site < N; ++site) {
		//hopping between nearest neighbors
		for (auto p = nearestNeigbors.begin_col(site);
				 p != nearestNeigbors.end_col(site); ++p) {
			tmat(*p, site) -= t;
		}
	}

	proptmat = computePropagator(dtau, tmat);
}

void DetHubbard::setupUdVStorage() {
	auto setup = [this](std::vector<UdV>& storage, Spin spinz) {
		auto svd = [](const nummat& mat) -> UdV {
			UdV res;
			arma::svd(res.U, res.d, res.V, mat, "standard");
			return res;
		};

		storage = std::vector<UdV>(m + 1);

		storage[0].U = arma::eye(N,N);
		storage[0].d = arma::ones(N);
		storage[0].V = arma::eye(N,N);
		storage[1] = svd(computeBmatNaive(1, 0, spinz));

		for (unsigned k = 1; k <= m - 1; ++k) {
			const nummat& U_k = storage[k].U;
			const nummat& d_k = storage[k].d;
			const nummat& V_k = storage[k].V;
			nummat B_kp1 = computeBmatNaive(k + 1, k, spinz);
			UdV UdV_temp = svd((B_kp1 * U_k) * arma::diagmat(d_k));
			storage[k+1].U = UdV_temp.U;
			storage[k+1].d = UdV_temp.d;
			storage[k+1].V = UdV_temp.V * V_k;
		}
	};

	setup(UdVStorageUp, Spin::Up);
	setup(UdVStorageDn, Spin::Down);
}


nummat DetHubbard::computePropagator(num scalar, nummat matrix) {
	using namespace arma;

	numvec eigval;
	nummat eigvec;
	eig_sym(eigval, eigvec, matrix);

	return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}


inline nummat DetHubbard::computeBmatNaive(unsigned n2, unsigned n1, Spin spinz,
		const intmat& arbitraryAuxfield) {
	using namespace arma;

	if (n2 == n1) {
		return eye(N, N);
	}

	assert(n2 > n1);
	assert(n2 <= m);
	//assert(n1 >= 0);

	int sign = static_cast<int>(spinz);

	//Propagator using the HS-field potential for the given timeslice
	auto singleTimeslicePropagator = [this, sign, arbitraryAuxfield](unsigned timeslice) -> nummat {
		return diagmat(exp(sign * alpha * arbitraryAuxfield.col(timeslice-1))) * proptmat;
	};

	nummat B = singleTimeslicePropagator(n2);

	for (unsigned n = n2 - 1; n >= n1 + 1; --n) {
		B *= singleTimeslicePropagator(n);
	}

	return B;
}

inline nummat DetHubbard::computeBmatNaive(unsigned n2, unsigned n1, Spin spinz) {
	return computeBmatNaive(n2, n1, spinz, auxfield);
}

inline nummat DetHubbard::computeGreenFunctionNaive(
		const nummat& bTau0, const nummat& bBetaTau) {
	return arma::inv(arma::eye(N,N) + bTau0 * bBetaTau);
}

inline nummat DetHubbard::computeGreenFunctionNaive(unsigned timeslice,
		Spin spinz) {
	//TODO: should use stored B-matrices, for the timeslices that have not changed
	return computeGreenFunctionNaive(computeBmatNaive(timeslice, 0, spinz),
			                    computeBmatNaive(m, timeslice, spinz));
}


num DetHubbard::weightRatioGenericNaive(const intmat& auxfieldBefore,
		const intmat& auxfieldAfter) {
	using namespace arma;

	num weightAfterUp  = det(eye(N,N) + computeBmatNaive(m, 0, Spin::Up, auxfieldAfter));
	num weightBeforeUp = det(eye(N,N) + computeBmatNaive(m, 0, Spin::Up, auxfieldBefore));
	num ratioUp = weightAfterUp / weightBeforeUp;

	num weightAfterDown  = det(eye(N,N) + computeBmatNaive(m, 0, Spin::Down, auxfieldAfter));
	num weightBeforeDown = det(eye(N,N) + computeBmatNaive(m, 0, Spin::Down, auxfieldBefore));
	num ratioDown = weightAfterDown / weightBeforeDown;

//	std::cout << weightafterup << " " << weightbeforeup << " " << weightafterdown << " " << weightbeforedown << "\n";

	return ratioUp * ratioDown;
}

inline num DetHubbard::weightRatioSingleFlip(unsigned site, unsigned timeslice) {
	using std::exp;
	//TODO: possibly precompute the exponential factors (auxfield is either +/- 1), would require an if though.

	//exponential factors
	num expUp   = exp(-2 * alpha * auxfield(site, timeslice-1));
	num expDown = exp( 2 * alpha * auxfield(site, timeslice-1));

	num ratioUp   = 1 + (expUp   - 1) * (1 - gUp(site,site, timeslice-1));
	num ratioDown = 1 + (expDown - 1) * (1 - gDn(site,site, timeslice-1));

//	std::cout << ratioUp << " " << ratioDown << std::endl;

	return ratioUp * ratioDown;
}

inline void DetHubbard::updateGreenFunctionsAfterFlip(unsigned site, unsigned timeslice) {
	auto update = [this, site](nummat& green, num expfactor) {
		const nummat& greenOld = green;		//reference
		nummat greenNew = green;			//copy
		const nummat oneMinusGreenOld = arma::eye(N,N) - greenOld;
		num divisor = 1 + (expfactor - 1) * (1 - greenOld(site, site));
		for (unsigned y = 0; y < N; ++y) {
			for (unsigned x = 0; x < N; ++x) {
				greenNew(x, y) -=  greenOld(x, site) * (1 - expfactor) *
						oneMinusGreenOld(site, y) / divisor;
			}
		}
		green = greenNew;
	};

	using std::exp;
	update(gUp.slice(timeslice-1), exp(-2 * alpha * auxfield(site, timeslice-1)));
	update(gDn.slice(timeslice-1), exp(+2 * alpha * auxfield(site, timeslice-1)));
}



