/*
 * detsdw.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include <cmath>
#include <numeric>
#include "detsdw.h"

const num PhiLow = 0.0;
const num PhiHigh = 2.0;

std::unique_ptr<DetSDW> createDetSDW(RngWrapper& rng, ModelParams pars) {
	//TODO: add checks
	return std::unique_ptr<DetSDW>(new DetSDW(rng, pars));
}

DetSDW::DetSDW(RngWrapper& rng_, const ModelParams& pars) :
		DetModelGC<1,cpx>(pars, L*L),
		rng(rng_),
		L(pars.L), N(L*L), r(pars.r), mu(pars.mu), c(1), u(1), lambda(1), //TODO: make these controllable by parameter
		spaceNeigh(L), timeNeigh(m),
		propK(), propKx(propK[XBAND]), propKy(propK[YBAND]),
		g(green[0]), gFwd(greenFwd[0]), gBwd(greenBwd[0]),
		phi1(N, m+1), phi2(N, m+1), phi3(N, m+1), phiCosh(N, m+1), phiSinh(N, m+1)
{
	setupRandomPhi();
	setupPropK();
	computeBmat[0] = [this](unsigned k2, unsigned k1) {
		return this->computeBmatSDW(k2, k1);
	};
	setupUdVStorage();
}

DetSDW::~DetSDW() {
}

unsigned DetSDW::getSystemN() const {
	return N;
}

MetadataMap DetSDW::prepareModelMetadataMap() const {
	MetadataMap meta;
	meta["model"] = "sdw";
	meta["timedisplaced"] = (timedisplaced ? "true" : "false");
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
	META_INSERT(r);
	META_INSERT(mu);
	META_INSERT(L);
	META_INSERT(d);
	META_INSERT(N);
	META_INSERT(beta);
	META_INSERT(m);
	META_INSERT(dtau);
	META_INSERT(s);
#undef META_INSERT
	return meta;
}

void DetSDW::measure() {
}

void DetSDW::setupRandomPhi() {
	for_each_timeslice( [this](unsigned k) {
		for_each_site( [this, k](unsigned site) {
			phi1(site, k) = rng.randRange(PhiLow, PhiHigh);
			phi2(site, k) = rng.randRange(PhiLow, PhiHigh);
			phi3(site, k) = rng.randRange(PhiLow, PhiHigh);
			num phiNorm = std::sqrt(std::pow(phi1(site, k), 2)
			+ std::pow(phi2(site, k), 2)
			+ std::pow(phi3(site, k), 2));
			phiCosh(site, k) = std::cosh(phiNorm);
			phiSinh(site, k) = std::sinh(phiNorm) / phiNorm;
		} );
	} );
}

void DetSDW::setupPropK() {
	std::array<std::array<num,z>, 2> t;
	t[XBAND][XPLUS] = t[XBAND][XMINUS] = -1.0;
	t[XBAND][YPLUS] = t[XBAND][YMINUS] =  0.5;
	t[YBAND][XPLUS] = t[YBAND][XMINUS] = -0.5;
	t[YBAND][YPLUS] = t[YBAND][YMINUS] =  1.0;

	for_each_band( [this, &t](unsigned band) {
		MatNum k = -mu * arma::eye(N,N);
		for_each_site( [this, band, &k, &t](unsigned site) {
			for (unsigned dir = 0; dir < z; ++dir) {
				unsigned neigh = spaceNeigh(dir, site);
				k(site, neigh) -= t[band][dir];
			}
		} );
		propK[band] = computePropagator(dtau, k);
	} );
}

MatCpx DetSDW::computeBmatSDW(unsigned k2, unsigned k1) const {
	using arma::eye; using arma::zeros; using arma::diagmat;
	if (k2 == k1) {
		return MatCpx(eye(4*N,4*N), zeros(4*N,4*N));
	}
	assert(k2 > k1);
	assert(k2 <= m);


	auto singleTimesliceProp = [this, N](unsigned k) {
		MatCpx result(4*N, 4*N);

		//submatrix view helper for a 4N*4N matrix
		auto block = [&result, N](unsigned row, unsigned col) {
			return result.submat(row * N, col * N,
					             row * (N + 1) - 1, col * (N + 1) - 1);
		};
		auto& kphi1 = phi1.col(k);
		auto& kphi2 = phi2.col(k);
		auto& kphi3 = phi3.col(k);
		auto& kphiCosh = phiCosh.col(k);
		auto& kphiSinh = phiSinh.col(k);
		//TODO: is this the best way to set the real and imaginary parts of a complex submatrix?
		block(0, 0) = MatCpx(diagmat(kphiCosh) * propKx, zeros(N,N));
		block(0, 1).zeros();
		block(0, 2) = MatCpx(diagmat(+kphi3 * kphiSinh) * propKy,
				zeros(N,N));
		block(0, 3) = MatCpx(diagmat( kphi1 * kphiSinh) * propKy,
				diagmat(-kphi2 * kphiSinh) * propKy);
		block(1, 0).zeros();
		block(1, 1) = block(0, 0);
		block(1, 2) = MatCpx(diagmat( kphi1 * kphiSinh) * propKy,
				diagmat(+kphi2 * kphiSinh) * propKy);
		block(1, 3) = MatCpx(diagmat(-kphi3 * kphiSinh) * propKy,
				zeros(N,N));
		block(2, 0) = MatCpx(diagmat(+kphi3 * kphiSinh) * propKx,
				zeros(N,N));
		block(2, 1) = MatCpx(diagmat( kphi1 * kphiSinh) * propKx,
				diagmat(-kphi2 * kphiSinh) * propKx);
		block(2, 2) = block(0, 0);
		block(2, 3).zeros();
		block(3, 0) = MatCpx(diagmat( kphi1 * kphiSinh) * propKx,
				diagmat(+kphi2 * kphiSinh) * propKx);
		block(3, 1) = MatCpx(diagmat(-kphi3 * kphiSinh) * propKy,
				zeros(N,N));
		block(3, 2).zeros();
		block(3, 3) = block(0, 0);

		return result;
	};

	MatCpx result = singleTimesliceProp(k2);

	for (unsigned k = k2 - 1; k > k1; --k) {
		result *= singleTimesliceProp(k);
	}

	return result;
}

void DetSDW::updateInSlice(unsigned timeslice) {
	for_each_site( [this, timeslice](unsigned site) {
		Phi newphi = proposeNewField(site, timeslice);
	});
}

DetSDW::Phi DetSDW::proposeNewField(unsigned site, unsigned timeslice) {
	(void) site; (void) timeslice;
	//TODO: make this smarter!
	return Phi{rng.randRange(PhiLow, PhiHigh),
			   rng.randRange(PhiLow, PhiHigh),
			   rng.randRange(PhiLow, PhiHigh)};
}

num DetSDW::deltaSPhi(unsigned site, unsigned timeslice, const Phi newphi) {
	const Phi oldphi = {phi1(site, timeslice),
					    phi2(site, timeslice),
					    phi3(site, timeslice)};

	num oldphiSq = arma::dot(oldphi, oldphi);
	num newphiSq = arma::dot(newphi, newphi);
	num delta1 = (2.0 / (8.0 * dtau*dtau * c*c) + 0.5*z + 0.5*r) * (newphiSq - oldphiSq);

	unsigned kEarlier = timeNeigh(ChainDir::MINUS, timeNeigh(ChainDir::MINUS, timeslice));
	Phi phiEarlier = Phi{phi1(site, kEarlier),
					     phi2(site, kEarlier),
					     phi3(site, kEarlier)};
	unsigned kLater = timeNeigh(ChainDir::PLUS, timeNeigh(ChainDir::PLUS, timeslice));
	Phi phiLater = Phi{phi1(site, kLater),
					   phi2(site, kLater),
					   phi3(site, kLater)};
	Phi phiTimeNeigh  = (1.0 / (8.0 * dtau*dtau * c*c)) * (phiEarlier + phiLater);
	Phi phiSpaceNeigh = 0.5 * std::accumulate(spaceNeigh.beginNeighbors(site),
								              spaceNeigh.endNeighbors(site),
								              Phi{0,0,0},
								              [this, timeslice] (Phi accum, unsigned neighSite) {
													return accum + Phi{phi1(neighSite, timeslice),
															     	   phi2(neighSite, timeslice),
															           phi3(neighSite, timeslice)};
	                                           }
											 );
	num delta2 = -2.0 * arma::dot((phiTimeNeigh + phiSpaceNeigh), (newphi - oldphi));

	num oldphiPow4 = oldphiSq * oldphiSq;
	num newphiPow4 = newphiSq * newphiSq;
	num delta3 = 0.25*u * (newphiPow4 - oldphiPow4);

	return delta1 + delta2 + delta3;
}





