/*
 * detsdw.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include <cmath>
#include <numeric>
#include <functional>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include <boost/assign/std/vector.hpp>    // 'operator+=()' for vectors
#pragma GCC diagnostic warning "-Wconversion"
#pragma GCC diagnostic warning "-Weffc++"
#include "observable.h"
#include "detsdw.h"

const num PhiLow = 0.0;
const num PhiHigh = 2.0;

std::unique_ptr<DetSDW> createDetSDW(RngWrapper& rng, ModelParams pars) {
	//TODO: add checks
	pars = updateTemperatureParameters(pars);
	return std::unique_ptr<DetSDW>(new DetSDW(rng, pars));
}

DetSDW::DetSDW(RngWrapper& rng_, const ModelParams& pars) :
		DetModelGC<1,cpx>(pars, 4 * pars.L*pars.L),
		rng(rng_),
		L(pars.L), N(L*L), r(pars.r), mu(pars.mu), c(1), u(1), lambda(1), //TODO: make these controllable by parameter
		spaceNeigh(L), timeNeigh(m),
		propK(), propKx(propK[XBAND]), propKy(propK[YBAND]),
		g(green[0]), gFwd(greenFwd[0]), gBwd(greenBwd[0]),
		phi0(N, m+1), phi1(N, m+1), phi2(N, m+1), phiCosh(N, m+1), phiSinh(N, m+1),
		normPhi()
{
	g = CubeCpx(4*N,4*N, m+1);
	if (pars.timedisplaced) {
		gFwd = CubeCpx(4*N,4*N, m+1);
		gBwd = CubeCpx(4*N,4*N, m+1);
	}
	setupRandomPhi();
	setupPropK();
	computeBmat[0] = [this](unsigned k2, unsigned k1) {
		return this->computeBmatSDW(k2, k1);
	};
	setupUdVStorage();

	using std::cref;
	using namespace boost::assign;
	obsScalar += ScalarObservable(cref(normPhi), "normPhi", "np");
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
	Phi meanPhi;
	meanPhi[0] = averageWholeSystem(phi0, 0.0);
	meanPhi[1] = averageWholeSystem(phi1, 0.0);
	meanPhi[2] = averageWholeSystem(phi2, 0.0);
	normPhi = arma::norm(meanPhi, 2);
}

void DetSDW::setupRandomPhi() {
	for_each_timeslice( [this](unsigned k) {
		for_each_site( [this, k](unsigned site) {
			phi0(site, k) = rng.randRange(PhiLow, PhiHigh);
			phi1(site, k) = rng.randRange(PhiLow, PhiHigh);
			phi2(site, k) = rng.randRange(PhiLow, PhiHigh);
			num phiNorm = std::sqrt(std::pow(phi0(site, k), 2)
			+ std::pow(phi1(site, k), 2)
			+ std::pow(phi2(site, k), 2));
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
					             (row + 1) * N - 1, (col + 1) * N - 1);
		};
		auto& kphi0 = phi0.col(k);
		auto& kphi1 = phi1.col(k);
		auto& kphi2 = phi2.col(k);
		auto& kphiCosh = phiCosh.col(k);
		auto& kphiSinh = phiSinh.col(k);
		//TODO: is this the best way to set the real and imaginary parts of a complex submatrix?
		block(0, 0) = MatCpx(diagmat(kphiCosh) * propKx, zeros(N,N));
		block(0, 1).zeros();
		block(0, 2) = MatCpx(diagmat(+kphi2 % kphiSinh) * propKy,
				zeros(N,N));
		block(0, 3) = MatCpx(diagmat( kphi0 % kphiSinh) * propKy,
				diagmat(-kphi1 % kphiSinh) * propKy);
		block(1, 0).zeros();
		block(1, 1) = block(0, 0);
		block(1, 2) = MatCpx(diagmat( kphi0 % kphiSinh) * propKy,
				diagmat(+kphi1 % kphiSinh) * propKy);
		block(1, 3) = MatCpx(diagmat(-kphi2 % kphiSinh) * propKy,
				zeros(N,N));
		block(2, 0) = MatCpx(diagmat(+kphi2 % kphiSinh) * propKx,
				zeros(N,N));
		block(2, 1) = MatCpx(diagmat( kphi0 % kphiSinh) * propKx,
				diagmat(-kphi1 % kphiSinh) * propKx);
		block(2, 2) = block(0, 0);
		block(2, 3).zeros();
		block(3, 0) = MatCpx(diagmat( kphi0 % kphiSinh) * propKx,
				diagmat(+kphi1 % kphiSinh) * propKx);
		block(3, 1) = MatCpx(diagmat(-kphi2 % kphiSinh) * propKy,
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

		num propSPhi = std::exp(-deltaSPhi(site, timeslice, newphi));

		//delta = e^(V_new)*e^(-V_old) - 1

		//compute non-zero elements of delta
		auto evMatrix = [](int sign, num kphi0, num kphi1,
						   num kphi2, num kphiCosh, num kphiSinh) -> MatCpx::fixed<4,4> {
			MatNum::fixed<4,4> ev_real;
			ev_real.diag().fill(kphiCosh);
			ev_real(0,1) = ev_real(1,0) = ev_real(2,3) = ev_real(3,2) = 0;
			ev_real(2,0) = ev_real(0,2) =  sign * kphi2 * kphiSinh;
			ev_real(2,1) = ev_real(0,3) =  sign * kphi0 * kphiSinh;
			ev_real(3,0) = ev_real(1,2) = -sign * kphi1 * kphiSinh;
			ev_real(3,1) = ev_real(1,3) = -sign * kphi2 * kphiSinh;

			MatCpx::fixed<4,4> ev;
			ev.set_real(ev_real);
			ev(0,3).imag(-sign * kphi1 * kphiSinh);
			ev(1,2).imag( sign * kphi1 * kphiSinh);
			ev(2,1).imag(-sign * kphi1 * kphiSinh);
			ev(3,0).imag( sign * kphi1 * kphiSinh);

			return ev;
		};
		MatCpx::fixed<4,4> emv = evMatrix(
				-1,
				phi0(site, timeslice), phi1(site, timeslice), phi2(site, timeslice),
				phiCosh(site, timeslice), phiSinh(site, timeslice)
				);
		num normnewphi = arma::norm(newphi,2);
		MatCpx::fixed<4,4> env = evMatrix(
				+1,
				newphi[0], newphi[1], newphi[2],
				std::cosh(normnewphi),
				std::sinh(normnewphi) / normnewphi
				);
		MatCpx::fixed<4,4> deltanonzero = env * emv;
		deltanonzero.diag() -= cpx(1.0, 0);

		SpMatCpx delta(4*N, 4*N);
		arma::uvec::fixed<4> idx = {site, site + N, site + 2*N, site + 3*N};
		//Armadilo lacks non-contiguous submatrix views for sparse matrices
		unsigned i = 0;
		for (auto col: idx) {
			unsigned j = 0;
			for (auto row: idx) {
				delta(row, col) = deltanonzero(j, i);
				++j;
			}
			++i;
		}

		//TODO: inefficient!
		static MatCpx eyeCpx = MatCpx(arma::eye(4*N, 4*N), arma::zeros(4*N, 4*N));
		MatCpx target = eyeCpx + delta * (eyeCpx - g.slice(timeslice));

		cpx weightRatio = arma::det(target);
		std::cout << weightRatio << std::endl;
		num propSFermion = weightRatio.real();

		num prop = propSPhi * propSFermion;

		if (prop > 1.0 or rng.rand01() < prop) {
			phi0(site, timeslice) = newphi[0];
			phi1(site, timeslice) = newphi[1];
			phi2(site, timeslice) = newphi[2];

			g.slice(timeslice) *= arma::inv(target);
		}
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
	//calculation currently based on 3-point formula for discretized time derivative

	const Phi oldphi = {phi0(site, timeslice),
					    phi1(site, timeslice),
					    phi2(site, timeslice)};

	num oldphiSq = arma::dot(oldphi, oldphi);
	num newphiSq = arma::dot(newphi, newphi);
	num delta1 = (2.0 / (8.0 * dtau*dtau * c*c) + 0.5*z + 0.5*r) * (newphiSq - oldphiSq);

	unsigned kEarlier = timeNeigh(ChainDir::MINUS, timeNeigh(ChainDir::MINUS, timeslice));
	Phi phiEarlier = Phi{phi0(site, kEarlier),
					     phi1(site, kEarlier),
					     phi2(site, kEarlier)};
	unsigned kLater = timeNeigh(ChainDir::PLUS, timeNeigh(ChainDir::PLUS, timeslice));
	Phi phiLater = Phi{phi0(site, kLater),
					   phi1(site, kLater),
					   phi2(site, kLater)};
	Phi phiTimeNeigh  = (1.0 / (8.0 * dtau*dtau * c*c)) * (phiEarlier + phiLater);
	Phi phiSpaceNeigh = 0.5 * std::accumulate(spaceNeigh.beginNeighbors(site),
								              spaceNeigh.endNeighbors(site),
								              Phi{0,0,0},
								              [this, timeslice] (Phi accum, unsigned neighSite) {
													return accum + Phi{phi0(neighSite, timeslice),
															     	   phi1(neighSite, timeslice),
															           phi2(neighSite, timeslice)};
	                                           }
											 );
	num delta2 = -2.0 * arma::dot((phiTimeNeigh + phiSpaceNeigh), (newphi - oldphi));

	num oldphiPow4 = oldphiSq * oldphiSq;
	num newphiPow4 = newphiSq * newphiSq;
	num delta3 = 0.25*u * (newphiPow4 - oldphiPow4);

	return delta1 + delta2 + delta3;
}





