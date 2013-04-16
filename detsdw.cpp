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
#include "exceptions.h"
#include "timing.h"

//initial values for field components chosen from this range:
const num PhiLow = -1;
const num PhiHigh = 1;
//adjustment of phiDelta:
const num InitialPhiDelta = 0.5;
const unsigned int AccRatioAdjustmentSamples = 100;
const num phiDeltaGrowFactor = 1.01;
const num phiDeltaShrinkFactor = 0.99;


std::unique_ptr<DetSDW> createDetSDW(RngWrapper& rng, ModelParams pars) {
	//TODO: add checks
	pars = updateTemperatureParameters(pars);

	//check parameters: passed all that are necessary
	using namespace boost::assign;
	std::vector<std::string> neededModelPars;
	neededModelPars += "mu", "L", "r", "accRatio";
	for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
		if (pars.specified.count(*p) == 0) {
			throw ParameterMissing(*p);
		}
	}

	if (pars.checkerboard and pars.L % 2 != 0) {
		throw ParameterWrong("Checker board decomposition only supported for even linear lattice sizes");
	}

#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x) 	{  					  						\
								IF_NOT_POSITIVE(x) {  					\
									throw ParameterWrong(#x, pars.x);	\
								}										\
							}
	CHECK_POSITIVE(L);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

	return std::unique_ptr<DetSDW>(new DetSDW(rng, pars));
}

DetSDW::DetSDW(RngWrapper& rng_, const ModelParams& pars) :
		DetModelGC<1,cpx>(pars, 4 * pars.L*pars.L),
		rng(rng_),
		checkerboard(pars.checkerboard),
		L(pars.L), N(L*L), r(pars.r), mu(pars.mu), c(1), u(1), lambda(1), //TODO: make these controllable by parameter
		spaceNeigh(L), timeNeigh(m),
		propK(), propKx(propK[XBAND]), propKy(propK[YBAND]),
		checkerEmKx(N,N), checkerEmKy(N,N), checkerEpKx(N,N), checkerEpKy(N,N),
		g(green[0]), gFwd(greenFwd[0]), gBwd(greenBwd[0]),
		phi0(N, m+1), phi1(N, m+1), phi2(N, m+1), phiCosh(N, m+1), phiSinh(N, m+1),
		phiDelta(InitialPhiDelta),
		targetAccRatio(pars.accRatio), lastAccRatio(0), accRatioRA(AccRatioAdjustmentSamples),
		normPhi(), phiSecond(), phiFourth(), binder(), sdwSusc(),
		kOcc(), kOccX(kOcc[XBAND]), kOccY(kOcc[YBAND]),
		kOccImag(), kOccXimag(kOccImag[XBAND]), kOccYimag(kOccImag[YBAND])
{
	g = CubeCpx(4*N,4*N, m+1);
	if (pars.timedisplaced) {
		gFwd = CubeCpx(4*N,4*N, m+1);
		gBwd = CubeCpx(4*N,4*N, m+1);
	}
	setupRandomPhi();

	if (checkerboard) {
		setupPropK_checkerboard();
		computeBmat[0] = [this](unsigned k2, unsigned k1) {
			return this->computeBmatSDW_checkerboard(k2, k1);
		};
	} else {
		setupPropK_direct();
		computeBmat[0] = [this](unsigned k2, unsigned k1) {
			return this->computeBmatSDW_direct(k2, k1);
		};
	}
	setupUdVStorage();

	using std::cref;
	using namespace boost::assign;
	obsScalar += ScalarObservable(cref(normPhi), "normPhi", "np"),
			ScalarObservable(cref(phiSecond), "phiSecond", "p2"),
			ScalarObservable(cref(phiFourth), "phiFourth", "p4"),
			ScalarObservable(cref(sdwSusc), "sdwSusceptibility", "sdwsusc"),
			ScalarObservable(cref(lastAccRatio), "accRatio", "ar");

	kOccX.zeros(N);
	kOccY.zeros(N);
	obsVector += VectorObservable(cref(kOccX), N, "kOccX", "nkx"),
			VectorObservable(cref(kOccY), N, "kOccY", "nky");
	kOccXimag.zeros(N);
	kOccYimag.zeros(N);
	obsVector += VectorObservable(cref(kOccXimag), N, "kOccXimag", "nkximag"),
			VectorObservable(cref(kOccYimag), N, "kOccYimag", "nkyimag");
}

DetSDW::~DetSDW() {
}

unsigned DetSDW::getSystemN() const {
	return N;
}

MetadataMap DetSDW::prepareModelMetadataMap() const {
	MetadataMap meta;
	meta["model"] = "sdw";
	meta["checkerboard"] = (checkerboard ? "true" : "false");
	meta["timedisplaced"] = (timedisplaced ? "true" : "false");
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
	META_INSERT(targetAccRatio);
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
	timing.start("sdw-measure");
	Phi meanPhi;
	meanPhi[0] = averageWholeSystem(phi0, 0.0);
	meanPhi[1] = averageWholeSystem(phi1, 0.0);
	meanPhi[2] = averageWholeSystem(phi2, 0.0);
	normPhi = arma::norm(meanPhi, 2);
//	phiSecond = averageWholeSystem( [this](unsigned i, unsigned k) {
//		using std::pow;
//		return pow(phi0(i, k), 2) + pow(phi1(i, k), 2) + pow(phi2(i, k), 2);
//	}, 0.0);
//	phiFourth = averageWholeSystem( [this](unsigned i, unsigned k) {
//		using std::pow;
//		return pow(pow(phi0(i, k), 2) + pow(phi1(i, k), 2) + pow(phi2(i, k), 2), 2);
//	}, 0.0);
	phiSecond = std::pow(normPhi, 2);
	phiFourth = std::pow(phiSecond, 2);

	//fermion occupation number
	kOccX.zeros(N);
	kOccY.zeros(N);
	kOccXimag.zeros(N);
	kOccYimag.zeros(N);
	static const num pi = M_PI;
#pragma omp parallel for
	for (unsigned ksitey = 0; ksitey < L; ++ksitey) {			//k-vectors
		num ky = 2 * pi * num(ksitey) / num(L);
		for (unsigned ksitex = 0; ksitex < L; ++ksitex) {
			num kx = 2 * pi * num(ksitex) / num(L);
			unsigned ksite = L*ksitey + ksitex;
			for (unsigned l = 1; l <= m; ++l) {					//timeslices
				for (unsigned jy = 0; jy < L; ++jy) {			//sites j
					for (unsigned jx = 0; jx < L; ++jx) {
						unsigned j = L*jy + jx;
						for (unsigned iy = 0; iy < L; ++iy) {	//sites i
							for (unsigned ix = 0; ix < L; ++ix) {
								unsigned i = L*iy + ix;
								cpx phase = std::exp(cpx(0, kx * (ix - jx) + ky * (iy - jy)));
								cpx diracDelta = cpx((i == j ? 1.0 : 0.0), 0);
								cpx greenEntryXBandSpinUp   = g.slice(l)(i, j);
								cpx greenEntryXBandSpinDown = g.slice(l)(i + N, j + N);
								cpx greenEntryYBandSpinUp   = g.slice(l)(i + 2*N, j + 2*N);
								cpx greenEntryYBandSpinDown = g.slice(l)(i + 3*N, j + 3*N);
								kOccX[ksite] += std::real(phase * ( diracDelta - greenEntryXBandSpinUp
																  + diracDelta - greenEntryXBandSpinDown));
								kOccY[ksite] += std::real(phase * ( diracDelta - greenEntryYBandSpinUp
																  + diracDelta - greenEntryYBandSpinDown));
								//imaginary parts should add up to zero..., but for now check:
								kOccXimag[ksite] += std::imag(phase * ( diracDelta - greenEntryXBandSpinUp
																  	  + diracDelta - greenEntryXBandSpinDown));
								kOccYimag[ksite] += std::imag(phase * ( diracDelta - greenEntryYBandSpinUp
																  	  + diracDelta - greenEntryYBandSpinDown));
							}
						}
					}
				}
			}
		}
	}
	using std::ref;
	for (VecNum& kocc : {ref(kOccX), ref(kOccY), ref(kOccXimag), ref(kOccYimag)}) {
		kocc /= num(m) * num(N);
	}

	//sdw-susceptibility
	sdwSusc = dtau * sumWholeSystem( [this](unsigned site, unsigned timeslice) {
											return phi0(site, timeslice) * phi0(0, m)
												 + phi1(site, timeslice) * phi1(0, m)
												 + phi2(site, timeslice) * phi2(0, m);
										},
									0.0);

	timing.stop("sdw-measure");
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
			phiCosh(site, k) = std::cosh(dtau * phiNorm);
			phiSinh(site, k) = std::sinh(dtau * phiNorm) / phiNorm;
		} );
	} );
}

void DetSDW::setupPropK_direct() {
	std::array<std::array<num,z>, 2> t;
	t[XBAND][XPLUS] = t[XBAND][XMINUS] = -1.0;
	t[XBAND][YPLUS] = t[XBAND][YMINUS] = -0.5;
	t[YBAND][XPLUS] = t[YBAND][XMINUS] =  0.5;
	t[YBAND][YPLUS] = t[YBAND][YMINUS] =  1.0;

	for_each_band( [this, &t](unsigned band) {
		MatNum k = -mu * arma::eye(N,N);
		for_each_site( [this, band, &k, &t](unsigned site) {
			for (unsigned dir = 0; dir < z; ++dir) {
				unsigned neigh = spaceNeigh(dir, site);
				k(site, neigh) -= t[band][dir];
			}
		} );
//		std::string name = std::string("k") + (band == XBAND ? "x" : band == YBAND ? "y" : "error");
//		debugSaveMatrix(k, name);
		propK[band] = computePropagator(dtau, k);
//		debugSaveMatrix(propK[band], "prop" + name);
	} );
}


MatCpx DetSDW::computeBmatSDW_direct(unsigned k2, unsigned k1) const {
	timing.start("computeBmatSDW_direct");
	using arma::eye; using arma::zeros; using arma::diagmat;
	if (k2 == k1) {
		return MatCpx(eye(4*N,4*N), zeros(4*N,4*N));
	}
	assert(k2 > k1);
	assert(k2 <= m);

	//compute the matrix e^(-dtau*V_k) * e^(-dtau*K)
	auto singleTimesliceProp = [this, N](unsigned k) {
		timing.start("singleTimesliceProp_direct");
		MatCpx result(4*N, 4*N);

		//submatrix view helper for a 4N*4N matrix
		auto block = [&result, N](unsigned row, unsigned col) {
			return result.submat(row * N, col * N,
					             (row + 1) * N - 1, (col + 1) * N - 1);
		};
		auto& kphi0 = phi0.col(k);
		auto& kphi1 = phi1.col(k);
		auto& kphi2 = phi2.col(k);
//		debugSaveMatrix(kphi0, "kphi0");
//		debugSaveMatrix(kphi1, "kphi1");
//		debugSaveMatrix(kphi2, "kphi2");
		auto& kphiCosh = phiCosh.col(k);
		auto& kphiSinh = phiSinh.col(k);
		//TODO: is this the best way to set the real and imaginary parts of a complex submatrix?
		//TODO: compare to using set_real / set_imag
		block(0, 0) = MatCpx(diagmat(kphiCosh) * propKx,
				             zeros(N,N));
		block(0, 1).zeros();
		block(0, 2) = MatCpx(diagmat(-kphi2 % kphiSinh) * propKy,
				             zeros(N,N));
		block(0, 3) = MatCpx(diagmat(-kphi0 % kphiSinh) * propKy,
				             diagmat(+kphi1 % kphiSinh) * propKy);
		block(1, 0).zeros();
		block(1, 1) = block(0, 0);
		block(1, 2) = MatCpx(diagmat(-kphi0 % kphiSinh) * propKy,
				             diagmat(-kphi1 % kphiSinh) * propKy);
		block(1, 3) = MatCpx(diagmat(+kphi2 % kphiSinh) * propKy,
				             zeros(N,N));
		block(2, 0) = MatCpx(diagmat(-kphi2 % kphiSinh) * propKx,
				             zeros(N,N));
		block(2, 1) = MatCpx(diagmat(-kphi0 % kphiSinh) * propKx,
				             diagmat(+kphi1 % kphiSinh) * propKx);
		block(2, 2) = MatCpx(diagmat(kphiCosh) * propKy,
				             zeros(N,N));
		block(2, 3).zeros();
		block(3, 0) = MatCpx(diagmat(-kphi0 % kphiSinh) * propKx,
				             diagmat(-kphi1 % kphiSinh) * propKx);
		block(3, 1) = MatCpx(diagmat(+kphi2 % kphiSinh) * propKx,
				             zeros(N,N));
		block(3, 2).zeros();
		block(3, 3) = block(2, 2);

//		debugSaveMatrix(arma::real(result), "emdtauVemdtauK_real");
//		debugSaveMatrix(arma::imag(result), "emdtauVemdtauK_imag");
		timing.stop("singleTimesliceProp_direct");
		return result;
	};

	MatCpx result = singleTimesliceProp(k2);

	for (unsigned k = k2 - 1; k > k1; --k) {
		result *= singleTimesliceProp(k);				// equivalent to: result = result * singleTimesliceProp(k);
	}

	timing.stop("computeBmatSDW_direct");

	return result;
}


void DetSDW::setupPropK_checkerboard() {
	//hopping constants:
	std::array<std::array<num,2>, 2> t;
	enum { HOR = 0, VER = 1 };
	t[XBAND][HOR] = -1.0;
	t[XBAND][VER] = -0.5;
	t[YBAND][HOR] = 0.5;
	t[YBAND][VER] = 1.0;

	//hopping bond matrices, x-band
	SpMatNum kx_a_hor = SpMatNum(N, N);
	SpMatNum kx_b_hor = SpMatNum(N, N);
	SpMatNum kx_a_ver = SpMatNum(N, N);
	SpMatNum kx_b_ver = SpMatNum(N, N);

	//hopping bond matrices, y-band
	SpMatNum ky_a_hor = SpMatNum(N, N);
	SpMatNum ky_b_hor = SpMatNum(N, N);
	SpMatNum ky_a_ver = SpMatNum(N, N);
	SpMatNum ky_b_ver = SpMatNum(N, N);

	auto xyToSite = [L](unsigned x, unsigned y) {
		return y * L + x;
	};

	//horizontal bonds
	for (unsigned y = 0; y < L; ++y) {
		for (unsigned x = 0; x < L; x += 2) {
			//sub board a, horizontal
			//=======================
			unsigned siteA = xyToSite(x, y);
			unsigned neighA = spaceNeigh(NeighDir::XPLUS, siteA);
			//x-band
			kx_a_hor(siteA, neighA) = 1;
			kx_a_hor(neighA, siteA) = 1;
			//y-band
			ky_a_hor(siteA, neighA) = 1;
			ky_a_hor(neighA, siteA) = 1;

			//sub board b, horizontal
			//=======================
			unsigned siteB = neighA;
			unsigned neighB = spaceNeigh(NeighDir::XPLUS, siteB);
			//x-band
			kx_b_hor(siteB, neighB) = 1;
			kx_b_hor(neighB, siteB) = 1;
			//y-band
			ky_b_hor(siteB, neighB) = 1;
			ky_b_hor(neighB, siteB) = 1;
		}
	}
	//vertical bonds
	for (unsigned x = 0; x < L; ++x) {
		for (unsigned y = 0; y < L; y += 2) {
			//sub board a, vertical
			//=====================
			unsigned siteA = xyToSite(x, y);
			unsigned neighA = spaceNeigh(NeighDir::YPLUS, siteA);
			//x-band
			kx_a_ver(siteA, neighA) = 1;
			kx_a_ver(neighA, siteA) = 1;
			//y-band
			ky_a_ver(siteA, neighA) = 1;
			ky_a_ver(neighA, siteA) = 1;

			//sub board b, vertical
			//=====================
			unsigned siteB = neighA;
			unsigned neighB = spaceNeigh(NeighDir::YPLUS, siteB);
			//x-band
			kx_b_ver(siteB, neighB) = 1;
			kx_b_ver(neighB, siteB) = 1;
			//y-band
			ky_b_ver(siteB, neighB) = 1;
			ky_b_ver(neighB, siteB) = 1;
		}
	}

	using std::sinh; using std::cosh; using std::exp; using arma::speye;

	SpMatNum em_kx_a_hor = -sinh(dtau * t[XBAND][HOR]) * kx_a_hor + cosh(dtau * t[XBAND][HOR]) * speye(N,N);
	SpMatNum em_kx_b_hor = -sinh(dtau * t[XBAND][HOR]) * kx_b_hor + cosh(dtau * t[XBAND][HOR]) * speye(N,N);
	SpMatNum em_kx_a_ver = -sinh(dtau * t[XBAND][VER]) * kx_a_ver + cosh(dtau * t[XBAND][VER]) * speye(N,N);
	SpMatNum em_kx_b_ver = -sinh(dtau * t[XBAND][VER]) * kx_b_ver + cosh(dtau * t[XBAND][VER]) * speye(N,N);

	checkerEmKx = exp(dtau * mu) * em_kx_a_hor * em_kx_a_ver * em_kx_b_hor * em_kx_b_ver;
	std::cout << "checkerEmKx: " << num(checkerEmKx.n_nonzero) / num(N*N) << std::endl;

	SpMatNum em_ky_a_hor = -sinh(dtau * t[YBAND][HOR]) * ky_a_hor + cosh(dtau * t[YBAND][HOR]) * speye(N,N);
	SpMatNum em_ky_b_hor = -sinh(dtau * t[YBAND][HOR]) * ky_b_hor + cosh(dtau * t[YBAND][HOR]) * speye(N,N);
	SpMatNum em_ky_a_ver = -sinh(dtau * t[YBAND][VER]) * ky_a_ver + cosh(dtau * t[YBAND][VER]) * speye(N,N);
	SpMatNum em_ky_b_ver = -sinh(dtau * t[YBAND][VER]) * ky_b_ver + cosh(dtau * t[YBAND][VER]) * speye(N,N);

	checkerEmKy = exp(dtau * mu) * em_ky_a_hor * em_ky_a_ver * em_ky_b_hor * em_ky_b_ver;
	std::cout << "checkerEmKy: " << num(checkerEmKy.n_nonzero) / num(N*N) << std::endl;

	SpMatNum ep_kx_a_hor = sinh(dtau * t[XBAND][HOR]) * kx_a_hor + cosh(dtau * t[XBAND][HOR]) * speye(N,N);
	SpMatNum ep_kx_b_hor = sinh(dtau * t[XBAND][HOR]) * kx_b_hor + cosh(dtau * t[XBAND][HOR]) * speye(N,N);
	SpMatNum ep_kx_a_ver = sinh(dtau * t[XBAND][VER]) * kx_a_ver + cosh(dtau * t[XBAND][VER]) * speye(N,N);
	SpMatNum ep_kx_b_ver = sinh(dtau * t[XBAND][VER]) * kx_b_ver + cosh(dtau * t[XBAND][VER]) * speye(N,N);

	checkerEpKx = exp(-dtau * mu) * ep_kx_a_hor * ep_kx_a_ver * ep_kx_b_hor * ep_kx_b_ver;
	std::cout << "checkerEpkx: " << num(checkerEpKx.n_nonzero) / num(N*N) << std::endl;

	SpMatNum ep_ky_a_hor = sinh(dtau * t[YBAND][HOR]) * ky_a_hor + cosh(dtau * t[YBAND][HOR]) * speye(N,N);
	SpMatNum ep_ky_b_hor = sinh(dtau * t[YBAND][HOR]) * ky_b_hor + cosh(dtau * t[YBAND][HOR]) * speye(N,N);
	SpMatNum ep_ky_a_ver = sinh(dtau * t[YBAND][VER]) * ky_a_ver + cosh(dtau * t[YBAND][VER]) * speye(N,N);
	SpMatNum ep_ky_b_ver = sinh(dtau * t[YBAND][VER]) * ky_b_ver + cosh(dtau * t[YBAND][VER]) * speye(N,N);

	checkerEpKy = exp(-dtau * mu) * ep_ky_a_hor * ep_ky_a_ver * ep_ky_b_hor * ep_ky_b_ver;
	std::cout << "checkerEpKy: " << num(checkerEpKy.n_nonzero) / num(N*N) << std::endl;
}


MatCpx DetSDW::computeBmatSDW_checkerboard(unsigned k2, unsigned k1) const {
	using arma::diagmat; using arma::eye; using arma::zeros;
	timing.start("computeBmatSDW_checkerboard");

	if (k2 == k1) {
		return MatCpx(eye(4*N,4*N), zeros(4*N,4*N));
	}
	assert(k2 > k1);
	assert(k2 <= m);

	//compute the matrix e^(-dtau*V_k) * e^(-dtau*K) * e^(dtau * mu)
	auto singleTimesliceProp = [this, N](unsigned k) {
//		timing.start("singleTimesliceProp_checkerboard");
		SpMatCpx result(4*N, 4*N);

		//submatrix view helper for a 4N*4N matrix //TODO: specify return value as a reference?
		auto block = [&result, N](unsigned row, unsigned col) {
			return result.submat(row * N, col * N,
					             (row + 1) * N - 1, (col + 1) * N - 1);
		};
		//Helper to compute the product of a diagonal matrix with a sparse matrix (NxN).
		//This is written with sparse matrixes without any empty columns in mind.
		auto diagTimesSpMat = [N](const VecNum& diag, const SpMatNum& spmat) {
			SpMatNum result(N,N);
			for (unsigned c = 0; c < N; ++c) {
				result.col(c) = diag % spmat.col(c);
			}
			return result;
		};
		//TODO: neaten this up... might be more efficient
//		auto blockAssignRealDiagTimesSpMat = [N, result&, block](unsigned blockRow, unsigned blockCol,
//				const VecNum& diag, const SpMatNum& spmat) -> void {
//			for (unsigned c = 0; c < N; ++c) {
//				block(blockRow, blockCol).col(c)) = diag % spmat.col(c), SpMatNum(N,N));
//			}
//		};


		const auto& kphi0 = phi0.col(k);
		const auto& kphi1 = phi1.col(k);
		const auto& kphi2 = phi2.col(k);
		const auto& kphiCosh = phiCosh.col(k);
		const auto& kphiSinh = phiSinh.col(k);

		block(0, 0) = SpMatCpx(diagTimesSpMat(kphiCosh, checkerEmKx),
							   SpMatNum(N,N));
		//block(0, 1) zero
		block(0, 2) = SpMatCpx(diagTimesSpMat(-kphi2 % kphiSinh, checkerEmKy),
				         	   SpMatNum(N,N));
		block(0, 3) = SpMatCpx(diagTimesSpMat(-kphi0 % kphiSinh, checkerEmKy),
				               diagTimesSpMat(+kphi1 % kphiSinh, checkerEmKy));
		//block(1, 0) zero
		block(1, 1) = block(0, 0);
		block(1, 2) = SpMatCpx(diagTimesSpMat(-kphi0 % kphiSinh, checkerEmKy),
				               diagTimesSpMat(-kphi1 % kphiSinh, checkerEmKy));
		block(2, 0) = SpMatCpx(diagTimesSpMat(-kphi2 % kphiSinh, checkerEmKx),
				               SpMatNum(N,N));
		block(2, 1) = SpMatCpx(diagTimesSpMat(-kphi0 % kphiSinh, checkerEmKx),
				               diagTimesSpMat(+kphi1 % kphiSinh, checkerEmKx));
		block(2, 2) = SpMatCpx(diagTimesSpMat(kphiCosh, checkerEmKy),
							   SpMatNum(N,N));
		//block(2, 3) zero
		block(3, 0) = SpMatCpx(diagTimesSpMat(-kphi0 % kphiSinh, checkerEmKx),
							   diagTimesSpMat(-kphi1 % kphiSinh, checkerEmKx));
		block(3, 1) = SpMatCpx(diagTimesSpMat(+kphi2 % kphiSinh, checkerEmKx),
				               SpMatNum(N,N));
		//block(3, 2) zero
		block(3, 3) = block(2, 2);

//		timing.stop("singleTimesliceProp_checkerboard");
		return result;
	};

	//TODO: try rvalue references, so I don't have to rely on named return value optimization being done

	SpMatCpx result = singleTimesliceProp(k2);

	for (unsigned k = k2 - 1; k > k1; --k) {
		result *= singleTimesliceProp(k);
	}

	timing.stop("computeBmatSDW_checkerboard");

	return MatCpx(result);
}






void DetSDW::updateInSlice(unsigned timeslice) {
	timing.start("sdw-updateInSlice");

	lastAccRatio = 0;
	for_each_site( [this, timeslice](unsigned site) {
		Phi newphi = proposeNewField(site, timeslice);

//		VecNum oldphi0 = phi0.col(timeslice);
//		VecNum oldphi1 = phi1.col(timeslice);
//		VecNum oldphi2 = phi2.col(timeslice);
//		debugSaveMatrix(oldphi0, "old_phi0");
//		debugSaveMatrix(oldphi1, "old_phi1");
//		debugSaveMatrix(oldphi2, "old_phi2");

//		VecNum newphi0 = phi0.col(timeslice);
//		VecNum newphi1 = phi1.col(timeslice);
//		VecNum newphi2 = phi2.col(timeslice);
//		newphi0[site] = newphi[0];
//		newphi1[site] = newphi[1];
//		newphi2[site] = newphi[2];
//		debugSaveMatrix(newphi0, "new_phi0");
//		debugSaveMatrix(newphi1, "new_phi1");
//		debugSaveMatrix(newphi2, "new_phi2");

		num dsphi = deltaSPhi(site, timeslice, newphi);
		num propSPhi = std::exp(-dsphi);
//		std::cout << propSPhi << std::endl;

		//delta = e^(-dtau*V_new)*e^(+dtau*V_old) - 1

		//compute non-zero elements of delta
		//evMatrix(): yield a 4x4 matrix containing the entries for the
		//current lattice site and time slice of e^(sign*dtau*V) with
		//given values of the field phi at that space-time location [and of
		//cosh(dtau*|phi|) and sinh(dtau*|phi|) / |phi|]
		auto evMatrix = [](int sign, num kphi0, num kphi1,
						   num kphi2, num kphiCosh, num kphiSinh) -> MatCpx::fixed<4,4> {
			MatNum::fixed<4,4> ev_real;
			ev_real.diag().fill(kphiCosh);
			ev_real(0,1) = ev_real(1,0) = ev_real(2,3) = ev_real(3,2) = 0;
			ev_real(2,0) = ev_real(0,2) =  sign * kphi2 * kphiSinh;
			ev_real(2,1) = ev_real(0,3) =  sign * kphi0 * kphiSinh;
			ev_real(3,0) = ev_real(1,2) =  sign * kphi0 * kphiSinh;
			ev_real(3,1) = ev_real(1,3) = -sign * kphi2 * kphiSinh;

			MatCpx::fixed<4,4> ev;
			ev.set_real(ev_real);
			ev(0,3).imag(-sign * kphi1 * kphiSinh);
			ev(1,2).imag( sign * kphi1 * kphiSinh);
			ev(2,1).imag(-sign * kphi1 * kphiSinh);
			ev(3,0).imag( sign * kphi1 * kphiSinh);

			return ev;
		};
		MatCpx::fixed<4,4> evOld = evMatrix(
				+1,
				phi0(site, timeslice), phi1(site, timeslice), phi2(site, timeslice),
				phiCosh(site, timeslice), phiSinh(site, timeslice)
				);
		num normnewphi = arma::norm(newphi,2);
		num coshnewphi = std::cosh(dtau * normnewphi);
		num sinhnewphi = std::sinh(dtau * normnewphi) / normnewphi;
		MatCpx::fixed<4,4> emvNew = evMatrix(
				-1,
				newphi[0], newphi[1], newphi[2],
				coshnewphi, sinhnewphi
				);
		MatCpx::fixed<4,4> deltanonzero = emvNew * evOld;
		deltanonzero.diag() -= cpx(1.0, 0);

		//****
		//Compute the determinant and inverse of I + Delta*(I - G)
		//based on Sherman-Morrison formula / Matrix-Determinant lemma
		//****

		//Delta*(I - G) is a sparse matrix containing just 4 rows:
		//site, site+N, site+2N, site+3N
		//Compute the values of these rows [O(N)]:
		std::array<VecCpx, 4> rows = {{VecCpx(4*N), VecCpx(4*N), VecCpx(4*N), VecCpx(4*N)}};
#pragma omp parallel for
		for (unsigned r = 0; r < 4; ++r) {
			for (unsigned col = 0; col < 4*N; ++col) {
				rows[r][col] = -deltanonzero(r,0) * g.slice(timeslice).col(col)[site];
			}
			rows[r][site] += deltanonzero(r,0);
			for (unsigned dc = 1; dc < 4; ++dc) {
				for (unsigned col = 0; col < 4*N; ++col) {
					rows[r][col] += -deltanonzero(r,dc) * g.slice(timeslice).col(col)[site + dc*N];
				}
				rows[r][site + dc*N] += deltanonzero(r,dc);
			}
		}

		// [I + Delta*(I - G)]^(-1) again is a sparse matrix
		// with four rows site, site+N, site+2N, site+3N
		// compute them iteratively, together with the determinant of
		// I + Delta*(I - G)
		// Apart from these rows, the remaining diagonal entries of
		// [I + Delta*(I - G)]^(-1) are 1
		//
		// before this loop rows[] holds the entries of Delta*(I - G),
		// after the loop rows[] holds the corresponding rows of [I + Delta*(I - G)]^(-1)
		cpx det = 1;
		for (unsigned l = 0; l < 4; ++l) {
			VecCpx row = rows[l];
			for (int k = l-1; k >= 0; --k) {
				row[site + k*N] = 0;
			}
			for (int k = l-1; k >= 0; --k) {
				row += rows[l][site + k*N] * rows[k];
			}
			cpx divisor = cpx(1.0, 0) + row[site + l*N];
			rows[l] = (-1.0/divisor) * row;
			rows[l][site + l*N] += 1;
			for (int k = l - 1; k >= 0; --k) {
				rows[k] -= (rows[k][site + l*N] / divisor) * row;
			}
			det *= divisor;
		}

		//****DEBUG
//		std::array<VecCpx, 4> invRows = {{rows[0], rows[1], rows[2], rows[3]}};
		//****END-DEBUG

		//****
		//DEBUG: This slow code was working before.
		//Compare results with the better-performing sherman-morrison code
//		SpMatCpx delta(4*N, 4*N);
//		arma::uvec::fixed<4> idx = {site, site + N, site + 2*N, site + 3*N};
//		//Armadilo lacks non-contiguous submatrix views for sparse matrices
//		unsigned i = 0;
//		for (auto col: idx) {
//			unsigned j = 0;
//			for (auto row: idx) {
//				delta(row, col) = deltanonzero(j, i);
//				++j;
//			}
//			++i;
//		}


//		MatCpx deltaDense(4*N,4*N);
//		deltaDense = delta;
//		debugSaveMatrix(MatNum(arma::real(deltaDense)), "delta_real");
//		debugSaveMatrix(MatNum(arma::imag(deltaDense)), "delta_imag");

//		//inefficient!
//		static MatCpx eyeCpx = MatCpx(arma::eye(4*N, 4*N), arma::zeros(4*N, 4*N));
//		MatCpx target = eyeCpx + delta * (eyeCpx - g.slice(timeslice));

//		debugSaveMatrix(MatNum(arma::real(target)), "target_real");
//		debugSaveMatrix(MatNum(arma::imag(target)), "target_imag");

//		cpx weightRatio = arma::det(target);
//		std::cout << weightRatio << std::endl;

//		std::cout << weightRatio << " vs. " << det << std::endl;
		//END DEBUG
		//****

		//****
		// DEBUG-CHECK Delta*(I - G) vs. rows
//		MatCpx check = delta * (eyeCpx - g.slice(timeslice));
//		for (unsigned r = 0; r < 4; ++r) {
//			check.row(site + r*N) -= rows[r].st();
//		}
//		std::cout << "Row check: " << arma::max(arma::max(arma::abs(check))) << std::endl;
//		debugSaveMatrix(MatNum(arma::real(check)), "check_real");
//		debugSaveMatrix(MatNum(arma::imag(check)), "check_imag");
//		exit(0);
		// END-DEBUG-CHECK
		//****

		//****
		//DEBUG-CHECK [I+Delta*(I - G)]^(-1) vs. invRows
//		MatCpx inv = arma::inv(target);
//		inv -= eyeCpx;
//		for (unsigned r = 0; r < 4; ++r) {
//			inv.row(site + r*N) -= invRows[r].st();
//			inv.row(site + r*N)[site + r*N] += 1;
//		}
//		std::cout << "inv-div: " << arma::max(arma::max(inv)) << std::endl;
		//END-DEBUG-CHECK
		//****

		num propSFermion = det.real();

		num prop = propSPhi * propSFermion;

		if (prop > 1.0 or rng.rand01() < prop) {
			//count accepted update
			lastAccRatio += 1.0;

//			num phisBefore = phiAction();
			phi0(site, timeslice) = newphi[0];
			phi1(site, timeslice) = newphi[1];
			phi2(site, timeslice) = newphi[2];
			phiCosh(site, timeslice) = coshnewphi;
			phiSinh(site, timeslice) = sinhnewphi;
//			num phisAfter = phiAction();
//			std::cout << std::scientific << dsphi << " vs. " << phisAfter << " - " << phisBefore << " = " <<
//					(phisAfter - phisBefore) << std::endl;

//			debugSaveMatrix(MatNum(arma::real(g.slice(timeslice))), "gslice_old_real");
//			debugSaveMatrix(MatNum(arma::imag(g.slice(timeslice))), "gslice_old_imag");
//			g.slice(timeslice) *= arma::inv(target);
//			debugSaveMatrix(MatNum(arma::real(g.slice(timeslice))), "gslice_new_real");
//			debugSaveMatrix(MatNum(arma::imag(g.slice(timeslice))), "gslice_new_imag");

			//****
			//DEBUG
//			MatCpx gPrimeRef = g.slice(timeslice) * arma::inv(target);
			//END DEBUG
			//****

			//compensate for already included diagonal entries of I in invRows
			rows[0][site] -= 1;
			rows[1][site + N] -= 1;
			rows[2][site + 2*N] -= 1;
			rows[3][site + 3*N] -= 1;
			//compute G' = G * [I + Delta*(I - G)]^(-1) = G * [I + invRows]
			MatCpx gTimesInvRows(4*N, 4*N);
			const auto& G = g.slice(timeslice);
#pragma omp parallel for
			for (unsigned col = 0; col < 4*N; ++col) {
				for (unsigned row = 0; row < 4*N; ++row) {
					gTimesInvRows(row, col) = G(row, site) * rows[0][col]
					                        + G(row, site + N) * rows[1][col]
					                        + G(row, site + 2*N) * rows[2][col]
					                        + G(row, site + 3*N) * rows[3][col]
					                        ;
				}
			}
			g.slice(timeslice) += gTimesInvRows;

			//****
			//DEBUG
//			std::cout << arma::max(arma::max(
//					(arma::abs(gPrimeRef - g.slice(timeslice))))) << std::endl;
			//END DEBUG
			//****
		}
	});
	lastAccRatio /= num(N);

	timing.stop("sdw-updateInSlice");
}

void DetSDW::updateInSliceThermalization(unsigned timeslice) {
	updateInSlice(timeslice);

	accRatioRA.addValue(lastAccRatio);
	if (accRatioRA.getSamplesAdded() % AccRatioAdjustmentSamples == 0) {
		num avgAccRatio = accRatioRA.get();
		if (avgAccRatio < targetAccRatio) {
			phiDelta *= phiDeltaShrinkFactor;
		} else if (avgAccRatio > targetAccRatio) {
			phiDelta *= phiDeltaGrowFactor;
		}
	}
}





DetSDW::Phi DetSDW::proposeNewField(unsigned site, unsigned timeslice) {
	(void) site; (void) timeslice;
	//TODO: make this smarter!

	Phi phi = { phi0(site, timeslice),
			    phi1(site, timeslice),
			    phi2(site, timeslice)
	          };

//	static int comp = 0;
//	num r = rng.randRange(-PhiDelta, +PhiDelta);
//	phi[comp] += r;
//	comp = (comp + 1) % 3;

	for (auto& comp: phi) {
		num r = rng.randRange(-phiDelta, +phiDelta);
		comp += r;
	}

	return phi;
}

num DetSDW::deltaSPhi(unsigned site, unsigned timeslice, const Phi newphi) {
	//switched to asymmetric numerical derivative
	using arma::dot;

	const Phi oldphi = {phi0(site, timeslice),
					    phi1(site, timeslice),
					    phi2(site, timeslice)};
	Phi phiDiff = newphi - oldphi;

	num oldphiSq = dot(oldphi, oldphi);
	num newphiSq = dot(newphi, newphi);
	num phiSqDiff = newphiSq - oldphiSq;

	num oldphiPow4 = oldphiSq * oldphiSq;
	num newphiPow4 = newphiSq * newphiSq;
	num phiPow4Diff = newphiPow4 - oldphiPow4;

	unsigned kEarlier = timeNeigh(ChainDir::MINUS, timeslice);
	Phi phiEarlier = Phi{phi0(site, kEarlier),
					     phi1(site, kEarlier),
					     phi2(site, kEarlier)};
	unsigned kLater = timeNeigh(ChainDir::PLUS, timeslice);
	Phi phiLater = Phi{phi0(site, kLater),
					   phi1(site, kLater),
					   phi2(site, kLater)};
	Phi phiTimeNeigh = phiLater + phiEarlier;

	Phi phiSpaceNeigh = std::accumulate(
			spaceNeigh.beginNeighbors(site),
			spaceNeigh.endNeighbors(site),
			Phi{0,0,0},
			[this, timeslice] (Phi accum, unsigned neighSite) {
				return accum + Phi{phi0(neighSite, timeslice),
					phi1(neighSite, timeslice),
					phi2(neighSite, timeslice)};
			}
	);


	num delta1 = (1.0 / (c * c * dtau)) * (phiSqDiff - dot(phiTimeNeigh, phiDiff));

	num delta2 = 0.5 * dtau * (z * phiSqDiff - 2.0 * dot(phiSpaceNeigh, phiDiff));

	num delta3 = dtau * (0.5 * r * phiSqDiff + 0.25 * u * phiPow4Diff);

	return delta1 + delta2 + delta3;
}


num DetSDW::phiAction() {
	//switched to asymmetric numerical derivative
	arma::field<Phi> phi(N, m+1);
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			phi(site, timeslice)[0] = phi0(site, timeslice);
			phi(site, timeslice)[1] = phi1(site, timeslice);
			phi(site, timeslice)[2] = phi2(site, timeslice);
		}
	}
	num action = 0;
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			Phi timeDerivative =
					(phi(site, timeslice) - phi(site, timeNeigh(ChainDir::MINUS, timeslice)))
					/ dtau;
			action += (dtau / (2.0 * c * c)) * arma::dot(timeDerivative, timeDerivative);

			//count only neighbors in PLUS-directions: no global overcounting of bonds
			Phi xneighDiff = phi(site, timeslice) -
					phi(spaceNeigh(XPLUS, site), timeslice);
			action += 0.5 * dtau * arma::dot(xneighDiff, xneighDiff);
			Phi yneighDiff = phi(site, timeslice) -
					phi(spaceNeigh(YPLUS, site), timeslice);
			action += 0.5 * dtau * arma::dot(yneighDiff, yneighDiff);

			num phisq = arma::dot(phi(site, timeslice), phi(site, timeslice));
			action += 0.5 * dtau * r * phisq;

			action += 0.25 * dtau * u * std::pow(phisq, 2);
		}
	}
	return action;
}


void DetSDW::thermalizationOver() {
	std::cout << "After thermalization: phiDelta = " << phiDelta << std::endl;
}
