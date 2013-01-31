/*
 * dethubbard.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: gerlach
 */

#include <vector>
#include <cmath>
#include <cassert>
#pragma GCC diagnostic ignored "-Weffc++"
#include <boost/assign/std/vector.hpp>    // 'operator+=()' for vectors
#pragma GCC diagnostic warning "-Weffc++"
#include "tools.h"
#include "exceptions.h"
#include "rngwrapper.h"
#include "dethubbard.h"
#include "parameters.h"

//using std::acosh;    //Intel compiler chokes with std::acosh


void debugSaveMatrix(const MatNum& matrix, const std::string& basename) {
	matrix.save(basename + ".csv", arma::csv_ascii);
	matrix.save(basename + ".txt", arma::arma_ascii);
}

void debugSaveMatrix(const MatInt& matrix, const std::string& basename) {
	matrix.save(basename + ".csv", arma::csv_ascii);
	matrix.save(basename + ".txt", arma::arma_ascii);
}


std::unique_ptr<DetHubbard> createDetHubbard(RngWrapper& rng, ModelParams pars) {
	//check parameters: passed all that are necessary
	using namespace boost::assign;
	std::vector<std::string> neededModelPars;
	neededModelPars += "t", "U", "mu", "L", "d", "beta", "s";
	for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
		if (pars.specified.count(*p) == 0) {
			throw ParameterMissing(*p);
		}
	}

	//check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x) 	{  					  						\
								IF_NOT_POSITIVE(x) {  					\
									throw ParameterWrong(#x, pars.x);	\
								}										\
							}
	CHECK_POSITIVE(L);
	CHECK_POSITIVE(d);
	CHECK_POSITIVE(beta);
	CHECK_POSITIVE(m);
	CHECK_POSITIVE(s);
	CHECK_POSITIVE(dtau);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

	//Special handling to allow passing either 'm' or 'dtau', but not both.
	//Also check that 's' is set correctly.
	if (pars.specified.count("dtau") != 0) {
		if (pars.specified.count("m")) {
			throw ParameterWrong("Only specify one of the parameters m and dtau");
		}
		if (pars.s * pars.dtau > pars.beta) {
			throw ParameterWrong("Parameters are incompatible: s * dtau > beta !");
		}
		unsigned n = unsigned(std::ceil(pars.beta / (pars.s * pars.dtau)));
		pars.m = pars.s * n;
	} else if (pars.specified.count("m") == 0) {
		throw ParameterMissing("m");
	}
	if (pars.m % pars.s != 0) {
		throw ParameterWrong("Parameters m=" + numToString(pars.m) + " and s=" + numToString(pars.s)
				+ " do not agree.");
	}

	return std::unique_ptr<DetHubbard>(new DetHubbard(rng, pars));
}



DetHubbard::DetHubbard(RngWrapper& rng_, const ModelParams& pars) :
		rng(rng_),
		t(pars.t), U(pars.U), mu(pars.mu), L(pars.L), d(pars.d),
		z(2*d), //coordination number: 2*d
		N(static_cast<unsigned>(uint_pow(L,d))),
		beta(pars.beta), m(pars.m), s(pars.s), n(m / s), dtau(beta/m),
		alpha(acosh(std::exp(dtau * U * 0.5))),
		nearestNeigbors(z, N),
		proptmat(N,N),
		auxfield(N, m+1),                //m+1 columns of N rows
		gUp(N,N,m+1), gDn(N,N,m+1),      //m+1 slices of N columns x N rows
		gFwdUp(N,N,m+1), gFwdDn(N,N,m+1),
		gBwdUp(N,N,m+1), gBwdDn(N,N,m+1),
		eye_UdV(N),
		UdVStorageUp(), UdVStorageDn(), lastSweepDir(SweepDirection::Up),
		occUp(), occDn(), occTotal(), eKinetic(), ePotential(), eTotal(),
		occDouble(), localMoment(), suscq0(), zcorr(),
		obsNames(), obsShorts(), obsValRefs(), obsCount(0),
		vecObsNames(), vecObsShorts(), vecObsValRefs(), vecObsCount(0)
{
	createNeighborTable();
	setupRandomAuxfield();
	setupPropTmat();
	setupUdVStorage();
	lastSweepDir = SweepDirection::Up;		//first sweep will be downwards

	using namespace boost::assign;         // bring operator+=() into scope
	using std::cref;
	obsNames += "occupationUp", "occupationDown", "totalOccupation",
			"doubleOccupation", "localMoment",
			"kineticEnergy", "potentialEnergy", "totalEnergy", "susceptibilityQ0";
	obsShorts += "nUp", "nDown", "n", "n2", "m^2", "e_t", "e_U", "e", "chi_q0";
	obsValRefs += cref(occUp), cref(occDn), cref(occTotal), cref(occDouble), cref(localMoment),
			cref(eKinetic), cref(ePotential), cref(eTotal), cref(suscq0);
	assert(obsNames.size() == obsShorts.size());
	assert(obsNames.size() == obsValRefs.size());
	obsCount = obsNames.size();

	vecObsNames += "spinzCorrelationFunction";
	vecObsShorts += "zcorr";
	vecObsValRefs += cref(zcorr);
	assert(vecObsNames.size() == vecObsShorts.size());
	assert(vecObsNames.size() == vecObsValRefs.size());
	vecObsCount = vecObsNames.size();

	zcorr.zeros(N);
}

DetHubbard::~DetHubbard() {
}


MetadataMap DetHubbard::prepareModelMetadataMap() {
	MetadataMap meta;
	meta["model"] = "hubbard";
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
	META_INSERT(t);
	META_INSERT(U);
	META_INSERT(mu);
	META_INSERT(L);
	META_INSERT(d);
	META_INSERT(N);
	META_INSERT(beta);
	META_INSERT(m);
	META_INSERT(dtau);
	META_INSERT(s);
	META_INSERT(alpha);
#undef META_INSERT
	return meta;
}

void DetHubbard::updateInSlice(unsigned timeslice) {
	// picking sites linearly: system seemed to alternate between two configurations
	// sweep after sweep
//	for (unsigned site = 0; site < N; ++site) {
//	std::cout << timeslice << " ";  //DEBUG
	for (unsigned count = 0; count < N; ++count) {
		unsigned site = rng.randInt(0, N-1);

		num ratio =  weightRatioSingleFlip(site, timeslice);

//		//DEBUG: comparison of weight ratio calculation
//		MatInt newAuxfield = auxfield;
//		newAuxfield(site, timeslice) *= -1;
//		num refRatio = weightRatioGenericNaive(auxfield, newAuxfield);
//		std::cout << ratio << " vs. " << refRatio << std::endl;

		assert(ratio > 0.0);

		//Metropolis
		if (ratio > 1.0 or rng.rand01() < ratio) {
			//			if (refRatio > 1 or rng.rand01() < refRatio) {
			//				std::cout << "acc" << '\n';
			//				auxfield = newAuxfield;

			updateGreenFunctionWithFlip(site, timeslice);
			auxfield(site, timeslice) *= -1;
		}
	}
}


void DetHubbard::sweepSimple() {
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		gUp.slice(timeslice) = computeGreenFunctionNaive(timeslice, Spin::Up);
		gDn.slice(timeslice) = computeGreenFunctionNaive(timeslice, Spin::Down);
		updateInSlice(timeslice);
	}
}

DetHubbard::UdV DetHubbard::svd(const MatNum& mat) {
	UdV result;
	MatNum V_transpose;
	arma::svd(result.U, result.d, V_transpose, mat, "standard");
	result.V = V_transpose.t();			//potentially it may be advisable to not do this generally
	return result;
}

unsigned DetHubbard::getSystemN() const {
	return N;
}

DetHubbard::MatNum4 DetHubbard::greenFromUdV_timedisplaced(
		const UdV& UdV_l, const UdV& UdV_r) const {
	//Ul vs Vl to be compatible with labeling in the notes
	const MatNum& Ul = UdV_l.V;   //!
	const VecNum& dl = UdV_l.d;
	const MatNum& Vl = UdV_l.U;   //!
	const MatNum& Ur = UdV_r.U;
	const VecNum& dr = UdV_r.d;
	const MatNum& Vr = UdV_r.V;

	//submatrix view helpers for 2*N x 2*N matrices
	//#define upleft(m) m.submat(0,0, N-1,N-1)
	//#define upright(m) m.submat(0,N, N-1,2*N-1)
	//#define downleft(m) m.submat(N,0, 2*N-1,N-1)
	//#define downright(m) m.submat(N,N, 2*N-1,2*N-1)
	auto upleft = [N](MatNum& m) {
		return m.submat(0,0, N-1,N-1);
	};
	auto upright = [N](MatNum& m) {
		return m.submat(0,N, N-1,2*N-1);
	};
	auto downleft = [N](MatNum& m) {
		return m.submat(N,0, 2*N-1,N-1);
	};
	auto downright = [N](MatNum& m) {
		return m.submat(N,N, 2*N-1,2*N-1);
	};

	MatNum temp(2*N,2*N);
	upleft(temp)    = arma::inv(Vr * Vl);
	upright(temp)   = arma::diagmat(dl);
	downleft(temp)  = arma::diagmat(-dr);
	downright(temp) = arma::inv(Ul * Ur);
	UdV tempUdV = svd(temp);

	MatNum left(2*N,2*N);
	upleft(left) = arma::inv(Vr);
	upright(left).zeros();
	downleft(left).zeros();
	downright(left) = arma::inv(Ul);

	MatNum right(2*N,2*N);
	upleft(right) = arma::inv(Vl);
	upright(right).zeros();
	downleft(right).zeros();
	downright(right) = arma::inv(Ur);

	MatNum result = (left * arma::inv(tempUdV.V)) * arma::diagmat(1.0 / tempUdV.d)
					* (arma::inv(tempUdV.U) * right);
	return MatNum4(upleft(result), upright(result),
				   downleft(result), downright(result));
}

MatNum DetHubbard::greenFromUdV(const UdV& UdV_l, const UdV& UdV_r) const {
	//variable names changed according to labeling in notes
	const MatNum& V_l = UdV_l.U;   //!
	const VecNum& d_l = UdV_l.d;
	const MatNum& U_l = UdV_l.V;   //!
	const MatNum& U_r = UdV_r.U;
	const VecNum& d_r = UdV_r.d;
	const MatNum& V_r = UdV_r.V;

	using arma::inv; using arma::diagmat; using arma::eye;

	UdV UdV_temp = svd( inv(U_l * U_r) + diagmat(d_r) * (V_r * V_l) * diagmat(d_l) );

	MatNum green = inv(UdV_temp.V * U_l) * diagmat(1.0 / UdV_temp.d) * inv(U_r * UdV_temp.U);

	return green;
}


void DetHubbard::setupUdVStorage() {
	eye_UdV.U = arma::eye(N,N);
	eye_UdV.d = arma::ones(N);
	eye_UdV.V = arma::eye(N,N);

	auto setup = [this](std::vector<UdV>& storage, Spin spinz) {
		storage = std::vector<UdV>(n + 1);

		storage[0] = eye_UdV;
		storage[1] = svd(computeBmatNaive(s, 0, spinz));

		for (unsigned l = 1; l <= n - 1; ++l) {
			const MatNum& U_l = storage[l].U;
			const VecNum& d_l = storage[l].d;
			const MatNum& V_l = storage[l].V;
			MatNum B_lp1 = computeBmatNaive(s*(l + 1), s*l, spinz);
			UdV UdV_temp = svd((B_lp1 * U_l) * arma::diagmat(d_l));
			storage[l+1].U = UdV_temp.U;
			storage[l+1].d = UdV_temp.d;
			storage[l+1].V = UdV_temp.V * V_l;
		}
	};

	setup(UdVStorageUp, Spin::Up);
	setup(UdVStorageDn, Spin::Down);

	lastSweepDir = SweepDirection::Up;
}


void DetHubbard::debugCheckBeforeSweepDown() {
	std::cout << "Before sweep down:\n";
	std::cout << "up: ";
	for (unsigned l = 0; l <= n; ++l) {
		UdV udv = UdVStorageUp[l];
		MatNum diff = computeBmatNaive(l*s, 0, Spin::Up) - udv.U * arma::diagmat(udv.d) * udv.V;
		std::cout << diff.max() << " ";
	}
	std::cout << "\n";
	std::cout << "down: ";
	for (unsigned l = 0; l <= n; ++l) {
		UdV udv = UdVStorageDn[l];
		MatNum diff = computeBmatNaive(l*s, 0, Spin::Down) - udv.U * arma::diagmat(udv.d) * udv.V;
		std::cout << diff.max() << " ";
	}
	std::cout << "\n\n";
}

void DetHubbard::debugCheckBeforeSweepUp() {
	std::cout << "Before sweep up:\n";
	std::cout << "up: ";
	for (unsigned l = 0; l <= n; ++l) {
		UdV udv = UdVStorageUp[l];
		MatNum diff = computeBmatNaive(m, l*s, Spin::Up) - udv.U * arma::diagmat(udv.d) * udv.V;
		std::cout << diff.max() << " ";
	}
	std::cout << "\n";
	std::cout << "down: ";
	for (unsigned l = 0; l <= n; ++l) {
		UdV udv = UdVStorageDn[l];
		MatNum diff = computeBmatNaive(m, l*s, Spin::Down) - udv.U * arma::diagmat(udv.d) * udv.V;
		std::cout << diff.max() << " ";
	}
	std::cout << "\n\n";
}

void DetHubbard::debugCheckGreenFunctions() {
	std::cout << "debugCheckGreenFunctions:\n";
	std::cout << "up: ";
	for (unsigned k = 1; k <= m; ++k) {
		MatNum green = gUp.slice(k);
		MatNum reldiff = (computeGreenFunctionNaive(k, Spin::Up) - green) / green;
		std::cout << reldiff.max() << " ";
	}
	std::cout << "\n";
	std::cout << "down: ";
	for (unsigned k = 1; k <= m; ++k) {
		MatNum green = gDn.slice(k);
		MatNum reldiff = (computeGreenFunctionNaive(k, Spin::Down) - green) / green;
		std::cout << reldiff.max() << " ";
	}
	std::cout << "\n\n";
}



void DetHubbard::sweep() {
	using std::tie; using std::ignore; using std::get;

	//compute the green function in timeslice s*(l-1) from scratch with the help
	//of the B-matrices computed before in the last up-sweep
	auto advanceDownGreen = [this](unsigned l, std::vector<UdV>& storage,
			CubeNum& green, CubeNum& greenFwd,
			CubeNum& greenBwd, Spin spinz) -> void {
		MatNum B_l = computeBmatNaive(s*l, s*(l - 1), spinz);

		//U_l, d_l, V_l correspond to B(beta,l*s*dtau) [set in the last step]
		const MatNum& U_l = storage[l].U;
		const VecNum& d_l = storage[l].d;
		const MatNum& V_l = storage[l].V;

		//UdV_L will correspond to B(beta,(l-1)*s*dtau)
		UdV UdV_L = svd(arma::diagmat(d_l) * (V_l * B_l));
		UdV_L.U = U_l * UdV_L.U;

		//UdV_R corresponds to B((l-1)*s*dtau,0) [set in last sweep]
		const UdV& UdV_R = storage[l - 1];

		unsigned next = s * (l - 1);
		tie(ignore, greenBwd.slice(next), greenFwd.slice(next), green.slice(next)) =
				greenFromUdV_timedisplaced(UdV_L, UdV_R);
		storage[l - 1] = UdV_L;
	};

	//compute the green function at k-1 by wrapping the one at k (accumulates rounding errors)
	auto wrapDownGreen = [this](unsigned k, CubeNum& green, Spin spinz) -> void {
		MatNum B_k = computeBmatNaive(k, k - 1, spinz);
		green.slice(k - 1) = arma::inv(B_k) * green.slice(k) * B_k;
	};

	//update the green function in timeslice s*(l+1) from scratch with the help
	//of B-matrices computed before
	auto advanceUpGreen = [this](unsigned l, const std::vector<UdV>& storage,
			CubeNum& green, CubeNum& greenFwd,
			CubeNum& greenBwd, Spin spinz) -> void {
		MatNum B_lp1 = computeBmatNaive(s*(l + 1), s*l, spinz);

		//The following is B(beta, (l+1)*s*dtau), valid from the last sweep
		const UdV& UdV_lp1 = storage[l + 1];

		//from the last step the following are B(l*s*dtau, 0):
		const MatNum& U_l = storage[l].U;
		const VecNum& d_l = storage[l].d;
		const MatNum& V_l = storage[l].V;

		//UdV_temp will be the new B((l+1)*s*dtau, 0):
		UdV UdV_temp = svd(((B_lp1 * U_l) * arma::diagmat(d_l)));
		UdV_temp.V *= V_l;

		unsigned next = s * (l + 1);
		tie(ignore, greenBwd.slice(next), greenFwd.slice(next), green.slice(next)) =
				greenFromUdV_timedisplaced(UdV_lp1, UdV_temp);

		//storage[l + 1] = UdV_temp;    //storage would be wrong after updateInSlice!
	};

	//Given B(l*s*dtau, 0) from the last step in the storage, compute
	//B((l+1)*s*dtau, 0) and put it into storage
	auto advanceUpUpdateStorage = [this](unsigned l, std::vector<UdV>& storage,
			Spin spinz) -> void {
		MatNum B_lp1 = computeBmatNaive(s*(l + 1), s*l, spinz);
		//from the last step the following are B(l*s*dtau, 0):
		const MatNum& U_l = storage[l].U;
		const VecNum& d_l = storage[l].d;
		const MatNum& V_l = storage[l].V;
		//the new B((l+1)*s*dtau, 0):
		storage[l+1] = svd(((B_lp1 * U_l) * arma::diagmat(d_l)));
		storage[l+1].V *= V_l;
	};

	//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors)
	auto wrapUpGreen = [this](unsigned k, CubeNum& green, Spin spinz) -> void {
		MatNum B_kp1 = computeBmatNaive(k + 1, k, spinz);
		green.slice(k + 1) = B_kp1 * green.slice(k) * arma::inv(B_kp1);
	};

	if (lastSweepDir == SweepDirection::Up) {
//		debugCheckBeforeSweepDown();
//		debugCheckGreenFunctions();
		//to compute green function for timeslice tau=beta:
		//we need VlDlUl = B(beta, beta) = I and UrDrVr = B(beta, 0).
		//The latter is given in storage slice m from the last sweep.
		tie(ignore, gBwdUp.slice(m), gFwdUp.slice(m), gUp.slice(m)) =
				greenFromUdV_timedisplaced(eye_UdV, UdVStorageUp[n]);
		tie(ignore, gBwdDn.slice(m), gFwdDn.slice(m), gDn.slice(m)) =
				greenFromUdV_timedisplaced(eye_UdV, UdVStorageDn[n]);
		UdVStorageUp[n] = eye_UdV;
		UdVStorageDn[n] = eye_UdV;
		for (unsigned l = n; l >= 1; --l) {
			updateInSlice(l*s);
			for (unsigned k = l*s - 1; k >= (l-1)*s + 1; --k) {
				wrapDownGreen(k + 1, gUp, Spin::Up);
				wrapDownGreen(k + 1, gDn, Spin::Down);
				updateInSlice(k);
			}
			//TODO: this will also compute the Green function at k=0, which technically is not necessary
			//but sensible for the following sweep up
			//TODO: alternatively just copy the k=m Green function to k=0  -- would that be up-to-date?
			advanceDownGreen(l, UdVStorageUp, gUp, gFwdUp, gBwdUp, Spin::Up);
			advanceDownGreen(l, UdVStorageDn, gDn, gFwdDn, gBwdDn, Spin::Down);
		}
		lastSweepDir = SweepDirection::Down;
	} else if (lastSweepDir == SweepDirection::Down) {
//		debugCheckGreenFunctions();
//		debugCheckBeforeSweepUp();
		//We need to have computed the Green function for time slice k=0 so that the first
		//wrap-up step is correct.
		for (unsigned k = 1; k <= s-1; ++k) {
			wrapUpGreen(k - 1, gUp, Spin::Up);
			wrapUpGreen(k - 1, gDn, Spin::Down);
			updateInSlice(k);
		}
		//set storage at k=0 to unity for the upcoming sweep:
		UdVStorageUp[0] = eye_UdV;
		UdVStorageDn[0] = eye_UdV;
		for (unsigned l = 1; l < n; ++l) {
			advanceUpGreen(l-1, UdVStorageUp, gUp, gFwdUp, gBwdUp, Spin::Up);
			advanceUpGreen(l-1, UdVStorageDn, gDn, gFwdDn, gBwdDn, Spin::Down);
			updateInSlice(l*s);
			advanceUpUpdateStorage(l - 1, UdVStorageUp, Spin::Up);
			advanceUpUpdateStorage(l - 1, UdVStorageDn, Spin::Down);
			for (unsigned k = l*s + 1; k <= l*s + (s-1); ++k) {
				wrapUpGreen(k - 1, gUp, Spin::Up);
				wrapUpGreen(k - 1, gDn, Spin::Down);
				updateInSlice(k);
			}
		}
		updateInSlice(n*s);
		advanceUpUpdateStorage(n - 1, UdVStorageUp, Spin::Up);
		advanceUpUpdateStorage(n - 1, UdVStorageDn, Spin::Down);
		lastSweepDir = SweepDirection::Up;
	}
//	std::cout << std::endl;		//DEBUG
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
	//FORMULA-TEST -- made no difference
//	num sum_doubleoccupancy = 0;
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for (unsigned site = 0; site < N; ++site) {
			//use diagonal elements of Green functions:
			sum_GiiUp += gUp(site, site, timeslice);
			sum_GiiDn += gDn(site, site, timeslice);
			sum_GiiUpDn += gUp(site, site, timeslice) * gDn(site, site, timeslice);
			//use nearest neighbor elements of Green functions:
			for (unsigned neighIndex = 0; neighIndex < z; ++neighIndex) {
				unsigned neigh = nearestNeigbors(neighIndex, site);
				sum_GneighUp += gUp(site, neigh, timeslice);
				sum_GneighDn += gDn(site, neigh, timeslice);
			}
			//FORMULA-TEST -- made no difference
//			sum_doubleoccupancy += (1 - gUp(site,site, timeslice))
//								 * (1 - gDn(site,site, timeslice));
		}
	}
	occUp = 1.0 - (1.0 / (N*m)) * sum_GiiUp;
	occDn = 1.0 - (1.0 / (N*m)) * sum_GiiDn;
	occTotal = occUp + occDn;

	//FORMULA-TEST -- made no difference
//	std::cout << (1.0 / (N*m)) * sum_doubleoccupancy;
	occDouble = 1.0 + (1.0 / (N*m)) * (sum_GiiUpDn - sum_GiiUp - sum_GiiDn);
//	std::cout << " vs. " << occDouble << std::endl;

	localMoment = occTotal - 2*occDouble;

//	ePotential = (U / (N*m)) * (sum_GiiUpDn + 0.5 * sum_GiiUp + 0.5 * sum_GiiDn);
	ePotential = U * ( 0.25 + (1.0 / (N*m)) * (sum_GiiUpDn - 0.5 * (sum_GiiUp + sum_GiiDn)) );

	//Note: chemical potential term included in kinetic energy:
	eKinetic   = (t / (N*m)) * (sum_GneighUp + sum_GneighDn) - mu * occTotal;
	eTotal = eKinetic + ePotential;

	//susceptibility
	auto sumTrace = [m](const CubeNum& green) {
		num sum = 0;
		for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
			sum += arma::trace(green.slice(timeslice));
		}
		return sum;
	};
	num sumTrGreenUp = sumTrace(gUp);
	num sumTrGreenDn = sumTrace(gDn);
	auto sumProdTrace = [m](const CubeNum& green1, const CubeNum& green2) {
		num sum = 0;
		for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
			sum += arma::trace(green1.slice(timeslice) * green2.slice(timeslice));
		}
		return sum;
	};
	num sumTrGreenDisplacedUp = sumProdTrace(gBwdUp, gFwdUp);
	num sumTrGreenDisplacedDn = sumProdTrace(gBwdDn, gFwdDn);
	num trGreenUp_0 = arma::trace(gUp.slice(m));			//g(beta) = g(0)
	num trGreenDn_0 = arma::trace(gDn.slice(m));
	suscq0 = dtau * (  (trGreenUp_0 - trGreenDn_0) * (sumTrGreenUp - sumTrGreenDn)
					 - (sumTrGreenDisplacedUp + sumTrGreenDisplacedDn)
					);

	// vector observables
	zcorr.zeros();
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		num gUp_00 = gUp(0,0, timeslice);
		num gDn_00 = gDn(0,0, timeslice);
		zcorr[0] += -2.0 * gUp_00 * gDn_00 + gUp_00 + gDn_00;
		for (unsigned siteJ = 1; siteJ < N; ++siteJ) {
			num gUp_0j = gUp(0,siteJ,     timeslice);
			num gDn_0j = gDn(0,siteJ,     timeslice);
			num gUp_jj = gUp(siteJ,siteJ, timeslice);
			num gDn_jj = gDn(siteJ,siteJ, timeslice);
			using std::pow;
			zcorr[siteJ] += gUp_00 * gUp_jj - gUp_00 * gDn_jj + gDn_00 * gDn_jj - gDn_00 * gUp_jj
					- pow(gUp_0j, 2) - pow(gDn_0j, 2);
		}
	}
	zcorr /= num(m);
}


unsigned DetHubbard::getNumberOfObservables() const {
	return obsCount;
}

num DetHubbard::obsNormalized(unsigned obsIndex) const {
	if (obsIndex < obsCount) {
		return obsValRefs[obsIndex];
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

unsigned DetHubbard::getNumberOfVectorObservables() const {
	return vecObsCount;
}

VecNum DetHubbard::vecObsNormalized(unsigned vecObsIndex) const {
	if (vecObsIndex < vecObsCount) {
		return vecObsValRefs[vecObsIndex];
	} else {
		throw WrongObsIndex(vecObsIndex, true);
	}
}

std::string DetHubbard::getVectorObservableName(unsigned vecObsIndex) const {
	if (vecObsIndex < vecObsCount) {
		return vecObsNames[vecObsIndex];
	} else {
		throw WrongObsIndex(vecObsIndex, true);
	}
}

std::string DetHubbard::getVectorObservableShort(unsigned vecObsIndex) const {
	if (vecObsIndex < vecObsCount) {
		return vecObsShorts[vecObsIndex];
	} else {
		throw WrongObsIndex(vecObsIndex, true);
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
				auxfield(site, timeslice) = +1;
			} else {
				auxfield(site, timeslice) = -1;
			}
		}
	}
}

void DetHubbard::setupPropTmat() {
	MatNum tmat = -mu * arma::eye(N, N);

	for (unsigned site = 0; site < N; ++site) {
		//hopping between nearest neighbors
		for (auto p = nearestNeigbors.begin_col(site);
				 p != nearestNeigbors.end_col(site); ++p) {
			tmat(*p, site) -= t;
		}
	}

	proptmat = computePropagator(dtau, tmat);
}

MatNum DetHubbard::computePropagator(num scalar, const MatNum& matrix) const {
	using namespace arma;

	VecNum eigval;
	MatNum eigvec;
	eig_sym(eigval, eigvec, matrix);

	return eigvec * diagmat(exp(-scalar * eigval)) * trans(eigvec);
}


inline MatNum DetHubbard::computeBmatNaive(unsigned k2, unsigned k1, Spin spinz,
		const MatInt& arbitraryAuxfield) const {
	using namespace arma;

	if (k2 == k1) {
		return eye(N, N);
	}

	assert(k2 > k1);
	assert(k2 <= m);
	//assert(n1 >= 0);

	num sign = num(int(spinz));

	//Propagator using the HS-field potential for the given timeslice
	auto singleTimeslicePropagator = [this, sign, arbitraryAuxfield](unsigned timeslice) -> MatNum {
		//the cast with conv_to is necessary here, else everything would result in integers!
		//-- an Armadillo bug IMHO
		return diagmat(exp(sign * alpha *
				conv_to<VecNum>::from(arbitraryAuxfield.col(timeslice)))) * proptmat;
	};

	MatNum B = singleTimeslicePropagator(k2);

	for (unsigned k = k2 - 1; k >= k1 + 1; --k) {
		B *= singleTimeslicePropagator(k);
	}

	return B;
}

inline MatNum DetHubbard::computeBmatNaive(unsigned n2, unsigned n1, Spin spinz) const {
	return computeBmatNaive(n2, n1, spinz, auxfield);
}

inline MatNum DetHubbard::computeGreenFunctionNaive(
		const MatNum& bTau0, const MatNum& bBetaTau) const {
	return arma::inv(arma::eye(N,N) + bTau0 * bBetaTau);
}

inline MatNum DetHubbard::computeGreenFunctionNaive(unsigned timeslice,
		Spin spinz) const {
	//TODO: should use stored B-matrices, for the timeslices that have not changed
	return computeGreenFunctionNaive(computeBmatNaive(timeslice, 0, spinz),
			                    computeBmatNaive(m, timeslice, spinz));
}


num DetHubbard::weightRatioGenericNaive(const MatInt& auxfieldBefore,
		const MatInt& auxfieldAfter) const {
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

inline num DetHubbard::weightRatioSingleFlip(unsigned site, unsigned timeslice) const {
	using std::exp;
	//TODO: possibly precompute the exponential factors (auxfield is either +/- 1), would require an if though.

	//exponential factors
	//results again do not seem to for the location of the -sign
	num expUp   = exp(-2.0 * alpha * num(auxfield(site, timeslice)));
	num expDown = exp( 2.0 * alpha * num(auxfield(site, timeslice)));

	num ratioUp   = 1.0 + (expUp   - 1.0) * (1.0 - gUp(site,site, timeslice));
	num ratioDown = 1.0 + (expDown - 1.0) * (1.0 - gDn(site,site, timeslice));

//	std::cout << ratioUp << " " << ratioDown << std::endl;

	return ratioUp * ratioDown;
}

inline void DetHubbard::updateGreenFunctionWithFlip(unsigned site, unsigned timeslice) {
	auto update = [this, site](MatNum& green, num deltaSite) {
		const MatNum& greenOld = green;		//reference
		MatNum greenNew = green;			//copy
		const MatNum oneMinusGreenOld = arma::eye(N,N) - greenOld;
		num divisor = 1.0 + deltaSite * oneMinusGreenOld(site, site);
		num greenFactor = deltaSite / divisor;
		for (unsigned y = 0; y < N; ++y) {
			for (unsigned x = 0; x < N; ++x) {
				greenNew(x, y) -=  greenOld(x, site) * greenFactor *
						oneMinusGreenOld(site, y);

				//experimental index swap below -- this seemed to be the choice in Santos 2003
//				greenNew(x, y) -=  greenOld(site, y) * greenFactor *
//						oneMinusGreenOld(x, site);
			}
		}
		green = greenNew;
	};

	using std::exp;
	update(gUp.slice(timeslice), exp(-2.0 * alpha * num(auxfield(site, timeslice))) - 1.0);
	update(gDn.slice(timeslice), exp(+2.0 * alpha * num(auxfield(site, timeslice))) - 1.0);
}



