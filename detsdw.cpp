/*
 * detsdw.cpp
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#include <cmath>
#include <numeric>
#include <functional>
#include <array>
#include <tuple>
#include <cassert>
#include <boost/assign/std/vector.hpp>    // 'operator+=()' for vectors
#include "observable.h"
#include "detsdw.h"
#include "exceptions.h"
#include "timing.h"
#include "checkarray.h"

#if defined(MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

//initial values for field components chosen from this range:
const num PhiLow = -1;
const num PhiHigh = 1;
//adjustment of phiDelta:
const num InitialPhiDelta = 0.5;
const uint32_t AccRatioAdjustmentSamples = 100;
const num phiDeltaGrowFactor = 1.01;
const num phiDeltaShrinkFactor = 0.99;

std::string cbmToString(CheckerboardMethod cbm) {
	switch (cbm) {
	case CB_NONE: return "NONE";
	case CB_SANTOS: return "santos";
	case CB_ASSAAD: return "assaad";
	case CB_ASSAAD_BERG: return "assaad_berg";
	default: return "INVALID_CHECKERBOARD_METHOD";
	}
}

std::unique_ptr<DetModel> createDetSDW(RngWrapper& rng, ModelParams pars) {
    //TODO: add checks
    pars = updateTemperatureParameters(pars);

    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "mu", "L", "r", "accRatio", "bc", "txhor", "txver", "tyhor",
    		"tyver", "rescale", "updateMethod", "repeatUpdateInSlice";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (pars.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

    std::string possibleBC[] = {"pbc", "apbc-x", "apbc-y", "apbc-xy"};
    bool bc_is_one_of_the_possible = false;
    for (const std::string& bc : possibleBC) {
        if (pars.bc == bc) bc_is_one_of_the_possible = true;
    }
    if (not bc_is_one_of_the_possible) {
        throw ParameterWrong("bc", pars.bc);
    }
    std::string possibleUpdateMethods[] = {"iterative", "woodbury", "delayed"};
    bool updateMethod_is_one_of_the_possible = false;
    for (const std::string& updateMethod : possibleUpdateMethods) {
        if (pars.updateMethod == updateMethod) updateMethod_is_one_of_the_possible = true;
    }
    if (not updateMethod_is_one_of_the_possible) {
        throw ParameterWrong("updateMethod", pars.updateMethod);
    }
    if (pars.specified.count("updateMethod") and pars.updateMethod == "delayed") {
    	if (not pars.specified.count("delaySteps")) {
    		throw ParameterMissing("delaySteps");
    	}
    	uint32_t N = std::pow(pars.L, 2);
    	if (pars.delaySteps <= 0 or pars.delaySteps > N) {
    		throw ParameterWrong("delaySteps", pars.delaySteps);
    	}
    }



    if (pars.checkerboard and pars.L % 2 != 0) {
        throw ParameterWrong("Checker board decomposition only supported for even linear lattice sizes");
    }
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x)   {                                           \
                                IF_NOT_POSITIVE(x) {                    \
                                    throw ParameterWrong(#x, pars.x);   \
                                }                                       \
                            }
    CHECK_POSITIVE(L);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

    CheckerboardMethod cbm = CB_NONE;
    if (pars.checkerboard) {
    	if (pars.checkerboardMethod == "santos") {
    		cbm = CB_SANTOS;
    	} else if (pars.checkerboardMethod == "assaad") {
    		cbm = CB_ASSAAD;
    	} else if (pars.checkerboardMethod == "assaad_berg") {
    		cbm = CB_ASSAAD_BERG;
    	} else {
    		throw ParameterWrong("checkerboardMethod", pars.checkerboardMethod);
    	}
    }

    //since pars is not a constant expression, we need this stupid if:
    if (pars.timedisplaced == true and cbm == CB_NONE) {
        return std::unique_ptr<DetModel>(new DetSDW<true,CB_NONE>(rng, pars));
    } else
    if (pars.timedisplaced == true and cbm == CB_SANTOS) {
        return std::unique_ptr<DetModel>(new DetSDW<true,CB_SANTOS>(rng, pars));
    } else
    if (pars.timedisplaced == true and cbm == CB_ASSAAD) {
    	return std::unique_ptr<DetModel>(new DetSDW<true,CB_ASSAAD>(rng, pars));
    } else
    if (pars.timedisplaced == true and cbm == CB_ASSAAD_BERG) {
    	return std::unique_ptr<DetModel>(new DetSDW<true,CB_ASSAAD_BERG>(rng, pars));
    } else
    if (pars.timedisplaced == false and cbm == CB_NONE) {
        return std::unique_ptr<DetModel>(new DetSDW<false,CB_NONE>(rng, pars));
    } else
    if (pars.timedisplaced == false and cbm == CB_SANTOS) {
        return std::unique_ptr<DetModel>(new DetSDW<false,CB_SANTOS>(rng, pars));
    } else
    if (pars.timedisplaced == false and cbm == CB_ASSAAD) {
    	return std::unique_ptr<DetModel>(new DetSDW<false,CB_ASSAAD>(rng, pars));
    } else
    if (pars.timedisplaced == false and cbm == CB_ASSAAD_BERG) {
    	return std::unique_ptr<DetModel>(new DetSDW<false,CB_ASSAAD_BERG>(rng, pars));
    }
    else {
        //this can't be reached
        //return 0;
        return std::unique_ptr<DetModel>();
    }
}

template<bool TD, CheckerboardMethod CB>
DetSDW<TD,CB>::DetSDW(RngWrapper& rng_, const ModelParams& pars) :
        DetModelGC<1,cpx,TD>(pars, 4 * pars.L*pars.L),
        eye4cpx(arma::eye(4,4), arma::zeros(4,4)),
        rng(rng_),
        checkerboard(pars.checkerboard),
        checkerboardMethod(pars.checkerboardMethod),
        L(pars.L), N(L*L), r(pars.r),
        txhor(pars.txhor), txver(pars.txver), tyhor(pars.tyhor), tyver(pars.tyver),
        mu(pars.mu),
        c(1), u(1), lambda(1), //TODO: make these controllable by parameter
        bc(PBC), updateMethod(ITERATIVE), delaySteps(pars.delaySteps),
        rescale(pars.rescale), rescaleInterval(pars.rescaleInterval),
        rescaleGrowthFactor(pars.rescaleGrowthFactor), rescaleShrinkFactor(pars.rescaleShrinkFactor),
        acceptedRescales(0), attemptedRescales(0),
        repeatUpdateInSlice(pars.repeatUpdateInSlice),
        hopHor(), hopVer(), sinhHopHor(), sinhHopVer(), coshHopHor(), coshHopVer(),
        sinhHopHorHalf(), sinhHopVerHalf(), coshHopHorHalf(), coshHopVerHalf(),
        spaceNeigh(L), timeNeigh(m),
        propK(), propKx(propK[XBAND]), propKy(propK[YBAND]),
        propK_half(), propKx_half(propK_half[XBAND]), propKy_half(propK_half[YBAND]),
        propK_half_inv(), propKx_half_inv(propK_half_inv[XBAND]), propKy_half_inv(propK_half_inv[YBAND]),
        g(green[0]), //gFwd(greenFwd[0]), gBwd(greenBwd[0]),
        phi0(N, m+1), phi1(N, m+1), phi2(N, m+1), phiCosh(N, m+1), phiSinh(N, m+1),
        phiDelta(InitialPhiDelta),
        targetAccRatioLocal(pars.accRatio), lastAccRatioLocal(0), accRatioLocalRA(AccRatioAdjustmentSamples),
        performedSweeps(0),
        normPhi(0), meanPhi(), normMeanPhi(0), sdwSusc(0),
        kOcc(), kOccX(kOcc[XBAND]), kOccY(kOcc[YBAND]),
//        kOccImag(), kOccXimag(kOccImag[XBAND]), kOccYimag(kOccImag[YBAND]),
        occ(), occX(occ[XBAND]), occY(occ[YBAND]),
//        occImag(), occXimag(occImag[XBAND]), occYimag(occImag[YBAND]),
        pairPlusMax(0.0), pairMinusMax(0.0), //pairPlusMaximag(0.0), pairMinusMaximag(0.0),
        pairPlus(), pairMinus(), //pairPlusimag(), pairMinusimag(),
        fermionEkinetic(0), fermionEcouple(0),// fermionEkinetic_imag(0), fermionEcouple_imag(0),
        dud(N, delaySteps)
{
	assert((pars.checkerboard and CB != CB_NONE) or (not pars.checkerboard and CB == CB_NONE));
	assert(not pars.checkerboard or (pars.checkerboardMethod == cbmToString(CB)));

    if (pars.bc == "pbc") {
        bc = PBC;
    } else if (pars.bc == "apbc-x") {
        bc = APBC_X;
    } else if (pars.bc == "apbc-y") {
        bc = APBC_Y;
    } else if (pars.bc == "apbc-xy") {
        bc = APBC_XY;
    } else {
        // "safe default"
        bc = PBC;
    }
    if (pars.updateMethod == "iterative") {
    	updateMethod = ITERATIVE;
    } else if (pars.updateMethod == "woodbury") {
    	updateMethod = WOODBURY;
    } else if (pars.updateMethod == "delayed") {
    	updateMethod = DELAYED;
    } else {
        // "safe default"
    	updateMethod = ITERATIVE;
    }
    setupRandomPhi();

    //hopping constants. These are the t_ij in sum_<i,j> -t_ij c^+_i c_j
    //So for actual calculations an additional minus-sign needs to be included.
    //In the case of anti-periodic boundaries between i and j, another extra minus-sign
    //must be added.
    hopHor[XBAND] = txhor;
    hopVer[XBAND] = txver;
    hopHor[YBAND] = tyhor;
    hopVer[YBAND] = tyver;
    //precalculate hyperbolic functions, used in checkerboard decomposition
    using std::sinh; using std::cosh;
    num dtauHere = dtau;                // to fix capture issues
    for_each_band( [this, dtauHere](Band band) {
        sinhHopHor[band] = sinh(-dtauHere * hopHor[band]);
        coshHopHor[band] = cosh(-dtauHere * hopHor[band]);
        sinhHopVer[band] = sinh(-dtauHere * hopVer[band]);
        coshHopVer[band] = cosh(-dtauHere * hopVer[band]);
        sinhHopHorHalf[band] = sinh(-0.5*dtauHere * hopHor[band]);
        coshHopHorHalf[band] = cosh(-0.5*dtauHere * hopHor[band]);
        sinhHopVerHalf[band] = sinh(-0.5*dtauHere * hopVer[band]);
        coshHopVerHalf[band] = cosh(-0.5*dtauHere * hopVer[band]);
    } );

    setupPropK();

    setupUdVStorage_skeleton(sdwComputeBmat(this));

    using std::cref;
    using namespace boost::assign;
    obsScalar += ScalarObservable(cref(normPhi), "normPhi", "np"),
    		ScalarObservable(cref(normMeanPhi), "normMeanPhi", "nmp"),
            ScalarObservable(cref(sdwSusc), "sdwSusceptibility", "sdwsusc"),
            ScalarObservable(cref(pairPlusMax), "pairPlusMax", "ppMax"),
            ScalarObservable(cref(pairMinusMax), "pairMinusMax", "pmMax"),
            //ScalarObservable(cref(pairPlusMaximag), "pairPlusMaximag", "ppMaximag"),
            //ScalarObservable(cref(pairMinusMaximag), "pairMinusMaximag", "pmMaximag"),
            ScalarObservable(cref(fermionEkinetic), "fermionEkinetic", "fEkin"),
            //ScalarObservable(cref(fermionEkinetic_imag), "fermionEkineticimag", "fEkinimag"),
            ScalarObservable(cref(fermionEcouple), "fermionEcouple", "fEcouple");
            //ScalarObservable(cref(fermionEcouple_imag), "fermionEcoupleimag", "fEcoupleimag");

    kOccX.zeros(N);
    kOccY.zeros(N);
    obsVector += VectorObservable(cref(kOccX), N, "kOccX", "nkx"),
            VectorObservable(cref(kOccY), N, "kOccY", "nky");
//    kOccXimag.zeros(N);
//    kOccYimag.zeros(N);
//    obsVector += VectorObservable(cref(kOccXimag), N, "kOccXimag", "nkximag"),
//            VectorObservable(cref(kOccYimag), N, "kOccYimag", "nkyimag");

    occX.zeros(N);
    occY.zeros(N);
    obsVector += VectorObservable(cref(occX), N, "occX", "nx"),
            VectorObservable(cref(occY), N, "occY", "ny");
//    occXimag.zeros(N);
//    occYimag.zeros(N);
//    obsVector += VectorObservable(cref(occXimag), N, "occXimag", "nximag"),
//            VectorObservable(cref(occYimag), N, "occYimag", "nyimag");

    //attention:
    // these do not have valid entries for site 0
    pairPlus.zeros(N);
    pairMinus.zeros(N);
//    pairPlusimag.zeros(N);
//    pairMinusimag.zeros(N);
    obsVector += VectorObservable(cref(pairPlus), N, "pairPlus", "pp"),
            VectorObservable(cref(pairMinus), N, "pairMinus", "pm");
//            VectorObservable(cref(pairPlusimag), N, "pairPlusimag", "ppimag"),
//            VectorObservable(cref(pairMinusimag), N, "pairMinusimag", "pmimag");
//
}

template<bool TD, CheckerboardMethod CB>
DetSDW<TD,CB>::~DetSDW() {
}

template<bool TD, CheckerboardMethod CB>
uint32_t DetSDW<TD,CB>::getSystemN() const {
    return N;
}

template<bool TD, CheckerboardMethod CB>
MetadataMap DetSDW<TD,CB>::prepareModelMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
    meta["model"] = "sdw";
    meta["checkerboard"] = (CB ? "true" : "false");
    if (CB) {
    	meta["checkerboardMethod"] = checkerboardMethod;
    }
    meta["updateMethod"] = updateMethodstr(updateMethod);
    if (updateMethod == DELAYED) {
    	META_INSERT(delaySteps);
    }
    meta["timedisplaced"] = (TD ? "true" : "false");
    if (bc == PBC) {
          meta["bc"] = "pbc";
    } else if (bc == APBC_X) {
          meta["bc"] = "apbc-x";
    } else if (bc == APBC_Y) {
          meta["bc"] = "apbc-y";
    } else if (bc == APBC_XY) {
          meta["bc"] = "apbc-xy";
    }
    META_INSERT(targetAccRatioLocal);
    META_INSERT(r);
    META_INSERT(txhor);
    META_INSERT(txver);
    META_INSERT(tyhor);
    META_INSERT(tyver);
    META_INSERT(mu);
    META_INSERT(L);
    META_INSERT(d);
    META_INSERT(N);
    META_INSERT(beta);
    META_INSERT(m);
    META_INSERT(dtau);
    META_INSERT(s);
    META_INSERT(rescale);
    if (rescale) {
    	META_INSERT(rescaleInterval);
    	META_INSERT(rescaleGrowthFactor);
    	META_INSERT(rescaleShrinkFactor);
    }
    META_INSERT(repeatUpdateInSlice);
#undef META_INSERT
    return meta;
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::initMeasurements() {
    timing.start("sdw-measure");

    //normphi, meanPhi
    normPhi = 0;
    meanPhi.zeros();
    normMeanPhi = 0;

    //sdw-susceptibility
    sdwSusc = 0;

    //fermion occupation number -- real space
    occX.zeros(N);
    occY.zeros(N);

    //fermion occupation number -- k-space
    kOccX.zeros(N);
    kOccY.zeros(N);

    //equal-time pairing-correlations
    pairPlus.zeros(N);
    pairMinus.zeros(N);

    // Fermionic energy contribution
    fermionEkinetic = 0;
    fermionEcouple = 0;

    timing.stop("sdw-measure");
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::measure(uint32_t timeslice) {
	timing.start("sdw-measure");

	MatCpx gshifted = shiftGreenSymmetric();

	//normphi, meanPhi, sdw-susceptibility
	for (uint32_t site = 0; site < N; ++site) {
		Phi phi_site;
		phi_site[0] = phi0(site, timeslice);
		phi_site[1] = phi1(site, timeslice);
		phi_site[2] = phi2(site, timeslice);

		meanPhi += phi_site;
		normPhi += arma::norm(phi_site, 2);

		//sdw-susceptibility will be calculated in finishMeasurements()
	}

	//fermion occupation number -- real space
	for (uint32_t i = 0; i < N; ++i) {
		occX[i] += std::real(gshifted(i, i) + gshifted(i+N, i+N));
		occY[i] += std::real(gshifted(i+2*N, i+2*N) + gshifted(i+3*N, i+3*N));
	}

	//fermion occupation number -- k-space
    static const num pi = M_PI;
    //offset k-components for antiperiodic bc
    num offset_x = 0.0;
    num offset_y = 0.0;
    if (bc == APBC_X or bc == APBC_XY) {
        offset_x = 0.5;
    }
    if (bc == APBC_Y or bc == APBC_XY) {
        offset_y = 0.5;
    }
    for (uint32_t ksite = 0; ksite < N; ++ksite) {
        uint32_t ksitey = ksite / L;
        uint32_t ksitex = ksite % L;
        num ky = -pi + (num(ksitey) + offset_y) * 2*pi / num(L);
        num kx = -pi + (num(ksitex) + offset_x) * 2*pi / num(L);

        for (uint32_t i = 0; i < N; ++i) {
            num iy = num(i / L);
            num ix = num(i % L);
            for (uint32_t j = 0; j  < N; ++j) {
                num jy = num(j / L);
                num jx = num(j % L);

                num argument = kx * (ix - jx) + ky * (iy - jy);
                cpx phase = std::exp(cpx(0, argument));

                cpx green_x_up   = gshifted(i, j);
                cpx green_x_down = gshifted(i + N, j + N);
                cpx green_y_up   = gshifted(i + 2*N, j + 2*N);
                cpx green_y_down = gshifted(i + 3*N, j + 3*N);

                cpx x_cpx = phase * (green_x_up + green_x_down);
                cpx y_cpx = phase * (green_y_up + green_y_down);

                kOccX[ksite] += std::real(x_cpx);
                kOccY[ksite] += std::real(y_cpx);
            }
        }
    }

    //equal-time pairing-correlations
    //-------------------------------

    //helper to access the green function for the current time-slice
    //(which used to be "l")
    // *1 is for the row index,
    // *2 is for the column index
    auto gl = [this, &gshifted](uint32_t site1, Band band1, Spin spin1,
    							uint32_t site2, Band band2, Spin spin2) -> cpx {
    	return gshifted(site1 + 2*N*band1 + N*spin1,
    			 	 	site2 + 2*N*band2 + N*spin2);
    };

    for (uint32_t i = 0; i < N; ++i) {
    	//            checkarray<std::tuple<uint32_t,uint32_t>, 2> sitePairs = {
    	//                    std::make_tuple(i, 0), std::make_tuple(0, i)
    	//            };
    	//compiler-compatibilty fix
    	std::tuple<uint32_t,uint32_t> sitePairs[2] = {
    		std::tuple<uint32_t,uint32_t>(i, 0),
    		std::tuple<uint32_t,uint32_t>(0, i)
    	};

    	cpx pairPlusCpx(0, 0);
    	cpx pairMinusCpx(0, 0);

    	for (auto sites : sitePairs) {
    		uint32_t siteA = std::get<0>(sites);
    		uint32_t siteB = std::get<1>(sites);

    		// the following two unwieldy sums have been evaluated with the Mathematica
    		// notebook pairing-corr.nb (and they match the terms calculated by hand on paper)
    		pairPlusCpx += cpx(-4.0, 0) * (
    				gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINDOWN) -
    				gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINUP) +
    				gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINDOWN) -
    				gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINUP) +
    				gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINDOWN) -
    				gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINUP) +
    				gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINDOWN) -
    				gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINUP)
    		);

    		pairMinusCpx += cpx(-4.0, 0) * (
    				gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINDOWN) -
    				gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINUP) -
    				gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINDOWN) +
    				gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINUP) -
    				gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINDOWN) +
    				gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINUP) +
    				gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINDOWN) -
    				gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINUP)
    		);
    	}

    	pairPlus[i] += std::real(pairPlusCpx);
    	//pairPlusimag[i] += std::imag(pairPlusCpx);
    	pairMinus[i] += std::real(pairMinusCpx);
    	//pairMinusimag[i] += std::imag(pairMinusCpx);
    }

    // Fermionic energy contribution
    // -----------------------------
    auto glij = [this, &gshifted](uint32_t site1, uint32_t site2, Band band, Spin spin) -> cpx {
    	return gshifted(site1 + 2*N*band + N*spin,
    			 	 	site2 + 2*N*band + N*spin);
    };
    for (uint32_t i = 0; i < N; ++i) {
    	//TODO: write in a nicer fashion using hopping-array as used in the checkerboard branch
    	Spin spins[] = {SPINUP, SPINDOWN};
    	for (auto spin: spins) {
    		cpx e = cpx(txhor,0) * glij(i, spaceNeigh(XPLUS, i), XBAND, spin)
          		  + cpx(txhor,0) * glij(i, spaceNeigh(XMINUS,i), XBAND, spin)
          		  + cpx(txver,0) * glij(i, spaceNeigh(YPLUS, i), XBAND, spin)
          		  + cpx(txver,0) * glij(i, spaceNeigh(YMINUS,i), XBAND, spin)
          		  + cpx(tyhor,0) * glij(i, spaceNeigh(XPLUS, i), YBAND, spin)
          		  + cpx(tyhor,0) * glij(i, spaceNeigh(XMINUS,i), YBAND, spin)
          		  + cpx(tyver,0) * glij(i, spaceNeigh(YPLUS, i), YBAND, spin)
          		  + cpx(tyver,0) * glij(i, spaceNeigh(YMINUS,i), YBAND, spin);
    		fermionEkinetic += std::real(e);
    		//fermionEkinetic_imag += std::imag(e);
    	}
    }
    for (uint32_t i = 0; i < N; ++i) {
    	auto glbs = [this, i, gshifted](Band band1, Spin spin1,
    									Band band2, Spin spin2) -> cpx {
    		return gshifted(i + 2*N*band1 + N*spin1,
    						i + 2*N*band2 + N*spin2);
    	};

    	//factors for different combinations of spins
    	//overall factor of -1 included
    	cpx up_up(-phi2(i,timeslice), 0);
    	cpx up_dn(-phi0(i,timeslice), +phi1(i,timeslice));
    	cpx dn_up(-phi0(i,timeslice), -phi1(i,timeslice));
    	cpx dn_dn(+phi2(i,timeslice), 0);

		cpx e = up_up * (glbs(XBAND, SPINUP, YBAND, SPINUP) +
						 glbs(YBAND, SPINUP, XBAND, SPINUP))
			  + up_dn * (glbs(XBAND, SPINUP, YBAND, SPINDOWN) +
						 glbs(YBAND, SPINUP, XBAND, SPINDOWN))
			  + dn_up * (glbs(XBAND, SPINDOWN, YBAND, SPINUP) +
						 glbs(YBAND, SPINDOWN, XBAND, SPINUP))
			  + dn_dn * (glbs(XBAND, SPINDOWN, YBAND, SPINDOWN) +
						 glbs(YBAND, SPINDOWN, XBAND, SPINDOWN));

    	fermionEcouple += std::real(e);
    	//fermionEcouple_imag += std::imag(e);
    }

	timing.stop("sdw-measure");
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::finishMeasurements() {
	//normphi, meanPhi, sdw-susceptibility
	normPhi /= num(N * m);
	meanPhi /= num(N * m);
	normMeanPhi = arma::norm(meanPhi, 2);

	Phi phi_0;
	phi_0[0] = phi0(0, m);
	phi_0[1] = phi1(0, m);
	phi_0[2] = phi2(0, m);
	sdwSusc = 0;
	for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
		for (uint32_t site = 0; site < N; ++site) {
			sdwSusc += ( phi_0[0] * phi0(site,timeslice)
					   + phi_0[1] * phi1(site,timeslice)
					   + phi_0[2] * phi2(site,timeslice)
					   );
		}
	}
	sdwSusc *= dtau;


	//fermion occupation number -- real space
    occX /= num(m * N);
    occY /= num(m * N);

	//fermion occupation number -- k-space
    for (uint32_t ksite = 0; ksite < N; ++ksite) {
        // add 2.0 and not 1.0 because spin is included
        kOccX[ksite] = 2.0 - kOccX[ksite] / num(m * N);
        kOccY[ksite] = 2.0 - kOccY[ksite] / num(m * N);
    }

    //equal-time pairing-correlations
    //-------------------------------
    pairPlus /= m;
    //pairPlusimag /= m;
    pairMinus /= m;
    //pairMinusimag /= m;
    // sites around the maximum range L/2, L/2
    static const uint32_t numSitesFar = 9;
    uint32_t sitesfar[numSitesFar] = {
    		coordsToSite(L/2 - 1, L/2 - 1), coordsToSite(L/2, L/2 - 1), coordsToSite(L/2 + 1, L/2 - 1),
    		coordsToSite(L/2 - 1, L/2),     coordsToSite(L/2, L/2),     coordsToSite(L/2 + 1, L/2),
    		coordsToSite(L/2 - 1, L/2 + 1), coordsToSite(L/2, L/2 + 1), coordsToSite(L/2 + 1, L/2 + 1)
    };
    pairPlusMax = 0;
    //pairPlusMaximag = 0;
    pairMinusMax = 0;
    //pairMinusMaximag = 0;
    for (uint32_t i : sitesfar) {
    	pairPlusMax += pairPlus[i];
    	//pairPlusMaximag += pairPlusimag[i];
    	pairMinusMax += pairMinus[i];
    	//pairMinusMaximag += pairMinusimag[i];
    }
    pairPlusMax /= numSitesFar;
    //pairPlusMaximag /= numSitesFar;
    pairMinusMax /= numSitesFar;
    //pairMinusMaximag /= numSitesFar;

    // Fermionic energy contribution
    // -----------------------------
    fermionEkinetic /= num(m*N);
    fermionEcouple /= num(m*N);
}


//template<bool TD, CheckerboardMethod CB>
//void DetSDW<TD,CB>::measure() {
//    timing.start("sdw-measure");
//
//    shiftGreenSymmetric();
//
//    Phi meanPhi;
//    meanPhi[0] = averageWholeSystem(phi0, 0.0);
//    meanPhi[1] = averageWholeSystem(phi1, 0.0);
//    meanPhi[2] = averageWholeSystem(phi2, 0.0);
//    normPhi = arma::norm(meanPhi, 2);
//
//
//    //fermion occupation number -- real space
//    //probably not very interesting data
//    occX.zeros(N);
//    occY.zeros(N);
////    occXimag.zeros(N);
////    occYimag.zeros(N);
//    for (uint32_t l = 1; l <= m; ++l) {
//        for (uint32_t i = 0; i < N; ++i) {
//            occX[i] += std::real(g.slice(l)(i, i) + g.slice(l)(i+N, i+N));
//            occY[i] += std::real(g.slice(l)(i+2*N, i+2*N) + g.slice(l)(i+3*N, i+3*N));
////            occXimag[i] += std::imag(g.slice(l)(i, i) + g.slice(l)(i+N, i+N));
////            occYimag[i] += std::imag(g.slice(l)(i+2*N, i+2*N) + g.slice(l)(i+3*N, i+3*N));
//        }
//    }
//    //not working with icpc 13.1:
////  using std::ref;
////  for (VecNum& occ : {ref(occX), ref(occY), ref(occXimag), ref(occYimag)}) {
////      occ /= num(m) * num(N);
////  }
//    occX /= num(m) * num(N);
//    occY /= num(m) * num(N);
////    occXimag /= num(m) * num(N);
////    occYimag /= num(m) * num(N);
//
//
//
//    //fermion occupation number -- k-space
//    static const num pi = M_PI;
//    //offset k-components for antiperiodic bc
//    num offset_x = 0.0;
//    num offset_y = 0.0;
//    if (bc == APBC_X or bc == APBC_XY) {
//        offset_x = 0.5;
//    }
//    if (bc == APBC_Y or bc == APBC_XY) {
//        offset_y = 0.5;
//    }
//    for (uint32_t ksite = 0; ksite < N; ++ksite) {
//        //try a slightly alternative approach..
//        uint32_t ksitey = ksite / L;
//        uint32_t ksitex = ksite % L;
//        num ky = -pi + (num(ksitey) + offset_y) * 2*pi / num(L);
//        num kx = -pi + (num(ksitex) + offset_x) * 2*pi / num(L);
//
//        kOccX[ksite] = 0.0;
//        kOccY[ksite] = 0.0;
////        kOccXimag[ksite] = 0.0;
////        kOccYimag[ksite] = 0.0;
//
//        for (uint32_t i = 0; i < N; ++i) {
//            num iy = num(i / L);
//            num ix = num(i % L);
//            for (uint32_t j = 0; j  < N; ++j) {
//                num jy = num(j / L);
//                num jx = num(j % L);
//
//                num argument = kx * (ix - jx) + ky * (iy - jy);
//                cpx phase = std::exp(cpx(0, argument));
//
//                for (uint32_t l = 1; l <= m; ++l) {
//                    cpx green_x_up   = g.slice(l)(i, j);
//                    cpx green_x_down = g.slice(l)(i + N, j + N);
//                    cpx green_y_up   = g.slice(l)(i + 2*N, j + 2*N);
//                    cpx green_y_down = g.slice(l)(i + 3*N, j + 3*N);
//
//                    cpx x_cpx = phase * (green_x_up + green_x_down);
//                    cpx y_cpx = phase * (green_y_up + green_y_down);
//
//                    kOccX[ksite] += std::real(x_cpx);
//                    kOccY[ksite] += std::real(y_cpx);
////                    kOccXimag[ksite] += std::imag(x_cpx);
////                    kOccYimag[ksite] += std::imag(y_cpx);
//                }
//            }
//        }
//
//        // add 2.0 and not 1.0 because spin is included
//        kOccX[ksite] = 2.0 - kOccX[ksite] / num(m * N);
//        kOccY[ksite] = 2.0 - kOccY[ksite] / num(m * N);
////        kOccXimag[ksite] =  -kOccXimag[ksite] / num(m * N);
////        kOccYimag[ksite] =  -kOccYimag[ksite] / num(m * N);
//    }
//
//    //sdw-susceptibility
//    uint32_t mm = m;
//    sdwSusc = dtau * sumWholeSystem( [this, mm](uint32_t site, uint32_t timeslice) {
//                                            return phi0(site, timeslice) * phi0(0, mm)
//                                                 + phi1(site, timeslice) * phi1(0, mm)
//                                                 + phi2(site, timeslice) * phi2(0, mm);
//                                        },
//                                    0.0);
//
//    //equal-time pairing-correlations
//    //-------------------------------
//    pairPlus.zeros(N);
//    pairMinus.zeros(N);
////    pairPlusimag.zeros(N);
////    pairMinusimag.zeros(N);
//    for (uint32_t l = 1; l <= m; ++l) {
//        //helper to access the green function
//        // *1 is for the row index,
//        // *2 is for the column index
//        auto gl = [this, l](uint32_t site1, Band band1, Spin spin1,
//                           uint32_t site2, Band band2, Spin spin2) -> cpx {
//            return g.slice(l)(site1 + 2*N*band1 + N*spin1,
//                              site2 + 2*N*band2 + N*spin2);
//        };
//
//        for (uint32_t i = 0; i < N; ++i) {
////            checkarray<std::tuple<uint32_t,uint32_t>, 2> sitePairs = {
////                    std::make_tuple(i, 0), std::make_tuple(0, i)
////            };
//        	//compiler-compatibilty fix
//        	std::tuple<uint32_t,uint32_t> sitePairs[2] = {
//        			std::tuple<uint32_t,uint32_t>(i, 0),
//        			std::tuple<uint32_t,uint32_t>(0, i)
//        	};
//
//            cpx pairPlusCpx(0, 0);
//            cpx pairMinusCpx(0, 0);
//
//            for (auto sites : sitePairs) {
//                uint32_t siteA = std::get<0>(sites);
//                uint32_t siteB = std::get<1>(sites);
//
//                // the following two unwieldy sums have been evaluated with the Mathematica
//                // notebook pairing-corr.nb (and they match the terms calculated by hand on paper)
//                pairPlusCpx += cpx(-4.0, 0) * (
//                        gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINDOWN) -
//                        gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINUP) +
//                        gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINDOWN) -
//                        gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINUP) +
//                        gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINDOWN) -
//                        gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINUP) +
//                        gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINDOWN) -
//                        gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINUP)
//                );
//
//                pairMinusCpx += cpx(-4.0, 0) * (
//                        gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINDOWN) -
//                        gl(siteA, XBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, XBAND, SPINUP) -
//                        gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINDOWN) +
//                        gl(siteA, XBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, XBAND, SPINUP, siteB, YBAND, SPINUP) -
//                        gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINDOWN) +
//                        gl(siteA, YBAND, SPINDOWN, siteB, XBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, XBAND, SPINUP) +
//                        gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINUP)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINDOWN) -
//                        gl(siteA, YBAND, SPINDOWN, siteB, YBAND, SPINDOWN)*gl(siteA, YBAND, SPINUP, siteB, YBAND, SPINUP)
//                );
//            }
//
//            pairPlus[i] += std::real(pairPlusCpx);
////            pairPlusimag[i] += std::imag(pairPlusCpx);
//            pairMinus[i] += std::real(pairMinusCpx);
////            pairMinusimag[i] += std::imag(pairMinusCpx);
//        }
//    }
//    pairPlus /= m;
////    pairPlusimag /= m;
//    pairMinus /= m;
////    pairMinusimag /= m;
//
//    // sites around the maximum range L/2, L/2
//    static const uint32_t numSitesFar = 9;
//    uint32_t sitesfar[numSitesFar] = {
//            coordsToSite(L/2 - 1, L/2 - 1), coordsToSite(L/2, L/2 - 1), coordsToSite(L/2 + 1, L/2 - 1),
//            coordsToSite(L/2 - 1, L/2),     coordsToSite(L/2, L/2),     coordsToSite(L/2 + 1, L/2),
//            coordsToSite(L/2 - 1, L/2 + 1), coordsToSite(L/2, L/2 + 1), coordsToSite(L/2 + 1, L/2 + 1)
//    };
//    pairPlusMax = 0;
////    pairPlusMaximag = 0;
//    pairMinusMax = 0;
////    pairMinusMaximag = 0;
//    for (uint32_t i : sitesfar) {
//        pairPlusMax += pairPlus[i];
////        pairPlusMaximag += pairPlusimag[i];
//        pairMinusMax += pairMinus[i];
////        pairMinusMaximag += pairMinusimag[i];
//    }
//    pairPlusMax /= numSitesFar;
////    pairPlusMaximag /= numSitesFar;
//    pairMinusMax /= numSitesFar;
////    pairMinusMaximag /= numSitesFar;
//
//
//    // Fermionic energy contribution
//    // -----------------------------
//    fermionEkinetic = 0;
////    fermionEkinetic_imag = 0;
//    for (uint32_t l = 1; l <= m; ++l) {
//        auto glij = [this, l](uint32_t site1, uint32_t site2, Band band, Spin spin) -> cpx {
//            return g.slice(l)(site1 + 2*N*band + N*spin,
//                              site2 + 2*N*band + N*spin);
//        };
//        for (uint32_t i = 0; i < N; ++i) {
//            //TODO: write in a nicer fashion using hopping-array as used in the checkerboard branch
//            Spin spins[] = {SPINUP, SPINDOWN};
//            for (auto spin: spins) {
//                cpx e = cpx(txhor,0) * glij(i, spaceNeigh(XPLUS, i), XBAND, spin)
//                      + cpx(txhor,0) * glij(i, spaceNeigh(XMINUS,i), XBAND, spin)
//                      + cpx(txver,0) * glij(i, spaceNeigh(YPLUS, i), XBAND, spin)
//                      + cpx(txver,0) * glij(i, spaceNeigh(YMINUS,i), XBAND, spin)
//                      + cpx(tyhor,0) * glij(i, spaceNeigh(XPLUS, i), YBAND, spin)
//                      + cpx(tyhor,0) * glij(i, spaceNeigh(XMINUS,i), YBAND, spin)
//                      + cpx(tyver,0) * glij(i, spaceNeigh(YPLUS, i), YBAND, spin)
//                      + cpx(tyver,0) * glij(i, spaceNeigh(YMINUS,i), YBAND, spin);
//                fermionEkinetic += std::real(e);
////                fermionEkinetic_imag += std::imag(e);
//            }
//        }
//    }
//    fermionEkinetic /= num(m*N);
////    fermionEkinetic_imag /= num(m*N);
//
//    fermionEcouple = 0;
////    fermionEcouple_imag = 0;
//    for (uint32_t l = 1; l <= m; ++l) {
//        for (uint32_t i = 0; i < N; ++i) {
//            auto glbs = [this, l,i](Band band1, Spin spin1,
//                                    Band band2, Spin spin2) -> cpx {
//                return g.slice(l)(i + 2*N*band1 + N*spin1,
//                                  i + 2*N*band2 + N*spin2);
//            };
//
//            //factors for different combinations of spins
//            //overall factor of -1 included
//            cpx up_up(-phi2(i,l), 0);
//            cpx up_dn(-phi0(i,l), +phi1(i,l));
//            cpx dn_up(-phi0(i,l), -phi1(i,l));
//            cpx dn_dn(+phi2(i,l), 0);
//
//            cpx e = up_up * (glbs(XBAND, SPINUP, YBAND, SPINUP) +
//                             glbs(YBAND, SPINUP, XBAND, SPINUP))
//                  + up_dn * (glbs(XBAND, SPINUP, YBAND, SPINDOWN) +
//                             glbs(YBAND, SPINUP, XBAND, SPINDOWN))
//                  + dn_up * (glbs(XBAND, SPINDOWN, YBAND, SPINUP) +
//                             glbs(YBAND, SPINDOWN, XBAND, SPINUP))
//                  + dn_dn * (glbs(XBAND, SPINDOWN, YBAND, SPINDOWN) +
//                             glbs(YBAND, SPINDOWN, XBAND, SPINDOWN));
//            fermionEcouple += std::real(e);
////            fermionEcouple_imag += std::imag(e);
//        }
//
//    }
//    fermionEcouple /= num(m*N);
////    fermionEcouple_imag /= num(m*N);
//
//    timing.stop("sdw-measure");
//}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::setupRandomPhi() {
    for (uint32_t k = 1; k <= m; ++k) {
        for (uint32_t site = 0; site < N; ++site) {
            phi0(site, k) = rng.randRange(PhiLow, PhiHigh);
            phi1(site, k) = rng.randRange(PhiLow, PhiHigh);
            phi2(site, k) = rng.randRange(PhiLow, PhiHigh);
            num phiNorm = std::sqrt(std::pow(phi0(site, k), 2)
                                    + std::pow(phi1(site, k), 2)
                                    + std::pow(phi2(site, k), 2));
            phiCosh(site, k) = std::cosh(dtau * phiNorm);
            phiSinh(site, k) = std::sinh(dtau * phiNorm) / phiNorm;
        }
    }
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::setupPropK() {
    checkarray<checkarray<num,z>, 2> t;
    t[XBAND][XPLUS] = t[XBAND][XMINUS] = hopHor[XBAND];
    t[XBAND][YPLUS] = t[XBAND][YMINUS] = hopVer[XBAND];
    t[YBAND][XPLUS] = t[YBAND][XMINUS] = hopHor[YBAND];
    t[YBAND][YPLUS] = t[YBAND][YMINUS] = hopVer[YBAND];

//  for (auto band : {XBAND, YBAND}) {
    Band bands[2] = {XBAND, YBAND};
    for (Band band : bands) {
        MatNum k = -mu * arma::eye(N,N);
        for (uint32_t site = 0; site < N; ++site) {
            for (uint32_t dir = 0; dir < z; ++dir) {
                uint32_t neigh = spaceNeigh(dir, site);
                num hop = t[band][dir];

                uint32_t siteY = site / L;
                uint32_t siteX = site % L;
                if (bc == APBC_X or bc == APBC_XY) {
                    if ((siteX == 0 and dir == XMINUS) or (siteX == L-1 and dir == XPLUS)) {
                        //crossing x-boundary
                        hop *= -1;
                    }
                }
                if (bc == APBC_Y or bc == APBC_XY) {
                    if ((siteY == 0 and dir == YMINUS) or (siteY == L-1 and dir == YPLUS)) {
                        //crossing y-boundary
                        hop *= -1;
                    }
                }

                k(site, neigh) -= hop;
            }
        }
        propK[band] = computePropagator(dtau, k);

        propK_half[band] = computePropagator(dtau / 2.0, k);
        propK_half_inv[band] = computePropagator(-dtau / 2.0, k);
    }
}


template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::computeBmatSDW(uint32_t k2, uint32_t k1) {
	if (CB == CB_NONE) {
		timing.start("computeBmatSDW_direct");
		using arma::eye; using arma::zeros; using arma::diagmat;
		if (k2 == k1) {
			return MatCpx(eye(4*N,4*N), zeros(4*N,4*N));
		}
		assert(k2 > k1);
		assert(k2 <= m);

		//compute the matrix e^(-dtau*V_k) * e^(-dtau*K)
		auto singleTimesliceProp = [this](uint32_t k) -> MatCpx {
			timing.start("singleTimesliceProp_direct");
			MatCpx result(4*N, 4*N);

			//submatrix view helper for a 4N*4N matrix
			auto block = [&result, this](uint32_t row, uint32_t col) {
				return result.submat(row * N, col * N,
						(row + 1) * N - 1, (col + 1) * N - 1);
			};
			const auto& kphi0 = phi0.col(k);
			const auto& kphi1 = phi1.col(k);
			const auto& kphi2 = phi2.col(k);
			//      debugSaveMatrix(kphi0, "kphi0");
			//      debugSaveMatrix(kphi1, "kphi1");
			//      debugSaveMatrix(kphi2, "kphi2");
			const auto& kphiCosh = phiCosh.col(k);
			const auto& kphiSinh = phiSinh.col(k);
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

			//      debugSaveMatrix(arma::real(result), "emdtauVemdtauK_real");
			//      debugSaveMatrix(arma::imag(result), "emdtauVemdtauK_imag");
			timing.stop("singleTimesliceProp_direct");
			return result;
		};

		MatCpx result = singleTimesliceProp(k2);

		for (uint32_t k = k2 - 1; k > k1; --k) {
			result *= singleTimesliceProp(k);               // equivalent to: result = result * singleTimesliceProp(k);
		}

		timing.stop("computeBmatSDW_direct");

		return result;
	}
	else {
		// use the checkerboard routines to compute B by left-multiplying to unity
		using arma::eye; using arma::zeros;
		if (k2 == k1) {
			return MatCpx(eye(4*N,4*N), zeros(4*N,4*N));
		}
		assert(k2 > k1);
		assert(k2 <= m);
		MatCpx unity(eye(4*N, 4*N), zeros(4*N, 4*N));
		return checkerboardLeftMultiplyBmat(unity, k2, k1);
	}
}



template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
											 const Matrix&, Band, int, bool) {
	throw GeneralError("CB_NONE makes no sense for the checkerboard multiplication routines");
	//TODO change things so this codepath is not needed
	return MatCpx();
}


//neigh == XNEIGH:
//   subgroup == 0:  bonds (2*i_x, i_y)--(2*i_x + 1, i_y)
//   subgroup == 1:  bonds (2*i_x + 1, i_y)--(2*i_x + 2, i_y)
//neigh == YNEIGH:
//   subgroup == 0:  bonds (i_x, 2*i_y)--(i_x, 2*i_y + 1)
//   subgroup == 1:  bonds (i_x, 2*i_y + 1)--(i_x, 2*i_y + 2)
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_santos_applyBondFactorsLeft(Matrix& result, const NeighDir neigh, const uint32_t subgroup, const num ch, const num sh) {
	assert(subgroup == 0 or subgroup == 1);
	assert(neigh == XPLUS or neigh == YPLUS);
	arma::Row<cpx> new_row_i(N);
	for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
		for (uint32_t i2 = 0; i2 < L; ++i2) {
			uint32_t i;
			switch (neigh) {
			case XPLUS:
				i = this->coordsToSite(i1, i2);
				break;
			case YPLUS:
				i = this->coordsToSite(i2, i1);
				break;
			default: //should not be reached
				break;
			}
			uint32_t j = spaceNeigh(neigh, i);
			//change rows i and j of result
			num b_sh = sh;
			if ((bc == APBC_X or bc == APBC_XY) and neigh == XPLUS and i1 == L-1) {
				//crossed antiperiodic boundary
				b_sh *= -1;
			}
			else if ((bc == APBC_Y or bc == APBC_XY) and neigh == YPLUS and i1 == L-1) {
				//crossed antiperiodic boundary
				b_sh *= -1;
			}
			new_row_i     = ch * result.row(i) + b_sh * result.row(j);
			result.row(j) = b_sh * result.row(i) + ch * result.row(j);
			result.row(i) = new_row_i;
		}
	}
}



// with sign = +/- 1, band = XBAND|YBAND: set R := E^(sign * dtau * K_band) * A
// using the method described in R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
											 const Matrix& A, Band band, int sign, bool invertedCbOrder) {
	MatCpx result = A;      //can't avoid this copy

	if (not invertedCbOrder) {
		cb_santos_applyBondFactorsLeft(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
		cb_santos_applyBondFactorsLeft(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
		cb_santos_applyBondFactorsLeft(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
		cb_santos_applyBondFactorsLeft(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
	} else {
		cb_santos_applyBondFactorsLeft(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
		cb_santos_applyBondFactorsLeft(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
		cb_santos_applyBondFactorsLeft(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
		cb_santos_applyBondFactorsLeft(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
	}
	return result;
}

//subgroup == 0: plaquettes A = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
//   i = (2m, 2n), with m,n integer, 2m < L, 2n < L
//   j = i + XPLUS
//   k = i + YPLUS
//   l = k + XPLUS
//subgroup == 1: plaquettes B = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
//   i = (2m+1, 2n+1), with m,n integer, 2m+1 < L, 2n+1 < L
//   j = i + XPLUS
//   k = i + YPLUS
//   l = k + XPLUS
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_assaad_applyBondFactorsLeft(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver) {
	assert(subgroup == 0 or subgroup == 1);
	arma::Row<cpx> new_row_i(N);
	arma::Row<cpx> new_row_j(N);
	arma::Row<cpx> new_row_k(N);
	for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
		for (uint32_t i2 = subgroup; i2 < L; i2 += 2) {
			uint32_t i = this->coordsToSite(i1, i2);
			uint32_t j = spaceNeigh(XPLUS, i);
			uint32_t k = spaceNeigh(YPLUS, i);
			uint32_t l = spaceNeigh(XPLUS, k);
			//change rows i,j,k,l of result
			const arma::Row<cpx>& ri = result.row(i);
			const arma::Row<cpx>& rj = result.row(j);
			const arma::Row<cpx>& rk = result.row(k);
			const arma::Row<cpx>& rl = result.row(l);
			num b_sh_hor = sh_hor;
			num b_sh_ver = sh_ver;
			if ((bc == APBC_X or bc == APBC_XY) and i1 == L-1) {
				//this plaquette has horizontal boundary crossing bonds and APBC
				b_sh_hor *= -1;
			}
			if ((bc == APBC_Y or bc == APBC_XY) and i2 == L-1) {
				//this plaquette has vertical boundary crossing bonds and APBC
				b_sh_ver *= -1;
			}
			new_row_i     = ch_hor*ch_ver*ri + ch_ver*b_sh_hor*rj + ch_hor*b_sh_ver*rk + b_sh_hor*b_sh_ver*rl;
			new_row_j     = ch_ver*b_sh_hor*ri + ch_hor*ch_ver*rj + b_sh_hor*b_sh_ver*rk + ch_hor*b_sh_ver*rl;
			new_row_k     = ch_hor*b_sh_ver*ri + b_sh_hor*b_sh_ver*rj + ch_hor*ch_ver*rk + ch_ver*b_sh_hor*rl;
			result.row(l) = b_sh_hor*b_sh_ver*ri + ch_hor*b_sh_ver*rj + ch_ver*b_sh_hor*rk + ch_hor*ch_ver*rl;
			result.row(i) = new_row_i;
			result.row(j) = new_row_j;
			result.row(k) = new_row_k;
		}
	}
}

// with sign = +/- 1, band = XBAND|YBAND: set R := E^(sign * dtau * K_band) * A
// using the method described in F. F. Assaad, in Quantum Simulations Complex Many-Body Syst. From Theory to Algorithms, edited by J. Grotendorst, D. Marx, and A. Muramatsu (FZ-Jülich, Jülich, Germany, 2002).
template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
											 const Matrix& A, Band band, int sign, bool invertedCbOrder) {
	MatCpx result = A;      //can't avoid this copy

	if (not invertedCbOrder) {
		cb_assaad_applyBondFactorsLeft(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
		cb_assaad_applyBondFactorsLeft(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
	} else {
		cb_assaad_applyBondFactorsLeft(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
		cb_assaad_applyBondFactorsLeft(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
	}
	return result;
}

// with sign = +/- 1, band = XBAND|YBAND: set R := E^(sign * dtau * K_band) * A
// using the symmetric checkerboard break up
template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
											 const Matrix& A, Band band, int sign, bool) {
	MatCpx result = A;      //can't avoid this copy

	// perform the multiplication e^(+-dtau K_1/2) e^(+-dtau K_0) e^(+-dtau K_a/2) X
	cb_assaad_applyBondFactorsLeft(result, 1, coshHopHorHalf[band], sign * sinhHopHorHalf[band],
			                       coshHopVerHalf[band], sign * sinhHopVerHalf[band]);
	cb_assaad_applyBondFactorsLeft(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
	cb_assaad_applyBondFactorsLeft(result, 1, coshHopHorHalf[band], sign * sinhHopHorHalf[band],
								   coshHopVerHalf[band], sign * sinhHopVerHalf[band]);
	return result;
}

// with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to A * E^(sign * dtau * K_band)
template<bool TD, CheckerboardMethod CB>
template <class Matrix> inline
MatCpx DetSDW<TD,CB>::cbLMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder) {
	return cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB>(),
								  A, band, sign, invertedCbOrder);
}




template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
											 const Matrix&, Band, int, bool) {
	throw GeneralError("CB_NONE makes no sense for the checkerboard multiplication routines");
	return MatCpx();
}


//neigh == XNEIGH:
//   subgroup == 0:  bonds (2*i_x, i_y)--(2*i_x + 1, i_y)
//   subgroup == 1:  bonds (2*i_x + 1, i_y)--(2*i_x + 2, i_y)
//neigh == YNEIGH:
//   subgroup == 0:  bonds (i_x, 2*i_y)--(i_x, 2*i_y + 1)
//   subgroup == 1:  bonds (i_x, 2*i_y + 1)--(i_x, 2*i_y + 2)
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_santos_applyBondFactorsRight(Matrix& result, const NeighDir neigh, const uint32_t subgroup, const num ch, const num sh) {
	assert(subgroup == 0 or subgroup == 1);
	assert(neigh == XPLUS or neigh == YPLUS);
	arma::Col<cpx> new_col_i(N);
	for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
		for (uint32_t i2 = 0; i2 < L; ++i2) {
			uint32_t i;
			switch (neigh) {
			case XPLUS:
				i = this->coordsToSite(i1, i2);
				break;
			case YPLUS:
				i = this->coordsToSite(i2, i1);
				break;
			default: //should not be reached
				break;
			}
			uint32_t j = spaceNeigh(neigh, i);
			//change columns i and j of result
			num b_sh = sh;
			if ((bc == APBC_X or bc == APBC_XY) and neigh == XPLUS and i1 == L-1) {
				//crossed antiperiodic boundary
				b_sh *= -1;
			}
			else if ((bc == APBC_Y or bc == APBC_XY) and neigh == YPLUS and i1 == L-1) {
				//crossed antiperiodic boundary
				b_sh *= -1;
			}
			new_col_i     = ch * result.col(i) + b_sh * result.col(j);
			result.col(j) = b_sh * result.col(i) + ch * result.col(j);
			result.col(i) = new_col_i;
		}
	}
}



// with sign = +/- 1, band = XBAND|YBAND: return A * E^(sign * dtau * K_band)
// using the method described in R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
template<bool TD, CheckerboardMethod CB>
template <class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
											 const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    //order reversed wrt cbLMultHoppingExp
    if (not invertedCbOrder) {
    	cb_santos_applyBondFactorsRight(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
    	cb_santos_applyBondFactorsRight(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
    	cb_santos_applyBondFactorsRight(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
    	cb_santos_applyBondFactorsRight(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
    } else {
    	cb_santos_applyBondFactorsRight(result, XPLUS, 0, coshHopHor[band], sign * sinhHopHor[band]);
    	cb_santos_applyBondFactorsRight(result, YPLUS, 0, coshHopVer[band], sign * sinhHopVer[band]);
    	cb_santos_applyBondFactorsRight(result, XPLUS, 1, coshHopHor[band], sign * sinhHopHor[band]);
    	cb_santos_applyBondFactorsRight(result, YPLUS, 1, coshHopVer[band], sign * sinhHopVer[band]);
    }

    return result;
}

//subgroup == 0: plaquettes A = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
//   i = (2m, 2n), with m,n integer, 2m < L, 2n < L
//   j = i + XPLUS
//   k = i + YPLUS
//   l = k + XPLUS
//subgroup == 1: plaquettes B = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
//   i = (2m+1, 2n+1), with m,n integer, 2m+1 < L, 2n+1 < L
//   j = i + XPLUS
//   k = i + YPLUS
//   l = k + XPLUS
template<bool TD, CheckerboardMethod CB>
template<class Matrix>
void DetSDW<TD,CB>::cb_assaad_applyBondFactorsRight(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver) {
	assert(subgroup == 0 or subgroup == 1);
	arma::Col<cpx> new_col_i(N);
	arma::Col<cpx> new_col_j(N);
	arma::Col<cpx> new_col_k(N);
	for (uint32_t i1 = subgroup; i1 < L; i1 += 2) {
		for (uint32_t i2 = subgroup; i2 < L; i2 += 2) {
			uint32_t i = this->coordsToSite(i1, i2);
			uint32_t j = spaceNeigh(XPLUS, i);
			uint32_t k = spaceNeigh(YPLUS, i);
			uint32_t l = spaceNeigh(XPLUS, k);
			//change cols i,j,k,l of result
			const arma::Col<cpx>& ci = result.col(i);
			const arma::Col<cpx>& cj = result.col(j);
			const arma::Col<cpx>& ck = result.col(k);
			const arma::Col<cpx>& cl = result.col(l);
			num b_sh_hor = sh_hor;
			num b_sh_ver = sh_ver;
			if ((bc == APBC_X or bc == APBC_XY) and i1 == L-1) {
				//this plaquette has horizontal boundary crossing bonds and APBC
				b_sh_hor *= -1;
			}
			if ((bc == APBC_Y or bc == APBC_XY) and i2 == L-1) {
				//this plaquette has vertical boundary crossing bonds and APBC
				b_sh_ver *= -1;
			}
			new_col_i     = ch_hor*ch_ver*ci + ch_ver*b_sh_hor*cj + ch_hor*b_sh_ver*ck + b_sh_hor*b_sh_ver*cl;
			new_col_j     = ch_ver*b_sh_hor*ci + ch_hor*ch_ver*cj + b_sh_hor*b_sh_ver*ck + ch_hor*b_sh_ver*cl;
			new_col_k     = ch_hor*b_sh_ver*ci + b_sh_hor*b_sh_ver*cj + ch_hor*ch_ver*ck + ch_ver*b_sh_hor*cl;
			result.col(l) = b_sh_hor*b_sh_ver*ci + ch_hor*b_sh_ver*cj + ch_ver*b_sh_hor*ck + ch_hor*ch_ver*cl;
			result.col(i) = new_col_i;
			result.col(j) = new_col_j;
			result.col(k) = new_col_k;
		}
	}
}

template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
											 const Matrix& A, Band band, int sign, bool invertedCbOrder) {
    MatCpx result = A;      //can't avoid this copy

    //order reversed wrt cbLMultHoppingExp
    if (not invertedCbOrder) {
    	cb_assaad_applyBondFactorsRight(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    	cb_assaad_applyBondFactorsRight(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    } else {
    	cb_assaad_applyBondFactorsRight(result, 1, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    	cb_assaad_applyBondFactorsRight(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
    }
    return result;
}

template<bool TD, CheckerboardMethod CB>
template<class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
											 const Matrix& A, Band band, int sign, bool) {
    MatCpx result = A;      //can't avoid this copy

    //order of matrix multiplications symmetric
	//perform the multiplication e^(+-dtau K_1/2) e^(+-dtau K_0) e^(+-dtau K_a/2) X
	cb_assaad_applyBondFactorsRight(result, 1, coshHopHorHalf[band], sign * sinhHopHorHalf[band],
			                 	 	coshHopVerHalf[band], sign * sinhHopVerHalf[band]);
	cb_assaad_applyBondFactorsRight(result, 0, coshHopHor[band], sign * sinhHopHor[band], coshHopVer[band], sign * sinhHopVer[band]);
	cb_assaad_applyBondFactorsRight(result, 1, coshHopHorHalf[band], sign * sinhHopHorHalf[band],
									coshHopVerHalf[band], sign * sinhHopVerHalf[band]);

	return result;
}

// with sign = +/- 1, band = XBAND|YBAND: return A * E^(sign * dtau * K_band)
template<bool TD, CheckerboardMethod CB>
template <class Matrix> inline
MatCpx DetSDW<TD,CB>::cbRMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder) {
	return cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB>(),
								  A, band, sign, invertedCbOrder);
}





template<bool TD, CheckerboardMethod CB> inline
MatCpx DetSDW<TD,CB>::leftMultiplyBk(const MatCpx& orig, uint32_t k) {
    //helper: submatrix block for a matrix
//  auto block = [this](const MatCpx& mat, uint32_t row, uint32_t col) {
//      return mat.submat( row * N, col * N,
//                        (row + 1) * N - 1, (col + 1) * N - 1);
//  };
#define block(mat,row,col) mat.submat( (row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)

    num muTerm = std::exp(dtau*mu);			// include chemical potential here

    const auto& kphi0 = phi0.col(k);
    const auto& kphi1 = phi1.col(k);
    const auto& kphi2 = phi2.col(k);
    const VecNum c = muTerm * phiCosh.col(k);   // cosh(dtau * |phi|)
    const auto& kphiSinh = phiSinh.col(k);  	// sinh(dtau * |phi|) / |phi|
    VecNum ax  =  muTerm * kphi2 % kphiSinh;
    VecNum max = -muTerm * kphi2 % kphiSinh;
    VecCpx b  {kphi0, -kphi1};
    VecCpx bc {kphi0, kphi1};
    VecCpx mbx  = muTerm * -b  % kphiSinh;
    VecCpx mbcx = muTerm * -bc % kphiSinh;

    MatCpx result(4*N, 4*N);

    for (uint32_t col = 0; col < 4; ++col) {
        using arma::diagmat;
        //only three terms each time because of zero blocks in the E^(-dtau*V) matrix
        block(result, 0, col) = diagmat(c)   * cbLMultHoppingExp(block(orig, 0, col), XBAND, -1, false)
                              + diagmat(max)  * cbLMultHoppingExp(block(orig, 2, col), YBAND, -1, false)
                              + diagmat(mbx)  * cbLMultHoppingExp(block(orig, 3, col), YBAND, -1, false);

        block(result, 1, col) = diagmat(c)   * cbLMultHoppingExp(block(orig, 1, col), XBAND, -1, false)
                              + diagmat(mbcx) * cbLMultHoppingExp(block(orig, 2, col), YBAND, -1, false)
                              + diagmat(ax) * cbLMultHoppingExp(block(orig, 3, col), YBAND, -1, false);

        block(result, 2, col) = diagmat(max)  * cbLMultHoppingExp(block(orig, 0, col), XBAND, -1, false)
                              + diagmat(mbx)  * cbLMultHoppingExp(block(orig, 1, col), XBAND, -1, false)
                              + diagmat(c)   * cbLMultHoppingExp(block(orig, 2, col), YBAND, -1, false);

        block(result, 3, col) = diagmat(mbcx) * cbLMultHoppingExp(block(orig, 0, col), XBAND, -1, false)
                              + diagmat(ax) * cbLMultHoppingExp(block(orig, 1, col), XBAND, -1, false)
                              + diagmat(c)   * cbLMultHoppingExp(block(orig, 3, col), YBAND, -1, false);
    }
#undef block
    return result;
}



template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::checkerboardLeftMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1) {
    assert(k2 > k1);
    assert(k2 <= m);

    MatCpx result = leftMultiplyBk(A, k1 + 1);

    for (uint32_t k = k1 + 2; k <= k2; ++k) {
        result = leftMultiplyBk(result, k);
    }

    //chemical potential terms included by leftMultiplyBk
    //result *= std::exp(+dtau * (k2 - k1) * mu);

    return result;
}


template<bool TD, CheckerboardMethod CB> inline
MatCpx DetSDW<TD,CB>::leftMultiplyBkInv(const MatCpx& orig, uint32_t k) {
    //helper: submatrix block for a matrix
//  auto block = [this](const MatCpx& mat, uint32_t row, uint32_t col) {
//      return mat.submat( row * N, col * N,
//                        (row + 1) * N - 1, (col + 1) * N - 1);
//  };
#define block(mat,row,col) mat.submat( (row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)

	num muTerm = std::exp(-dtau*mu);			// include chemical potential here

    const auto& kphi0 = phi0.col(k);
    const auto& kphi1 = phi1.col(k);
    const auto& kphi2 = phi2.col(k);
    const VecNum c = muTerm * phiCosh.col(k);     // cosh(dtau * |phi|)
    const auto& kphiSinh = phiSinh.col(k);  	  // sinh(dtau * |phi|) / |phi|
    VecNum ax  =  muTerm * kphi2 % kphiSinh;
    VecNum max = -muTerm * kphi2 % kphiSinh;
    VecCpx b  {kphi0, -kphi1};
    VecCpx bc {kphi0, +kphi1};
    VecCpx bx  = muTerm * b  % kphiSinh;
    VecCpx bcx = muTerm * bc % kphiSinh;

    MatCpx result(4*N, 4*N);

    for (uint32_t col = 0; col < 4; ++col) {
        using arma::diagmat;
        //only three terms each time because of zero blocks in the E^(dtau*V) matrix
        block(result, 0, col) = cbLMultHoppingExp(diagmat(c)   * block(orig, 0, col), XBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(ax)  * block(orig, 2, col), XBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(bx)  * block(orig, 3, col), XBAND, +1, true);

        block(result, 1, col) = cbLMultHoppingExp(diagmat(c)   * block(orig, 1, col), XBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(bcx) * block(orig, 2, col), XBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(max) * block(orig, 3, col), XBAND, +1, true);

        block(result, 2, col) = cbLMultHoppingExp(diagmat(ax)  * block(orig, 0, col), YBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(bx)  * block(orig, 1, col), YBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(c)   * block(orig, 2, col), YBAND, +1, true);

        block(result, 3, col) = cbLMultHoppingExp(diagmat(bcx) * block(orig, 0, col), YBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(max) * block(orig, 1, col), YBAND, +1, true)
                              + cbLMultHoppingExp(diagmat(c)   * block(orig, 3, col), YBAND, +1, true);
    }
#undef block
    return result;
}


template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::checkerboardLeftMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1) {
    assert(k2 > k1);
    assert(k2 <= m);

//    MatCpx result = leftMultiplyBkInv(A, k1 + 1);
//
//    for (uint32_t k = k1 + 2; k <= k2; ++k) {
//        result = leftMultiplyBkInv(result, k);
//    }

    MatCpx result = leftMultiplyBkInv(A, k2);

    for (uint32_t k = k2 - 1; k >= k1 + 1; --k) {
    	result = leftMultiplyBkInv(result, k);
    }

    //chemical potential terms already included
    //result *= std::exp(-dtau * (k2 - k1) * mu);

    return result;
}

template<bool TD, CheckerboardMethod CB> inline
MatCpx DetSDW<TD,CB>::rightMultiplyBk(const MatCpx& orig, uint32_t k) {
    //helper: submatrix block for a matrix
//  auto block = [this](const MatCpx& mat, uint32_t row, uint32_t col) {
//      return mat.submat( row * N, col * N,
//                        (row + 1) * N - 1, (col + 1) * N - 1);
//  };
#define block(mat,row,col) mat.submat( (row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)

	num muTerm = std::exp(dtau*mu);			// include chemical potential here

    const auto& kphi0 = phi0.col(k);
    const auto& kphi1 = phi1.col(k);
    const auto& kphi2 = phi2.col(k);
    VecNum c = muTerm * phiCosh.col(k);         // cosh(dtau * |phi|)
    const auto& kphiSinh = phiSinh.col(k);      // sinh(dtau * |phi|) / |phi|
    VecNum ax  =  muTerm * kphi2 % kphiSinh;
    VecNum max = -muTerm * kphi2 % kphiSinh;
    VecCpx b  {kphi0, -kphi1};
    VecCpx bc {kphi0, kphi1};
    VecCpx mbx  = muTerm * -b  % kphiSinh;
    VecCpx mbcx = muTerm * -bc % kphiSinh;

    MatCpx result(4*N, 4*N);

    for (uint32_t row = 0; row < 4; ++row) {
            using arma::diagmat;
            //only three terms each time because of zero blocks in the E^(-dtau*V) matrix
            block(result, row, 0) = cbRMultHoppingExp(block(orig, row, 0) * diagmat(c),    XBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 2) * diagmat(max),  XBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 3) * diagmat(mbcx), XBAND, -1, false);

            block(result, row, 1) = cbRMultHoppingExp(block(orig, row, 1) * diagmat(c),    XBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 2) * diagmat(mbx),  XBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 3) * diagmat(ax),   XBAND, -1, false);

            block(result, row, 2) = cbRMultHoppingExp(block(orig, row, 0) * diagmat(max),  YBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 1) * diagmat(mbcx), YBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 2) * diagmat(c),    YBAND, -1, false);

            block(result, row, 3) = cbRMultHoppingExp(block(orig, row, 0) * diagmat(mbx),  YBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 1) * diagmat(ax),   YBAND, -1, false)
                                  + cbRMultHoppingExp(block(orig, row, 3) * diagmat(c),    YBAND, -1, false);
    }

#undef block
    return result;
}

template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::checkerboardRightMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1) {
    assert(k2 > k1);
    assert(k2 <= m);

    MatCpx result = rightMultiplyBk(A, k2);

    for (uint32_t k = k2 - 1; k >= k1 +1; --k) {
        result = rightMultiplyBk(result, k);
    }

    //chemical potential terms included above
    //result *= std::exp(+dtau * (k2 - k1) * mu);

    return result;
}

template<bool TD, CheckerboardMethod CB> inline
MatCpx DetSDW<TD,CB>::rightMultiplyBkInv(const MatCpx& orig, uint32_t k) {
    //helper: submatrix block for a matrix
//  auto block = [this](const MatCpx& mat, uint32_t row, uint32_t col) {
//      return mat.submat( row * N, col * N,
//                        (row + 1) * N - 1, (col + 1) * N - 1);
//  };
#define block(mat,row,col) mat.submat( (row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)

	num muTerm = std::exp(-dtau*mu);			// include chemical potential here

    const auto& kphi0 = phi0.col(k);
    const auto& kphi1 = phi1.col(k);
    const auto& kphi2 = phi2.col(k);
    const VecNum c = muTerm * phiCosh.col(k);         // cosh(dtau * |phi|)
    const auto& kphiSinh = phiSinh.col(k);  		  // sinh(dtau * |phi|) / |phi|
    VecNum ax  =  muTerm * kphi2 % kphiSinh;
    VecNum max = -muTerm * kphi2 % kphiSinh;
    VecCpx b  {kphi0, -kphi1};
    VecCpx bc {kphi0, kphi1};
    VecCpx bx  = muTerm * b  % kphiSinh;
    VecCpx bcx = muTerm * bc % kphiSinh;

    MatCpx result(4*N, 4*N);

    for (uint32_t row = 0; row < 4; ++row) {
        using arma::diagmat;
        //only three terms each time because of zero blocks in the E^(+dtau*V) matrix
        block(result, row, 0) = cbRMultHoppingExp(block(orig, row, 0), XBAND, +1, true) * diagmat(c)
                              + cbRMultHoppingExp(block(orig, row, 2), YBAND, +1, true) * diagmat(ax)
                              + cbRMultHoppingExp(block(orig, row, 3), YBAND, +1, true) * diagmat(bcx);

        block(result, row, 1) = cbRMultHoppingExp(block(orig, row, 1), XBAND, +1, true) * diagmat(c)
                              + cbRMultHoppingExp(block(orig, row, 2), YBAND, +1, true) * diagmat(bx)
                              + cbRMultHoppingExp(block(orig, row, 3), YBAND, +1, true) * diagmat(max);

        block(result, row, 2) = cbRMultHoppingExp(block(orig, row, 0), XBAND, +1, true) * diagmat(ax)
                              + cbRMultHoppingExp(block(orig, row, 1), XBAND, +1, true) * diagmat(bcx)
                              + cbRMultHoppingExp(block(orig, row, 2), YBAND, +1, true) * diagmat(c);

        block(result, row, 3) = cbRMultHoppingExp(block(orig, row, 0), XBAND, +1, true) * diagmat(bx)
                              + cbRMultHoppingExp(block(orig, row, 1), XBAND, +1, true) * diagmat(max)
                              + cbRMultHoppingExp(block(orig, row, 3), YBAND, +1, true) * diagmat(c);
    }
    return result;
#undef block
}

template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::checkerboardRightMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1) {
    assert(k2 > k1);
    assert(k2 <= m);

//    MatCpx result = rightMultiplyBkInv(A, k2);
//
//    for (uint32_t k = k2 - 1; k >= k1 +1; --k) {
//        result = rightMultiplyBkInv(result, k);
//    }

    MatCpx result = rightMultiplyBkInv(A, k1 + 1);

    for (uint32_t k = k1 + 2; k <= k2; ++k) {
    	result = rightMultiplyBkInv(result, k);
    }

    //chemical potential terms included above
    //result *= std::exp(-dtau * (k2 - k1) * mu);

    return result;
}






template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::updateInSlice(uint32_t timeslice) {
    timing.start("sdw-updateInSlice");

    for (uint32_t rep = 0; rep < repeatUpdateInSlice; ++rep) {
    	switch(updateMethod) {
    	case ITERATIVE:
    		updateInSlice_iterative(timeslice);
    		break;
    	case WOODBURY:
    		updateInSlice_woodbury(timeslice);
    		break;
    	case DELAYED:
    		updateInSlice_delayed(timeslice);
    		break;
    	}
	}

    if (rescale) {
    	if (performedSweeps % rescaleInterval == 0) {
    		num rnd = rng.rand01();
    		if (rnd <= 0.5) {
    			attemptGlobalRescaleMove(timeslice, rescaleGrowthFactor);
    		} else {
    			attemptGlobalRescaleMove(timeslice, rescaleShrinkFactor);
    		}
    	}
    }

    timing.stop("sdw-updateInSlice");
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::updateInSlice_iterative(uint32_t timeslice) {
    lastAccRatioLocal = 0;
    for (uint32_t site = 0; site < N; ++site) {
        Phi newphi = proposeNewField(site, timeslice);

//      VecNum oldphi0 = phi0.col(timeslice);
//      VecNum oldphi1 = phi1.col(timeslice);
//      VecNum oldphi2 = phi2.col(timeslice);
//      debugSaveMatrix(oldphi0, "old_phi0");
//      debugSaveMatrix(oldphi1, "old_phi1");
//      debugSaveMatrix(oldphi2, "old_phi2");

//      VecNum newphi0 = phi0.col(timeslice);
//      VecNum newphi1 = phi1.col(timeslice);
//      VecNum newphi2 = phi2.col(timeslice);
//      newphi0[site] = newphi[0];
//      newphi1[site] = newphi[1];
//      newphi2[site] = newphi[2];
//      debugSaveMatrix(newphi0, "new_phi0");
//      debugSaveMatrix(newphi1, "new_phi1");
//      debugSaveMatrix(newphi2, "new_phi2");

        num dsphi = deltaSPhi(site, timeslice, newphi);
        num probSPhi = std::exp(-dsphi);
//      std::cout << probSPhi << std::endl;

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
        checkarray<VecCpx, 4> rows;
        for (uint32_t r = 0; r < 4; ++r) {
        	//TODO: Here are some unnecessary operations: deltanonzero contains many repeated
        	//elements, and even some zeros
            rows[r] = VecCpx(4*N);
            for (uint32_t col = 0; col < 4*N; ++col) {
                rows[r][col] = -deltanonzero(r,0) * g.col(col)[site];
            }
            rows[r][site] += deltanonzero(r,0);
            for (uint32_t dc = 1; dc < 4; ++dc) {
                for (uint32_t col = 0; col < 4*N; ++col) {
                    rows[r][col] += -deltanonzero(r,dc) * g.col(col)[site + dc*N];
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
        for (uint32_t l = 0; l < 4; ++l) {
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
//      checkarray<VecCpx, 4> invRows = {{rows[0], rows[1], rows[2], rows[3]}};
        //****END-DEBUG

        //****
        //DEBUG: This slow code was working before.
        //Compare results with the better-performing sherman-morrison code
//      SpMatCpx delta(4*N, 4*N);
//      arma::uvec::fixed<4> idx = {site, site + N, site + 2*N, site + 3*N};
//      //Armadilo lacks non-contiguous submatrix views for sparse matrices
//      uint32_t i = 0;
//      for (auto col: idx) {
//          uint32_t j = 0;
//          for (auto row: idx) {
//              delta(row, col) = deltanonzero(j, i);
//              ++j;
//          }
//          ++i;
//      }


//      MatCpx deltaDense(4*N,4*N);
//      deltaDense = delta;
//      debugSaveMatrix(MatNum(arma::real(deltaDense)), "delta_real");
//      debugSaveMatrix(MatNum(arma::imag(deltaDense)), "delta_imag");

//      //inefficient!
//      static MatCpx eyeCpx = MatCpx(arma::eye(4*N, 4*N), arma::zeros(4*N, 4*N));
//      MatCpx target = eyeCpx + delta * (eyeCpx - g.slice(timeslice));

//      debugSaveMatrix(MatNum(arma::real(target)), "target_real");
//      debugSaveMatrix(MatNum(arma::imag(target)), "target_imag");

//      cpx weightRatio = arma::det(target);
//      std::cout << weightRatio << std::endl;

//      std::cout << weightRatio << " vs. " << det << std::endl;
        //END DEBUG
        //****

        //****
        // DEBUG-CHECK Delta*(I - G) vs. rows
//      MatCpx check = delta * (eyeCpx - g.slice(timeslice));
//      for (uint32_t r = 0; r < 4; ++r) {
//          check.row(site + r*N) -= rows[r].st();
//      }
//      std::cout << "Row check: " << arma::max(arma::max(arma::abs(check))) << std::endl;
//      debugSaveMatrix(MatNum(arma::real(check)), "check_real");
//      debugSaveMatrix(MatNum(arma::imag(check)), "check_imag");
//      exit(0);
        // END-DEBUG-CHECK
        //****

        //****
        //DEBUG-CHECK [I+Delta*(I - G)]^(-1) vs. invRows
//      MatCpx inv = arma::inv(target);
//      inv -= eyeCpx;
//      for (uint32_t r = 0; r < 4; ++r) {
//          inv.row(site + r*N) -= invRows[r].st();
//          inv.row(site + r*N)[site + r*N] += 1;
//      }
//      std::cout << "inv-div: " << arma::max(arma::max(inv)) << std::endl;
        //END-DEBUG-CHECK
        //****

        num probSFermion = det.real();

        //DEBUG: determinant computation from new routine updateInSlice_woodbury:
//        MatCpx::fixed<4,4> g_sub;
//        for (uint32_t a = 0; a < 4; ++a) {
//        	for (uint32_t b = 0; b < 4; ++b) {
//        		g_sub(a,b) = g(site + a*N, site + b*N);
//        	}
//        }
//        MatCpx::fixed<4,4> M = eye4cpx + (eye4cpx - g_sub) * deltanonzero;						//!
//        std::cout << "det: " << (probSFermion - arma::det(M).real()) / probSFermion << "\n";
        //END-DEBUG: relative difference: 0 or at most ~E-16 --> results are equal


        num prob = probSPhi * probSFermion;

        if (prob > 1.0 or rng.rand01() < prob) {
            //count accepted update
            lastAccRatioLocal += 1.0;

//          num phisBefore = phiAction();
            phi0(site, timeslice) = newphi[0];
            phi1(site, timeslice) = newphi[1];
            phi2(site, timeslice) = newphi[2];
            phiCosh(site, timeslice) = coshnewphi;
            phiSinh(site, timeslice) = sinhnewphi;
//          num phisAfter = phiAction();
//          std::cout << std::scientific << dsphi << " vs. " << phisAfter << " - " << phisBefore << " = " <<
//                  (phisAfter - phisBefore) << std::endl;

//          debugSaveMatrix(MatNum(arma::real(g.slice(timeslice))), "gslice_old_real");
//          debugSaveMatrix(MatNum(arma::imag(g.slice(timeslice))), "gslice_old_imag");
//          g.slice(timeslice) *= arma::inv(target);
//          debugSaveMatrix(MatNum(arma::real(g.slice(timeslice))), "gslice_new_real");
//          debugSaveMatrix(MatNum(arma::imag(g.slice(timeslice))), "gslice_new_imag");

            //****
            //DEBUG
//          MatCpx gPrimeRef = g.slice(timeslice) * arma::inv(target);
            //END DEBUG
            //****


            //DEBUG
            //compare with a full-force evaluation of G*[I + Delta*(I - G)]^{-1}
//            MatCpx delta(4*N, 4*N);
//            delta.fill(cpx(0,0));
//            arma::uvec::fixed<4> idx = {site, site + N, site + 2*N, site + 3*N};
//            uint32_t i = 0;
//            for (auto col: idx) {
//            	uint32_t j = 0;
//            	for (auto row: idx) {
//            		delta(row, col) = deltanonzero(j, i);
//            		++j;
//            	}
//            	++i;
//            }
//            MatCpx g_new_ref = g * arma::inv(
//            		arma::eye(4*N,4*N) + delta*(arma::eye(4*N,4*N) - g));
            //END DEBUG


            //DEBUG
            //Compare green's function updated with the method of updateInSlice_woodbury
            //with this one:
//            MatCpx g_woodbury = g;		//copy
//            MatCpx::fixed<4,4> g_woodbury_sub;
//            for (uint32_t a = 0; a < 4; ++a) {
//            	for (uint32_t b = 0; b < 4; ++b) {
//            		g_woodbury_sub(a,b) = g_woodbury(site + a*N, site + b*N);
//            	}
//            }
//            MatCpx::fixed<4,4> M = eye4cpx + (eye4cpx - g_woodbury_sub) * deltanonzero;	//!!
//            MatCpx mat_V(4, 4*N);
//            for (uint32_t r = 0; r < 4; ++r) {
//            	mat_V.row(r) = g_woodbury.row(site + r*N);
//            	mat_V(r, site + r*N) -= 1.0;
//            }
//            MatCpx g_woodbury_times_mat_U(4*N, 4);
//            for (uint32_t c = 0; c < 4; ++c) {
//            	g_woodbury_times_mat_U.col(c) = g_woodbury.col(site + c*N); //!!
//            }
//            g_woodbury_times_mat_U = g_woodbury_times_mat_U * deltanonzero;
//            g_woodbury += (g_woodbury_times_mat_U) * (arma::inv(M) * mat_V);
            //END DEBUG

            //DEBUG
            //Compare green's function updated with a large matrix woodbury formula
            //with this one:
//            MatCpx delta(4*N, 4*N);
//            delta.fill(cpx(0,0));
//            arma::uvec::fixed<4> idx = {site, site + N, site + 2*N, site + 3*N};
//            uint32_t i = 0;
//            for (auto col: idx) {
//            	uint32_t j = 0;
//            	for (auto row: idx) {
//            		delta(row, col) = deltanonzero(j, i);
//            		++j;
//            	}
//            	++i;
//            }
//            MatCpx mat_U_large = -delta;
//            MatCpx mat_V_large = arma::eye(4*N,4*N) - g;
//            MatCpx g_woodbury_large = g * (arma::eye(4*N,4*N) +
//            		mat_U_large *
//            		arma::inv(arma::eye(4*N,4*N) - mat_V_large * mat_U_large) *
//            		mat_V_large);
            //END DEBUG


            //DEBUG
            //Compare green's function updated with a reasonably sized matrix woodbury formula
            //with this one:
//            MatCpx mat_U_reas(4*N,4);
//            mat_U_reas.fill(cpx(0,0));
//            for (uint32_t k = 0; k < 4; ++k) {
//            	for (uint32_t l = 0; l < 4; ++l) {
//            		mat_U_reas(site + k*N, l) = -deltanonzero(k, l);
//            	}
//            }
//            MatCpx mat_V_reas(4,4*N);
//            for (uint32_t l = 0; l < 4; ++l) {
//            	mat_V_reas.row(l) = arma::conv_to<MatCpx>::from((arma::eye(4*N,4*N) - g)).row(site + l*N);
//            }
//            MatCpx::fixed<4,4> g_woodbury_reas_sub;
//            for (uint32_t a = 0; a < 4; ++a) {
//            	for (uint32_t b = 0; b < 4; ++b) {
//            		g_woodbury_reas_sub(a,b) = g(site + a*N, site + b*N);
//            	}
//            }
//            MatCpx::fixed<4,4> M_rev = eye4cpx + (eye4cpx - g_woodbury_reas_sub) * deltanonzero;
//            //MatCpx g_woodbury_reas = g * (arma::eye(4*N,4*N) +
//            //		mat_U_reas *
//            //		arma::inv(eye4cpx - mat_V_reas * mat_U_reas) *
//            //		mat_V_reas);
//            MatCpx g_woodbury_reas = g * (arma::eye(4*N,4*N) +
//            		mat_U_reas *
//            		arma::inv(M_rev) *
//            		mat_V_reas);
            //END DEBUG


            //compensate for already included diagonal entries of I in invRows
            rows[0][site] -= 1;
            rows[1][site + N] -= 1;
            rows[2][site + 2*N] -= 1;
            rows[3][site + 3*N] -= 1;
            //compute G' = G * [I + Delta*(I - G)]^(-1) = G * [I + invRows]
            // [O(N^2)]
            MatCpx gTimesInvRows(4*N, 4*N);
            const auto& G = g;
            for (uint32_t col = 0; col < 4*N; ++col) {
                for (uint32_t row = 0; row < 4*N; ++row) {
                    gTimesInvRows(row, col) = G(row, site) * rows[0][col]
                                            + G(row, site + N) * rows[1][col]
                                            + G(row, site + 2*N) * rows[2][col]
                                            + G(row, site + 3*N) * rows[3][col]
                                            ;
                }
            }
            g += gTimesInvRows;

            //DEBUG
//            std::cout << "ref g_mean: " << arma::mean(arma::mean(arma::abs(((g - g_new_ref) / g)))) << "\n";
//            std::cout << "ref g_max: " << arma::max(arma::max(arma::abs(((g - g_new_ref) / g)))) << "\n";
            //compare with a full-force evaluation of G*[I + Delta*(I - G)]^{-1}

            //DEBUG
            //Compare green's function updated with an large woodbury formula
            //with this one:
//            std::cout << "woodbury large g_mean: " << arma::mean(arma::mean(arma::abs(((g - g_woodbury_large) / g)))) << "\n";
//            std::cout << "woodbury large g_max: " << arma::max(arma::max(arma::abs(((g - g_woodbury_large) / g)))) << "\n";
            //END DEBUG

            //DEBUG
            //Compare green's function updated with a reasonably sized matrix woodbury formula
            //with this one:
//            std::cout << "woodbury reas g_mean: " << arma::mean(arma::mean(arma::abs(((g - g_woodbury_reas) / g)))) << "\n";
//            std::cout << "woodbury reas g_max: " << arma::max(arma::max(arma::abs(((g - g_woodbury_reas) / g)))) << "\n";
            //END DEBUG


            //DEBUG
            //Compare green's function updated with the method of updateInSlice_woodbury
            //with this one:
//            std::cout << "woodbury g_mean: " << arma::mean(arma::mean(arma::abs(((g - g_woodbury) / g)))) << "\n";
//            std::cout << "woodbury g_max: " << arma::max(arma::max(arma::abs(((g - g_woodbury) / g)))) << "\n";
            //END DEBUG

            //****
            //DEBUG
//          std::cout << arma::max(arma::max(
//                  (arma::abs(gPrimeRef - g.slice(timeslice))))) << std::endl;
            //END DEBUG
            //****
        }
    }
    lastAccRatioLocal /= num(N);
}


template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::updateInSlice_woodbury(uint32_t timeslice) {
    lastAccRatioLocal = 0;
    for (uint32_t site = 0; site < N; ++site) {
        Phi newphi = proposeNewField(site, timeslice);
        num dsphi = deltaSPhi(site, timeslice, newphi);
        num probSPhi = std::exp(-dsphi);
        //delta = e^(-dtau*V_new)*e^(+dtau*V_old) - 1

        //TODO: put calculation of deltanonzero into separate function

        //compute non-zero elements of delta
        // deltanonzero is \Delta^i from the notes
        //
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

        //Compute the 4x4 submatrix of G that corresponds to the site i
        //g_sub = g[i::N, i::N]
        MatCpx::fixed<4,4> g_sub;
        for (uint32_t a = 0; a < 4; ++a) {
        	for (uint32_t b = 0; b < 4; ++b) {
        		g_sub(a,b) = g(site + a*N, site + b*N);
        	}
        }

        //the determinant ratio for the spin update is given by the determinant
        //of the following matrix M
        MatCpx::fixed<4,4> M = eye4cpx + (eye4cpx - g_sub) * deltanonzero;

        num probSFermion = arma::det(M).real();

        num prob = probSPhi * probSFermion;

        if (prob > 1.0 or rng.rand01() < prob) {
            //count accepted update
            lastAccRatioLocal += 1.0;

            phi0(site, timeslice) = newphi[0];
            phi1(site, timeslice) = newphi[1];
            phi2(site, timeslice) = newphi[2];
            phiCosh(site, timeslice) = coshnewphi;
            phiSinh(site, timeslice) = sinhnewphi;

            //update g

            MatCpx mat_V(4, 4*N);
            for (uint32_t r = 0; r < 4; ++r) {
            	mat_V.row(r) = g.row(site + r*N);
            	mat_V(r, site + r*N) -= 1.0;
            }

            //TODO: is it a good idea to do this copy? or would it be better to
            //compute the product directly with a non-contiguous subview?
            MatCpx g_times_mat_U(4*N, 4);
            for (uint32_t c = 0; c < 4; ++c) {
            	g_times_mat_U.col(c) = g.col(site + c*N);
            }
            g_times_mat_U = g_times_mat_U * deltanonzero;

            g += (g_times_mat_U) * (arma::inv(M) * mat_V);
        }
    }
    lastAccRatioLocal /= num(N);
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::updateInSlice_delayed(uint32_t timeslice) {
	lastAccRatioLocal = 0;

	auto getX = [this](uint32_t step) {
		return dud.X.cols(4*step, 4*step + 3);
	};
	auto getY = [this](uint32_t step) {
		return dud.Y.rows(4*step, 4*step + 3);
	};

	auto take4rows = [this](MatCpx& target, const MatCpx& source, uint32_t for_site) {
		for (uint32_t r = 0; r < 4; ++r) {
			target.row(r) = source.row(for_site + r*N);
		}
	};
	auto take4cols = [this](MatCpx& target, const MatCpx& source, uint32_t for_site) {
		for (uint32_t c = 0; c < 4; ++c) {
			target.col(c) = source.col(for_site + c*N);
		}
	};

	uint32_t site = 0;
	while (site < N) {
		uint32_t delayStepsNow = std::min(delaySteps, N - site);
		dud.X.set_size(4*N, 4*delayStepsNow);
		dud.Y.set_size(4*delayStepsNow, 4*N);
		uint32_t j = 0;
		while (j < delayStepsNow and site < N) {
			Phi newphi = proposeNewField(site, timeslice);
			num dsphi = deltaSPhi(site, timeslice, newphi);
			num probSPhi = std::exp(-dsphi);

			MatCpx::fixed<4,4> deltanonzero = get_deltanonzero(newphi, timeslice, site);
			//TODO: die naechsten drei Rechnungen werden auch schon in get_deltanonzero
			//durchgefuehrt...
			num normnewphi = arma::norm(newphi,2);
			num coshnewphi = std::cosh(dtau * normnewphi);
			num sinhnewphi = std::sinh(dtau * normnewphi) / normnewphi;

			take4rows(dud.Rj, g, site);
			for (uint32_t l = 0; l < j; ++l) {
				take4rows(dud.tempBlock, getX(l), site);
				dud.Rj += dud.tempBlock * getY(l);
			}

			take4cols(dud.Sj, dud.Rj, site);

			dud.Mj = eye4cpx - dud.Sj * deltanonzero + deltanonzero;
			num probSFermion = arma::det(dud.Mj).real();

			num prob = probSPhi * probSFermion;
			if (prob > 1.0 or rng.rand01() < prob) {
				//count accepted update
				lastAccRatioLocal += 1.0;

				phi0(site, timeslice) = newphi[0];
				phi1(site, timeslice) = newphi[1];
				phi2(site, timeslice) = newphi[2];
				phiCosh(site, timeslice) = coshnewphi;
				phiSinh(site, timeslice) = sinhnewphi;

				//we need Cj only to update X
				take4cols(dud.Cj, g, site);
				for (uint32_t l = 0; l < j; ++l) {
					take4cols(dud.tempBlock, getY(l), site);
					dud.Cj += getX(l) * dud.tempBlock;
				}
				//Rj is now Rj - \Id_j, for updating Y
				for (uint32_t rc = 0; rc < 4; ++rc) {
					uint32_t entry = site + rc * N;
					dud.Rj(rc, entry) -= cpx(1.0, 0.0);
				}

				//update X and Y
				getX(j) = dud.Cj * deltanonzero;
				getY(j) = arma::inv(dud.Mj) * dud.Rj;
				//count successful delayed update
				j += 1;
			}
			++site;
		}
		if (j > 0) {
			if (j < delayStepsNow) {
				dud.X.resize(4*N, 4*j);
				dud.Y.resize(4*j, 4*N);
			}
			//carry out the delayed updates of the Green's function
			g += dud.X*dud.Y;
		}
	}

	lastAccRatioLocal /= num(N);
}

template<bool TD, CheckerboardMethod CB>
MatCpx::fixed<4,4> DetSDW<TD,CB>::get_deltanonzero(Phi newphi, uint32_t timeslice, uint32_t site) {
	//delta = e^(-dtau*V_new)*e^(+dtau*V_old) - 1

	//compute non-zero elements of delta
	// deltanonzero is \Delta^i from the notes
	//
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
	return deltanonzero;
}



template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::updateInSliceThermalization(uint32_t timeslice) {
    updateInSlice(timeslice);

    accRatioLocalRA.addValue(lastAccRatioLocal);
    if (accRatioLocalRA.getSamplesAdded() % AccRatioAdjustmentSamples == 0) {
        num avgAccRatio = accRatioLocalRA.get();
        if (avgAccRatio < targetAccRatioLocal) {
            phiDelta *= phiDeltaShrinkFactor;
        } else if (avgAccRatio > targetAccRatioLocal) {
            phiDelta *= phiDeltaGrowFactor;
        }
    }
}


template<bool TD, CheckerboardMethod CB>
inline void DetSDW<TD,CB>::attemptGlobalRescaleMove(uint32_t timeslice, num factor) {
	timing.start("sdw-attemptGlobalRescaleMove");

	//see hand-written notes and Ipython notebook sdw-rescale-move to understand these formulas

	//original fields
	//TODO: unnecessary copies
	const VecNum a  {phi2.col(timeslice)};
	const VecCpx b  {phi0.col(timeslice), -phi1.col(timeslice)};
	const VecCpx bc {phi0.col(timeslice), +phi1.col(timeslice)};
	const VecNum x  {phiSinh.col(timeslice)};
	const VecNum c  {phiCosh.col(timeslice)};

	//rescaled fields
	const VecNum rphi0 {factor * phi0.col(timeslice)};
	const VecNum rphi1 {factor * phi1.col(timeslice)};
	const VecNum rphi2 {factor * phi2.col(timeslice)};
	const VecNum ra  {rphi2};
	const VecCpx rb  {rphi0, -rphi1};
	const VecCpx rbc {rphi0, +rphi1};
	using arma::pow; using arma::sqrt; using arma::sinh; using arma::cosh;
	const VecNum rnorm { sqrt(pow(rphi0,2) + pow(rphi1,2) + pow(rphi2,2)) };
	const VecNum rx    { sinh(dtau * rnorm) / rnorm };
	const VecNum rc    { cosh(dtau * rnorm) };

	// 1) Calculate Delta = exp(-dtau V(a',b',c'))*exp(+dtau V(a,b,c)) - 1
	// Delta is setup by 4x4 blocks of size NxN, each being diagonal.
	// 4 blocks are zero, apart from that there are 5 different blocks:
	const VecNum delta_a  { rc % a % x - ra % rx % c };
	const VecNum delta_ma { -delta_a };
	const VecNum delta_c  { rc % c - ra % rx % a % x - rx % arma::real(rb % bc) % x - arma::ones<VecNum>(N) };	//Note: rb % bc will result in a purely real result
	// the block diagonals that are complex are stored with real and imaginary parts
	// separated:
	const VecNum delta_b_r  { rc % arma::real(b) % x - arma::real(rb) % rx % c };
	const VecNum delta_b_i  { rc % arma::imag(b) % x - arma::imag(rb) % rx % c };
	const VecNum delta_bc_r { delta_b_r };
	const VecNum delta_bc_i { -delta_b_i };

	// real part of matrix represented by 4x4 array of pointers to our vectors
	using std::array; using std::cref;
	array< array<const VecNum*, 4>, 4> delta_r;
	delta_r[0][0] = &delta_c;
	delta_r[0][1] = 0;
	delta_r[0][2] = &delta_a;
	delta_r[0][3] = &delta_b_r;
	delta_r[1][0] = 0;
	delta_r[1][1] = &delta_c;
	delta_r[1][2] = &delta_bc_r;
	delta_r[1][3] = &delta_ma;
	delta_r[2][0] = &delta_a;
	delta_r[2][1] = &delta_b_r;
	delta_r[2][2] = &delta_c;
	delta_r[2][3] = 0;
	delta_r[3][0] = &delta_bc_r;
	delta_r[3][1] = &delta_ma;
	delta_r[3][2] = 0;
	delta_r[3][3] = &delta_c;

	// imaginary part of matrix: only the antidiagonal blocks
	array< array<const VecNum*, 4>, 4> delta_i {{}};		//init with nulls
	delta_i[0][3] = &delta_b_i;
	delta_i[1][2] = &delta_bc_i;
	delta_i[2][1] = &delta_b_i;
	delta_i[3][0] = &delta_bc_i;

	// 2) Compute the matrix M = I + Delta * (I - G(timeslice))
	MatCpx oneMinusG { arma::eye(4*N,4*N) - g };
	MatCpx M { arma::eye(4*N,4*N), arma::zeros(4*N,4*N) };
#define block(matrix, row, col) matrix.submat((row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)
	//real parts
	for (uint32_t row = 0; row < 4; ++row) {
		for (uint32_t col = 0; col < 4; ++col) {
			//skip the zero blocks of delta:
			uint32_t skip_i;
			switch (row) {
			case 0: skip_i = 1; break;
			case 1: skip_i = 0; break;
			case 2: skip_i = 3; break;
			case 3: skip_i = 2; break;
			}
			uint32_t start_i = 0;
			for (uint32_t i = start_i; i < 4; ++i) {
				if (i == skip_i) continue;
				block(M,row,col) += arma::diagmat(*(delta_r[row][i])) *
						block(oneMinusG, i, col);
			}
		}
	}
	//imaginary parts
	for (uint32_t row = 0; row < 4; ++row) {
		VecCpx temp { arma::zeros<VecNum>(N), *(delta_i[row][3 - row]) };   //antidiagonal: col = 3 - row
		for (uint32_t i = 0; i < 4; ++i) {
			block(M,row,i) += arma::diagmat(temp) * block(oneMinusG, 3 - row, i);
		}
	}
#undef Mblock

	// 3) Compute probability of accepting the global rescale move
	num probFermion = arma::det(M).real();
	num probBoson = std::exp(-deltaSPhiGlobalRescale(timeslice, factor));
	num prob = probFermion * probBoson;

	if (prob > 1.0 or rng.rand01() < prob) {
		//count accepted update
		++acceptedRescales;

		phi0.col(timeslice) = rphi0;
		phi1.col(timeslice) = rphi1;
		phi2.col(timeslice) = rphi2;

		using arma::trans; using arma::solve;
		g = trans(solve(trans(M), trans(g)));
		//TODO: the three transpositions here bug me
	}

	++attemptedRescales;

	timing.stop("sdw-attemptGlobalRescaleMove");
}


template<bool TD, CheckerboardMethod CB>
num DetSDW<TD,CB>::deltaSPhiGlobalRescale(uint32_t timeslice, num factor) {
	using std::pow;
	num delta1 = 0;
	for (uint32_t site_i = 0; site_i < N; ++site_i) {
		uint32_t neighSites[] = {spaceNeigh(XPLUS, site_i), spaceNeigh(YPLUS, site_i)};   //for Icpc
		for (uint32_t site_j : neighSites) {
			delta1 += pow(phi0(site_i, timeslice) - phi0(site_j, timeslice), 2)
					+ pow(phi1(site_i, timeslice) - phi1(site_j, timeslice), 2)
					+ pow(phi2(site_i, timeslice) - phi2(site_j, timeslice), 2);
		}
	}

	num delta2 = 0;
	for (uint32_t site_i = 0; site_i < N; ++site_i) {
		delta2 += pow(phi0(site_i, timeslice), 2)
				+ pow(phi1(site_i, timeslice), 2)
				+ pow(phi2(site_i, timeslice), 2);
	}

	num delta3 = 0;
	for (uint32_t site_i = 0; site_i < N; ++site_i) {
		delta3 += pow(pow(phi0(site_i, timeslice), 2)
					+ pow(phi1(site_i, timeslice), 2)
				    + pow(phi2(site_i, timeslice), 2), 2);
	}

	num delta4 = 0;
	uint32_t timeslicePlus = timeNeigh(ChainDir::PLUS, timeslice);
	uint32_t timesliceMinus = timeNeigh(ChainDir::MINUS, timeslice);
	for (uint32_t site_i = 0; site_i < N; ++site_i) {
		delta4 += (pow(factor,2) - 1.0) * (
				      pow(phi0(site_i, timeslice), 2)
					+ pow(phi1(site_i, timeslice), 2)
					+ pow(phi2(site_i, timeslice), 2)
				);
		delta4 -= (factor - 1.0) * (
					  phi0(site_i, timeslice) * (phi0(site_i, timesliceMinus) + phi0(site_i, timeslicePlus))
					+ phi1(site_i, timeslice) * (phi1(site_i, timesliceMinus) + phi1(site_i, timeslicePlus))
					+ phi2(site_i, timeslice) * (phi2(site_i, timesliceMinus) + phi2(site_i, timeslicePlus))
				);
	}

	num delta = (dtau/2.0) * (pow(factor,2) - 1.0) * delta1
			  + (dtau*r/2.0) * (pow(factor,2) - 1.0) * delta2
			  + (dtau*u/4.0) * (pow(factor,4) - 1.0) * delta3
			  + (1.0/(c*dtau)) * delta4;

	return delta;
}



template<bool TD, CheckerboardMethod CB>
typename DetSDW<TD,CB>::Phi DetSDW<TD,CB>::proposeNewField(uint32_t site, uint32_t timeslice) {
    (void) site; (void) timeslice;
    //TODO: make this smarter!

    Phi phi;
    phi[0] = phi0(site, timeslice);
    phi[1] = phi1(site, timeslice);
    phi[2] = phi2(site, timeslice);


//  static int comp = 0;
//  num r = rng.randRange(-PhiDelta, +PhiDelta);
//  phi[comp] += r;
//  comp = (comp + 1) % 3;

    for (auto& phi_comp: phi) {
        num r = rng.randRange(-phiDelta, +phiDelta);
        phi_comp += r;
    }

    return phi;
}

template<bool TD, CheckerboardMethod CB>
num DetSDW<TD,CB>::deltaSPhi(uint32_t site, uint32_t timeslice, const Phi newphi) {
    //switched to asymmetric numerical derivative
    using arma::dot;

    Phi oldphi;
    oldphi[0] = phi0(site, timeslice);
    oldphi[1] = phi1(site, timeslice);
    oldphi[2] = phi2(site, timeslice);

    Phi phiDiff = newphi - oldphi;

    num oldphiSq = dot(oldphi, oldphi);
    num newphiSq = dot(newphi, newphi);
    num phiSqDiff = newphiSq - oldphiSq;

    num oldphiPow4 = oldphiSq * oldphiSq;
    num newphiPow4 = newphiSq * newphiSq;
    num phiPow4Diff = newphiPow4 - oldphiPow4;

    uint32_t kEarlier = timeNeigh(ChainDir::MINUS, timeslice);
    Phi phiEarlier;
    phiEarlier[0] = phi0(site, kEarlier);
    phiEarlier[1] = phi1(site, kEarlier);
    phiEarlier[2] = phi2(site, kEarlier);
    uint32_t kLater = timeNeigh(ChainDir::PLUS, timeslice);
    Phi phiLater;
    phiLater[0] = phi0(site, kLater);
    phiLater[1] = phi1(site, kLater);
    phiLater[2] = phi2(site, kLater);
    Phi phiTimeNeigh = phiLater + phiEarlier;

    Phi phiZero;
    phiZero[0] = phiZero[1] = phiZero[2] = 0;
    Phi phiSpaceNeigh = std::accumulate(
            spaceNeigh.beginNeighbors(site),
            spaceNeigh.endNeighbors(site),
            phiZero,
            [this, timeslice] (Phi accum, uint32_t neighSite) -> Phi{
                accum[0] += phi0(neighSite, timeslice);
                accum[1] += phi1(neighSite, timeslice);
                accum[2] += phi2(neighSite, timeslice);
                return accum;
            }
    );


    num delta1 = (1.0 / (c * c * dtau)) * (phiSqDiff - dot(phiTimeNeigh, phiDiff));

    num delta2 = 0.5 * dtau * (z * phiSqDiff - 2.0 * dot(phiSpaceNeigh, phiDiff));

    num delta3 = dtau * (0.5 * r * phiSqDiff + 0.25 * u * phiPow4Diff);

    return delta1 + delta2 + delta3;
}


template<bool TD, CheckerboardMethod CB>
num DetSDW<TD,CB>::phiAction() {
    //switched to asymmetric numerical derivative
    arma::field<Phi> phi(N, m+1);
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for (uint32_t site = 0; site < N; ++site) {
            phi(site, timeslice)[0] = phi0(site, timeslice);
            phi(site, timeslice)[1] = phi1(site, timeslice);
            phi(site, timeslice)[2] = phi2(site, timeslice);
        }
    }
    num action = 0;
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for (uint32_t site = 0; site < N; ++site) {
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


template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::thermalizationOver() {
    std::cout << "After thermalization: phiDelta = " << phiDelta << '\n'
              << "recent local accRatio = " << accRatioLocalRA.get()
              << std::endl;
    if (rescale) {
    	num ratio = num(acceptedRescales) / num(attemptedRescales);
    	std::cout << "Global rescale move acceptance ratio = " << ratio
    			  << std::endl;
    }
}



template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::sweepSimple(bool takeMeasurements) {
    sweepSimple_skeleton(takeMeasurements,
    					 sdwComputeBmat(this),
    					 [this]() {this->initMeasurements();},
    					 [this](uint32_t timeslice) {this->measure(timeslice);},
    					 [this]() {this->finishMeasurements();});
    ++performedSweeps;
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::sweepSimpleThermalization() {
    sweepSimpleThermalization_skeleton(sdwComputeBmat(this));
    ++performedSweeps;
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::sweep(bool takeMeasurements) {
    sweep_skeleton(takeMeasurements,
    			   sdwLeftMultiplyBmat(this), sdwRightMultiplyBmat(this),
                   sdwLeftMultiplyBmatInv(this), sdwRightMultiplyBmatInv(this),
                   [this]() {this->initMeasurements();},
                   [this](uint32_t timeslice) {this->measure(timeslice);},
                   [this]() {this->finishMeasurements();});
    ++performedSweeps;
}

template<bool TD, CheckerboardMethod CB>
void DetSDW<TD,CB>::sweepThermalization() {
    sweepThermalization_skeleton(sdwLeftMultiplyBmat(this), sdwRightMultiplyBmat(this),
                                 sdwLeftMultiplyBmatInv(this), sdwRightMultiplyBmatInv(this));
    ++performedSweeps;
}


template<bool TD, CheckerboardMethod CB>
MatCpx DetSDW<TD,CB>::shiftGreenSymmetric() {
	typedef arma::subview<cpx> SubMatCpx;			//don't do references or const-references of this type
	if (CB == CB_NONE) {
		//non-checkerboard
		return shiftGreenSymmetric_impl(
						//rightMultiply
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input * propK_half_inv[band];
			},
			//leftMultiply
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = propK_half[band] * input;
			}
		);
	}
	// unclear to me: why do I need to explicitly qualify 'this->' in the lambdas?
	else if (CB == CB_SANTOS) {
		return shiftGreenSymmetric_impl(
						//rightMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies e^{+ dtau K^band_b / 2} e^{+ dtau K^band_a / 2}
			// to the right of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;			//copy
				this->cb_santos_applyBondFactorsRight(output, YPLUS, 1, coshHopVerHalf[band], +sinhHopVerHalf[band]);
				this->cb_santos_applyBondFactorsRight(output, XPLUS, 1, coshHopHorHalf[band], +sinhHopHorHalf[band]);
				this->cb_santos_applyBondFactorsRight(output, YPLUS, 0, coshHopVerHalf[band], +sinhHopVerHalf[band]);
				this->cb_santos_applyBondFactorsRight(output, XPLUS, 0, coshHopHorHalf[band], +sinhHopHorHalf[band]);
			},
			//leftMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies e^{- dtau K^band_a / 2} e^{- dtau K^band_b / 2}
			// to the left of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;			//copy
				this->cb_santos_applyBondFactorsLeft(output, XPLUS, 0, coshHopHorHalf[band], -sinhHopHorHalf[band]);
				this->cb_santos_applyBondFactorsLeft(output, YPLUS, 0, coshHopVerHalf[band], -sinhHopVerHalf[band]);
				this->cb_santos_applyBondFactorsLeft(output, XPLUS, 1, coshHopHorHalf[band], -sinhHopHorHalf[band]);
				this->cb_santos_applyBondFactorsLeft(output, YPLUS, 1, coshHopVerHalf[band], -sinhHopVerHalf[band]);
			}
		);
	}
	else if (CB == CB_ASSAAD) {
		return shiftGreenSymmetric_impl(
						//rightMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies [Input] * e^{+ dtau K^band_b / 2} e^{+ dtau K^band_a / 2}
			// to the right of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;      //copy
				this->cb_assaad_applyBondFactorsRight(output, 1, coshHopHorHalf[band], +sinhHopHorHalf[band],
																 coshHopVerHalf[band], +sinhHopVerHalf[band]);
				this->cb_assaad_applyBondFactorsRight(output, 0, coshHopHorHalf[band], +sinhHopHorHalf[band],
														   	   	 coshHopVerHalf[band], +sinhHopVerHalf[band]);
			},
			//leftMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies e^{- dtau K^band_b / 2} e^{- dtau K^band_a / 2} * [Input]
			// to the left of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;      //copy
				this->cb_assaad_applyBondFactorsLeft(output, 0, coshHopHorHalf[band], -sinhHopHorHalf[band],
											 	 	 	  	  	coshHopVerHalf[band], -sinhHopVerHalf[band]);
				this->cb_assaad_applyBondFactorsLeft(output, 1, coshHopHorHalf[band], -sinhHopHorHalf[band],
																coshHopVerHalf[band], -sinhHopVerHalf[band]);
			}
		);
	}
	else if (CB == CB_ASSAAD_BERG) {
		return shiftGreenSymmetric_impl(
			//rightMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies [Input] * e^{+ dtau K^band_b / 2} e^{+ dtau K^band_a / 2}
			// to the right of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;      //copy
				this->cb_assaad_applyBondFactorsRight(output, 1, coshHopHorHalf[band], +sinhHopHorHalf[band],
																 coshHopVerHalf[band], +sinhHopVerHalf[band]);
				this->cb_assaad_applyBondFactorsRight(output, 0, coshHopHorHalf[band], +sinhHopHorHalf[band],
														   	   	 coshHopVerHalf[band], +sinhHopVerHalf[band]);
			},
			//leftMultiply
			// output and input are NxN blocks of a complex matrix
			// this effectively multiplies e^{- dtau K^band_a / 2} e^{- dtau K^band_b / 2} * [Input]
			// to the left of input and stores the result in output
			[this](SubMatCpx output, SubMatCpx input, Band band) -> void {
				output = input;      //copy
				this->cb_assaad_applyBondFactorsLeft(output, 1, coshHopHorHalf[band], -sinhHopHorHalf[band],
											 	 	 	  	  	coshHopVerHalf[band], -sinhHopVerHalf[band]);
				this->cb_assaad_applyBondFactorsLeft(output, 0, coshHopHorHalf[band], -sinhHopHorHalf[band],
																coshHopVerHalf[band], -sinhHopVerHalf[band]);
			}
		);
	}
}

//RightMultiply and LeftMultiply should be functors for complex matrix subviews,
//that take parameters (output, input, [BAND]).  Armadillo submatrix views apparently do not have
//any const correctness, and passing them by reference makes no sense (they are rich
//references in a sense)
template<bool TD, CheckerboardMethod CB>
template<class RightMultiply, class LeftMultiply>
MatCpx DetSDW<TD,CB>::shiftGreenSymmetric_impl(RightMultiply rightMultiply, LeftMultiply leftMultiply) {
    //submatrix view helper for a 4N*4N matrix
#define block(matrix, row, col) matrix.submat((row) * N, (col) * N, ((row) + 1) * N - 1, ((col) + 1) * N - 1)
	MatCpx tempG(4*N, 4*N);
	const MatCpx& oldG = g;
	//multiply e^(dtau/2 K) from the right
	for (uint32_t row = 0; row < 4; ++row) {
		//block(tempG, row, 0) = block(oldG, row, 0) * propKx_half_inv;
		rightMultiply(block(tempG, row, 0), block(oldG, row, 0), XBAND);
		rightMultiply(block(tempG, row, 1), block(oldG, row, 1), XBAND);
		rightMultiply(block(tempG, row, 2), block(oldG, row, 2), YBAND);
		rightMultiply(block(tempG, row, 3), block(oldG, row, 3), YBAND);
	}
	//multiply e^(-dtau/2 K) from the left
	MatCpx newG(4*N, 4*N);
	for (uint32_t col = 0; col < 4; ++col) {
		//block(newG, 0, col) = propKx_half * block(tempG, 0, col);
		leftMultiply(block(newG, 0, col), block(tempG, 0, col), XBAND);
		leftMultiply(block(newG, 1, col), block(tempG, 1, col), XBAND);
		leftMultiply(block(newG, 2, col), block(tempG, 2, col), YBAND);
		leftMultiply(block(newG, 3, col), block(tempG, 3, col), YBAND);
	}
#undef block
	return newG;
}



//explicit template instantiations:
template class DetSDW<true,CB_NONE>;
template class DetSDW<false,CB_NONE>;
template class DetSDW<true,CB_SANTOS>;
template class DetSDW<false,CB_SANTOS>;
template class DetSDW<true,CB_ASSAAD>;
template class DetSDW<false,CB_ASSAAD>;
template class DetSDW<true,CB_ASSAAD_BERG>;
template class DetSDW<false,CB_ASSAAD_BERG>;

