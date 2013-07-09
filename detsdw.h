/*
 * detsdw.h
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#ifndef DETSDW_H_
#define DETSDW_H_

#include <tuple>
#include <vector>
#include <complex>
#include <omp.h>
#include "rngwrapper.h"
#include "detmodel.h"
#include "neighbortable.h"
#include "RunningAverage.h"

typedef std::complex<num> cpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<cpx> SpMatCpx;
typedef arma::Col<cpx> VecCpx;
typedef arma::Cube<cpx> CubeCpx;

class SerializeContentsKey;

std::unique_ptr<DetModel> createDetSDW(RngWrapper& rng, ModelParams pars);

// template parameters: evaluate time-displaced Green functions? do a checker-board decomposition?
template <bool TimeDisplaced, bool CheckerBoard>
class DetSDW: public DetModelGC<1, cpx, TimeDisplaced> {
	DetSDW(RngWrapper& rng, const ModelParams& pars	);
public:
	friend std::unique_ptr<DetModel> createDetSDW(RngWrapper& rng, ModelParams pars);
	virtual ~DetSDW();
	virtual uint32_t getSystemN() const;
	virtual MetadataMap prepareModelMetadataMap() const;
    virtual void measure();

    virtual void thermalizationOver();

    virtual void sweep();
    virtual void sweepThermalization();
    virtual void sweepSimple();
    virtual void sweepSimpleThermalization();
protected:
    typedef DetModelGC<1, cpx, TimeDisplaced> Base;
    // stupid C++ weirdness forces us to explicitly "import" these protected base
    // class member variables:
    // (see: http://stackoverflow.com/questions/11405/gcc-problem-using-a-member-of-a-base-class-that-depends-on-a-template-argument )
    using Base::dtau;
    using Base::m;
    using Base::green;
    using Base::greenFwd;
    using Base::greenBwd;
    using Base::UdVStorage;
    using Base::lastSweepDir;
    using Base::obsScalar;
    using Base::obsVector;
    using Base::obsKeyValue;
    using Base::beta;
    using Base::s;

	RngWrapper& rng;

	static const uint32_t d = 2;
	static const uint32_t z = 2*d;
	const bool checkerboard;
	const uint32_t L;
	const uint32_t N;
	const num r;
	const num txhor;
	const num txver;
	const num tyhor;
	const num tyver;
	const num mu;
	const num c;
	const num u;
	const num lambda;

        enum Band {XBAND = 0, YBAND = 1};
	enum Spin {SPINUP = 0, SPINDOWN = 1};
	enum BC_Type { PBC, APBC_X, APBC_Y, APBC_XY };
	BC_Type bc;
        
	//hopping constants for XBAND and YBAND
	checkarray<num,2> hopHor;
	checkarray<num,2> hopVer;
	// sinh|cosh(dtau * hop..)
	checkarray<num,2> sinhHopHor;
	checkarray<num,2> sinhHopVer;
	checkarray<num,2> coshHopHor;
	checkarray<num,2> coshHopVer;
	PeriodicSquareLatticeNearestNeighbors spaceNeigh;
	PeriodicChainNearestNeighbors<1> timeNeigh;

	checkarray<MatNum, 2> propK;
	MatNum& propKx;
	MatNum& propKy;

	//for shifting green functions to obtain equivalency of symmetric Trotter decomposition
	//[checker board decomposition could be applied alternatively]
	std::array<MatNum, 2> propK_half;		//factor of -dtau/2 in exponential
	MatNum& propKx_half;
	MatNum& propKy_half;
	std::array<MatNum, 2> propK_half_inv;	//factor of +dtau/2 in exponential
	MatNum& propKx_half_inv;
	MatNum& propKy_half_inv;

	CubeCpx& g;
	CubeCpx& gFwd;
	CubeCpx& gBwd;

	//three component sdw-order parameter,
	//column indexes timeslice, row indexes site
	MatNum phi0;
	MatNum phi1;
	MatNum phi2;
	//evaluation of element-wise functions of phi:
	MatNum phiCosh;			// cosh(dtau * |phi|)
	MatNum phiSinh;			// sinh(dtau * |phi|) / |phi|


	num phiDelta;		//MC step size for field components
	//used to adjust phiDelta
	num targetAccRatio;
	num lastAccRatio;
	RunningAverage accRatioRA;

	//Observables:
	num normPhi;		//magnitude of averaged field
	num sdwSusc;		//spin-density-wave susceptibility

	checkarray<VecNum, 2> kOcc;		//Fermion occupation number in momentum space for x/y-band; site-index: k-vectors
	VecNum& kOccX;
	VecNum& kOccY;
	checkarray<VecNum, 2> kOccImag;
	VecNum& kOccXimag;
	VecNum& kOccYimag;

	checkarray<VecNum, 2> occ;     //Fermion occupation number in Real space for x/y-band; indexed by site
	VecNum& occX;
	VecNum& occY;
	checkarray<VecNum, 2> occImag;
	VecNum& occXimag;
	VecNum& occYimag;

	//pairing correlations near maximum range x_max = (L/2, L/2)
	num pairPlusMax;
	num pairMinusMax;
	num pairPlusMaximag;
	num pairMinusMaximag;
	//and between site 0 and any other site
	VecNum pairPlus;
	VecNum pairMinus;
	VecNum pairPlusimag;
	VecNum pairMinusimag;

	//fermion energy
	num fermionEkinetic;		//kinetic
	num fermionEkinetic_imag;
	num fermionEcouple;			//coupling
	num fermionEcouple_imag;

	//these for_each functions don't really work well with the class template
    template<typename Callable>
    void for_each_band(Callable func) {
    	func(XBAND);
    	func(YBAND);
    }
    template<typename Callable>
    void for_each_site(Callable func) {
    	for (uint32_t site = 0; site < N; ++site) {
    		func(site);
    	}
    }
    template<typename Callable>
    void for_each_timeslice(Callable func) {
    	for (uint32_t k = 1; k <= m; ++k) {
    		func(k);
    	}
    }
    template<typename CallableSiteTimeslice, typename V>
    V sumWholeSystem(CallableSiteTimeslice f, V init) {
#pragma omp parallel for reduction(+:init)
    	for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
    		for (uint32_t site = 0; site < N; ++site) {
				init += f(site, timeslice);
			}
    	}
    	return init;
    }
    template<typename CallableSiteTimeslice, typename V>
    V averageWholeSystem(CallableSiteTimeslice f, V init) {
    	V sum = sumWholeSystem(f, init);
    	return sum / num(m * N);
    }

    uint32_t coordsToSite(uint32_t x, uint32_t y) {
    	return y*L + x;
    }

    void setupRandomPhi();
    void setupPropK();			//compute e^(-dtau*K..) matrices by diagonalization
//	void setupPropK_direct();
//	void setupPropK_checkerboard();


    //the following are template function to allow applying them
    //to submatrices as well
    //"cb" = "checkerboard"
    //Reminder: These do not include chemical potential
    // with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to E^(sign * dtau * K_band) * A
    template <class Matrix>
    MatCpx cbLMultHoppingExp(const Matrix& A, Band band, int sign);
    // with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to A * E^(sign * dtau * K_band)
    template <class Matrix>
    MatCpx cbRMultHoppingExp(const Matrix& A, Band band, int sign);

    //the following take a 4Nx4N matrix A and effectively multiply B(k2,k1)
    //or its inverse to the left or right of it and return the result
    MatCpx checkerboardLeftMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardRightMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardLeftMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardRightMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1);

    //helpers for the checkerboardMultiplyFunctions
    MatCpx rightMultiplyBk(const MatCpx& orig, uint32_t k);		//multiply B(k,k-1) from right to orig, return result
    MatCpx rightMultiplyBkInv(const MatCpx& orig, uint32_t k);	//multiply B(k,k-1)^-1 from right to orig, return result
    MatCpx leftMultiplyBk(const MatCpx& orig, uint32_t k);		//multiply B(k,k-1) from left to orig, return result
    MatCpx leftMultiplyBkInv(const MatCpx& orig, uint32_t k);	//multiply B(k,k-1)^-1 from left to orig, return result


	MatCpx computeBmatSDW(uint32_t k2, uint32_t k1) const;			//compute B-matrix using dense matrix products

	virtual void updateInSlice(uint32_t timeslice);
	//this one does some adjusting of the box size from which new fields are chosen:
	virtual void updateInSliceThermalization(uint32_t timeslice);

	//functions used by updateInSlice:
	typedef VecNum::fixed<3> Phi;		//value of the three-component field at a single site and timeslice
	Phi proposeNewField(uint32_t site, uint32_t timeslice);
	num deltaSPhi(uint32_t site, uint32_t timeslice, Phi newphi);

	//compute the total value of the action associated with the field phi
	num phiAction();

	//wrappers to use to instantiate template functions of the base class
	struct sdwComputeBmat {
		DetSDW<TimeDisplaced,CheckerBoard>* parent;
		sdwComputeBmat(DetSDW<TimeDisplaced,CheckerBoard>* parent) :
			parent(parent)
		{ }
		MatCpx operator()(uint32_t gc, uint32_t k2, uint32_t k1) {
			(void)gc;
			assert(gc == 0);
			return parent->computeBmatSDW(k2, k1);
		}
	};

	struct sdwLeftMultiplyBmat {
		DetSDW<TimeDisplaced,CheckerBoard>* parent;
		sdwLeftMultiplyBmat(DetSDW<TimeDisplaced,CheckerBoard>* parent) :
			parent(parent)
		{ }
		MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
			(void)gc;
			assert(gc == 0);
			if (CheckerBoard) {
				return parent->checkerboardLeftMultiplyBmat(mat, k2, k1);
			} else {
				return parent->computeBmatSDW(k2, k1) * mat;
			}
		}
	};

	struct sdwRightMultiplyBmat {
		DetSDW<TimeDisplaced,CheckerBoard>* parent;
		sdwRightMultiplyBmat(DetSDW<TimeDisplaced,CheckerBoard>* parent) :
			parent(parent)
		{ }
		MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
			(void)gc;
			assert(gc == 0);
			if (CheckerBoard) {
				return parent->checkerboardRightMultiplyBmat(mat, k2, k1);
			} else {
				return mat * parent->computeBmatSDW(k2, k1);
			}
		}
	};

	struct sdwLeftMultiplyBmatInv {
		DetSDW<TimeDisplaced,CheckerBoard>* parent;
		sdwLeftMultiplyBmatInv(DetSDW<TimeDisplaced,CheckerBoard>* parent) :
			parent(parent)
		{ }
		MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
			(void)gc;
			assert(gc == 0);
			if (CheckerBoard) {
				return parent->checkerboardLeftMultiplyBmatInv(mat, k2, k1);
			} else {
				return arma::inv(parent->computeBmatSDW(k2, k1)) * mat;
			}
		}
	};

	struct sdwRightMultiplyBmatInv {
		DetSDW<TimeDisplaced,CheckerBoard>* parent;
		sdwRightMultiplyBmatInv(DetSDW<TimeDisplaced,CheckerBoard>* parent) :
			parent(parent)
		{ }
		MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
			(void)gc;
			assert(gc == 0);
			if (CheckerBoard) {
				return parent->checkerboardRightMultiplyBmatInv(mat, k2, k1);
			} else {
				return mat * arma::inv(parent->computeBmatSDW(k2, k1));
			}
		}
	};


public:
    // only functions that can pass the key to these functions have access
    // -- in this way access is granted only to select DetQMC methods
    template<class Archive>
    void saveContents(SerializeContentsKey const &sck, Archive &ar) {
    	DetModelGC<1,cpx>::saveContents(sck, ar);			//base class
    	serializeContentsCommon(sck, ar);
    }

    //after loadContents() a sweep must be performed before any measurements are taken:
    //else the green function would not be in a valid state
    template<class Archive>
    void loadContents(SerializeContentsKey const &sck, Archive &ar) {
    	DetModelGC<1,cpx>::loadContents(sck, ar);			//base class
    	serializeContentsCommon(sck, ar);
    	//the fields now have a valid state, update UdV-storage to start
    	//sweeping again
    	setupUdVStorage_skeleton(sdwComputeBmat(this));
    	//now: lastSweepDir == SweepDirection::Up --> the next sweep will be downwards
    }

    template<class Archive>
    void serializeContents(SerializeContentsKey const& sck, Archive& ar) {
    	Base::serializeContents(sck, ar);			//base class
		ar & phi0 & phi1 & phi2;
		ar & phiCosh & phiSinh;
		ar & phiDelta & targetAccRatio & lastAccRatio;
		ar & accRatioRA;
		ar & normPhi & sdwSusc;
		ar & kOcc & kOccImag;
    }
};

#endif /* DETSDW_H_ */
