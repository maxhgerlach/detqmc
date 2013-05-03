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

class DetSDW;
std::unique_ptr<DetSDW> createDetSDW(RngWrapper& rng, ModelParams pars);


class DetSDW: public DetModelGC<1, cpx> {
	DetSDW(RngWrapper& rng, const ModelParams& pars	);
public:
	friend std::unique_ptr<DetSDW> createDetSDW(RngWrapper& rng, ModelParams pars);
	virtual ~DetSDW();
	virtual unsigned getSystemN() const;
	virtual MetadataMap prepareModelMetadataMap() const;
    virtual void measure();

    virtual void thermalizationOver();
protected:
	RngWrapper& rng;

	static const unsigned d = 2;
	static const unsigned z = 2*d;
	const unsigned L;
	const unsigned N;
	const num r;
	const num mu;
	const num c;
	const num u;
	const num lambda;

	PeriodicSquareLatticeNearestNeighbors spaceNeigh;
	PeriodicChainNearestNeighbors<1> timeNeigh;

	enum {XBAND = 0, YBAND = 1};
	std::array<MatNum, 2> propK;
	MatNum& propKx;
	MatNum& propKy;

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
	num phiSecond;		//second moment of averaged field
	num phiFourth;		//fourth moment of averaged field
	num binder;			//Binder cumulant: 1 - 3*phiFourth / (5 * phiSecond**2)
	num sdwSusc;		//spin-density-wave susceptibility
	std::array<VecNum, 2> kOcc;		//Fermion occupation number in momentum space for x/y-band; site-index: k-vectors
	VecNum& kOccX;
	VecNum& kOccY;
	std::array<VecNum, 2> kOccImag;
	VecNum& kOccXimag;
	VecNum& kOccYimag;

	std::array<VecNum, 2> occ;     //Fermion occupation number in Real space for x/y-band; indexed by site
	VecNum& occX;
	VecNum& occY;
	std::array<VecNum, 2> occImag;
	VecNum& occXimag;
	VecNum& occYimag;

    template<typename Callable>
    void for_each_band(Callable func) {
    	func(XBAND);
    	func(YBAND);
    }
    template<typename Callable>
    void for_each_site(Callable func) {
    	for (unsigned site = 0; site < N; ++site) {
    		func(site);
    	}
    }
    template<typename Callable>
    void for_each_timeslice(Callable func) {
    	for (unsigned k = 1; k <= m; ++k) {
    		func(k);
    	}
    }
    template<typename CallableSiteTimeslice, typename V>
    V sumWholeSystem(CallableSiteTimeslice f, V init) {
#pragma omp parallel for reduction(+:init)
    	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
    		for (unsigned site = 0; site < N; ++site) {
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

    void setupRandomPhi();
	void setupPropK();

	MatCpx computeBmatSDW(unsigned k2, unsigned k1) const;

	virtual void updateInSlice(unsigned timeslice);
	//this one does some adjusting of the box size from which new fields are chosen:
	virtual void updateInSliceThermalization(unsigned timeslice);

	//functions used by updateInSlice:
	typedef VecNum::fixed<3> Phi;		//value of the three-component field at a single site and timeslice
	Phi proposeNewField(unsigned site, unsigned timeslice);
	num deltaSPhi(unsigned site, unsigned timeslice, Phi newphi);

	//compute the total value of the action associated with the field phi
	num phiAction();

public:
    // only functions that can pass the key to this function have access
    // -- in this way access is granted only to DetQMC::serializeContents
    template<class Archive>
    void serializeContents(SerializeContentsKey const &sck, Archive &ar) {
    	DetModelGC<1,cpx>::serializeContents(sck, ar);			//base class
		ar & phi0 & phi1 & phi2;
		ar & phiCosh & phiSinh;
		ar & phiDelta & targetAccRatio & lastAccRatio;
		ar & accRatioRA;
		ar & normPhi & phiSecond & phiFourth & binder & sdwSusc;
		ar & kOcc & kOccImag;
    }
};

#endif /* DETSDW_H_ */
