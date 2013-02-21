/*
 * detsdw.h
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#ifndef DETSDW_H_
#define DETSDW_H_

#include <array>
#include <vector>
#include <complex>
#include "rngwrapper.h"
#include "detmodel.h"
#include "neighbortable.h"

typedef std::complex<num> cpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::Col<cpx> VecCpx;
typedef arma::Cube<cpx> CubeCpx;

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
	MatNum phi1;
	MatNum phi2;
	MatNum phi3;
	//evaluation of element-wise functions of phi:
	MatNum phiCosh;			// cosh(|phi|)
	MatNum phiSinh;			// sinh(|phi|) / |phi|


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

    void setupRandomPhi();
	void setupPropK();

	MatCpx computeBmatSDW(unsigned k2, unsigned k1) const;

	virtual void updateInSlice(unsigned timeslice);
};

#endif /* DETSDW_H_ */
