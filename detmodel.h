/*
 * detmodel.h
 *
 *  Created on: Feb 18, 2013
 *      Author: gerlach
 */

#ifndef DETMODEL_H_
#define DETMODEL_H_


#include <array>
#include <functional>
#include <utility>
#include <memory>
#include <vector>
#include <tuple>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include <armadillo>
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"

#include "parameters.h"
#include "observable.h"
#include "udv.h"
#include "metadata.h"
#include "timing.h"

typedef std::complex<double> cpx;

typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;
typedef arma::Mat<int> MatInt;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<num> SpMatNum;
typedef std::tuple<MatNum,MatNum,MatNum,MatNum> MatNum4;
typedef UdV<num> UdVnum;


//base class for a model to be simulated by determinantal quantum Monte Carlo
//

//purely abstract base class
class DetModel {
public:
	virtual ~DetModel() { };
	virtual unsigned getSystemN() const = 0;

	//Create a MetadataMap describing the parameters of the
	//simulated model
	virtual MetadataMap prepareModelMetadataMap() const = 0;

	//perform measurements of all observables
    virtual void measure() = 0;

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
	virtual std::vector<ScalarObservable> getScalarObservables() = 0;
	virtual std::vector<VectorObservable> getVectorObservables() = 0;
	virtual std::vector<KeyValueObservable> getKeyValueObservables() = 0;

    //perform a sweep updating the auxiliary field with costly re-computations
    //of Green functions from scratch
    virtual void sweepSimple() = 0;
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepSimpleThermalization() = 0;


    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    virtual void sweep() = 0;
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepThermalization() = 0;

    //notify DetModel that thermalization has finished
    //do nothing by default
    virtual void thermalizationOver() {
    }
};


//GreenComponents is the number of independent sectors of the Green's function,
//e.g. in the Hubbard model it is 2 for spin up and spin down

template<unsigned GreenComponents, typename ValueType = num>
class DetModelGC : public DetModel {
protected:
	DetModelGC(const ModelParams& pars, unsigned greenComponentSize);
public:
	virtual ~DetModelGC()
	{ }

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
	virtual std::vector<ScalarObservable> getScalarObservables();
	virtual std::vector<VectorObservable> getVectorObservables();
	virtual std::vector<KeyValueObservable> getKeyValueObservables();

    //perform a sweep updating the auxiliary field with costly re-computations
    //of Green functions from scratch
    virtual void sweepSimple();
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepSimpleThermalization();

    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    virtual void sweep();
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepThermalization();
protected:
    typedef ValueType V;
    typedef arma::Mat<ValueType> MatV;
    typedef arma::Col<ValueType> VecV;
    typedef arma::Cube<ValueType> CubeV;
    typedef UdV<ValueType> UdVV;
    typedef std::tuple<MatV,MatV,MatV,MatV> MatV4;

    //update the auxiliary field and the green function in the single timeslice
	virtual void updateInSlice(unsigned timeslice) = 0;
	//separate function to be called during thermalization, by default just do the
	//same; a derived class may override this to introduce an adaptive behavior
	virtual void updateInSliceThermalization(unsigned timeslice) {
		updateInSlice(timeslice);
	}

    //functions to compute B-matrices for the different Green function sectors
    typedef std::function<MatV(unsigned k2, unsigned k1)> FuncComputeBmat;
    std::array<FuncComputeBmat, GreenComponents> computeBmat;

    //functions that compute Green functions from UdV-decomposed matrices L/R
    //for a single timeslice
    //and update the member variables green -- and if desired -- greenFwd and greenBwd
    typedef std::function<void(unsigned targetSlice, const UdVV& UdV_L, const UdVV& UdV_R
    		                   )> FuncUpdateGreenFunctionUdV;
    std::array<FuncUpdateGreenFunctionUdV, GreenComponents> updateGreenFunctionUdV;
	//Given B(beta, tau) = V_l d_l U_l and B(tau, 0) = U_r d_r V_r
	//calculate a tuple of four NxN matrices (a,b,c,d) with
	// a = G(0), b = -(1-G(0))*B^(-1)(tau,0), c = B(tau,0)*G(0), d = G(tau)
	//b is the backward time-displaced Green function; c the forward time-
	//displaced Green function; d is the equal-time Green function
    //todo: get rid of this MatNum4 business
    MatV4 greenFromUdV_timedisplaced(const UdVV& UdV_l, const UdVV& UdV_r) const;
	//use a faster method that does not yield information about the time-displaced
	//Green functions
    MatV greenFromUdV(const UdVV& UdV_l, const UdVV& UdV_r) const;

    //for each greenComponent call a function with the greenComponent as a parameter
    template<typename Callable>
    void for_each_gc(Callable func) {
    	for (unsigned gc = 0; gc < GreenComponents; ++gc) {
    		func(gc);
    	}
    }

    //only call *after* computeBmat[] is valid, i.e. in a derived class:
    void setupUdVStorage();

    //helpers for sweep(), sweepThermalization():
	void advanceDownGreen(unsigned l, unsigned greenComponent);
	void wrapDownGreen_timedisplaced(unsigned k, unsigned greenComponent);
	void wrapDownGreen(unsigned k, unsigned greenComponent);
	void advanceUpGreen(unsigned l, unsigned greenComponent);
	void advanceUpUpdateStorage(unsigned l, unsigned greenComponent);
	void wrapUpGreen_timedisplaced(unsigned k, unsigned greenComponent);
	void wrapUpGreen(unsigned k, unsigned greenComponent);
	//functions to do the wrapping are set at runtime (depending on whether timedisplaced routines are used or not)
	typedef std::function<void(unsigned k, unsigned greenComponent)> FuncWrapGreen;
	FuncWrapGreen wrapUp;
	FuncWrapGreen wrapDown;
	//these receive as a template parameter the function to call for updates in a slice
	template <class CallableUpdateInSlice> void sweepUp(CallableUpdateInSlice funcUpdateInSlice);
	template <class CallableUpdateInSlice> void sweepDown(CallableUpdateInSlice funcUpdateInSlice);

	//Green component size, e.g. sz == N for the Hubbard model
	unsigned sz;

    //some simulation parameters are already relevant for member functions implemented
    //in this base class, the rest will only be used in derived classes
	bool timedisplaced;
	const num beta;		//inverse temperature
	const unsigned m;	//number of imaginary time discretization steps (time slices) beta*m=dtau
	const unsigned s;	//interval between time slices where the Green-function is calculated from scratch
	const unsigned n;	//number of time slices where the Green-function is calculated from scratch n*s*dtau=beta
	const num dtau;     // beta / m


    //equal-imaginary-time and time-displaced Green's functions
	//slices indexed k=0..m correspond to time slices at dtau*k,
	//which are then indexed by sites in row and column.
	//Most code, however, only uses timeslices k >= 1 ! Don't rely on g*.slice(0)
	//being valid.
	//The Green functions for k=0 are conceptually equal to those for k=m.
    std::array<CubeV, GreenComponents> green;
    std::array<CubeV, GreenComponents> greenFwd;
    std::array<CubeV, GreenComponents> greenBwd;

    //The UdV-instances in UdVStorage will not move around much after setup, so storing
	//the (rather big) objects in the vector is fine
    UdVV eye_UdV;
	std::array<std::vector<UdVV>, GreenComponents> UdVStorage;

    enum class SweepDirection: int {Up = 1, Down = -1};
	SweepDirection lastSweepDir;

	//observable handling -- these contain information about observables (such as their names)
	//as well as reference to their current value, which will be shared with simulation management
	//in a different class. The values reference there are to be updated here in the replica class.
	std::vector<ScalarObservable> obsScalar;
	std::vector<VectorObservable> obsVector;
	std::vector<KeyValueObservable> obsKeyValue;
};


template<unsigned GC, typename V>
DetModelGC<GC,V>::DetModelGC(const ModelParams& pars, unsigned greenComponentSize) :
	computeBmat{}, updateGreenFunctionUdV{}, wrapUp(), wrapDown(),
	sz(greenComponentSize),
	timedisplaced(pars.timedisplaced), beta(pars.beta), m(pars.m), s(pars.s), n(m / s), dtau(beta/m),
	green{}, greenFwd{}, greenBwd{}, eye_UdV(sz), UdVStorage{},
	lastSweepDir(SweepDirection::Up),
	obsScalar(), obsVector(), obsKeyValue()
{
	if (timedisplaced) {
		wrapUp = [this](unsigned k, unsigned gc) { this->wrapUpGreen_timedisplaced(k, gc); };
		wrapDown = [this](unsigned k, unsigned gc) { this->wrapDownGreen_timedisplaced(k, gc); };
		for_each_gc( [this](unsigned gc) {
			updateGreenFunctionUdV[gc] = [this, gc](unsigned targetSlice,
					                                const UdVV& UdV_L, const UdVV& UdV_R) {
				std::tie(std::ignore, greenBwd[gc].slice(targetSlice),
						 greenFwd[gc].slice(targetSlice), green[gc].slice(targetSlice)) =
					this->greenFromUdV_timedisplaced(UdV_L, UdV_R);
			};
		} );
	} else {
		wrapUp = [this](unsigned k, unsigned gc) { this->wrapUpGreen(k, gc); };
		wrapDown = [this](unsigned k, unsigned gc) { this->wrapDownGreen(k, gc); };
		for_each_gc( [this](unsigned gc) {
			updateGreenFunctionUdV[gc] = [this, gc](unsigned targetSlice,
					                                const UdVV& UdV_L, const UdVV& UdV_R) {
				green[gc].slice(targetSlice) = this->greenFromUdV(UdV_L, UdV_R);
			};
		} );
	}
}

template<unsigned GC, typename V>
std::vector<ScalarObservable> DetModelGC<GC,V>::getScalarObservables() {
	return obsScalar;
}

template<unsigned GC, typename V>
std::vector<VectorObservable> DetModelGC<GC,V>::getVectorObservables() {
	return obsVector;
}

template<unsigned GC, typename V>
std::vector<KeyValueObservable> DetModelGC<GC,V>::getKeyValueObservables() {
	return obsKeyValue;
}





template<unsigned GC, typename V>
void DetModelGC<GC,V>::setupUdVStorage() {
	timing.start("setupUdVStorage");
	auto setup = [this](unsigned gc) {
		std::vector<UdVV>& storage = UdVStorage[gc];
		storage = std::vector<UdVV>(n + 1);

		storage[0] = eye_UdV;
		storage[1] = udvDecompose(computeBmat[gc](s, 0));

		for (unsigned l = 1; l <= n - 1; ++l) {
			const MatV& U_l = storage[l].U;
			const VecV& d_l = storage[l].d;
			const MatV& V_l = storage[l].V;
			MatV B_lp1 = computeBmat[gc](s*(l + 1), s*l);
			UdVV UdV_temp = udvDecompose<V>((B_lp1 * U_l) * arma::diagmat(d_l));
			storage[l+1].U = UdV_temp.U;
			storage[l+1].d = UdV_temp.d;
			storage[l+1].V = UdV_temp.V * V_l;
		}
	};
	for_each_gc(setup);
	lastSweepDir = SweepDirection::Up;
	timing.stop("setupUdVStorage");
}


//warning: the thermalization version below is almost a copy of this
template<unsigned GC, typename V>
void DetModelGC<GC,V>::sweepSimple() {
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for_each_gc( [this, timeslice](unsigned gc) {
			green[gc].slice(timeslice) =
					arma::inv(arma::eye(sz,sz) + computeBmat[gc](timeslice, 0) *
					                              computeBmat[gc](m, timeslice));
		});
		updateInSlice(timeslice);
	}
}

//warning: this is almost a copy of sweepSimple() defined above
template<unsigned GC, typename V>
void DetModelGC<GC,V>::sweepSimpleThermalization() {
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for_each_gc( [this, timeslice](unsigned gc) {
			green[gc].slice(timeslice) =
					arma::inv(arma::eye(sz,sz) + computeBmat[gc](timeslice, 0) *
					                              computeBmat[gc](m, timeslice));
		});
		updateInSliceThermalization(timeslice);
	}
}



template<unsigned GC, typename V>
typename DetModelGC<GC,V>::MatV DetModelGC<GC,V>::greenFromUdV(const UdVV& UdV_l, const UdVV& UdV_r) const {
	timing.start("greenFromUdV");
	//variable names changed according to labeling in notes
	const MatV& V_l = UdV_l.U;   //!
	const VecV& d_l = UdV_l.d;
	const MatV& U_l = UdV_l.V;   //!
	const MatV& U_r = UdV_r.U;
	const VecV& d_r = UdV_r.d;
	const MatV& V_r = UdV_r.V;

	using arma::inv; using arma::diagmat; using arma::eye;

	UdVV UdV_temp = udvDecompose<V>( inv(U_l * U_r) + diagmat(d_r) * (V_r * V_l) * diagmat(d_l) );

	MatV green = inv(UdV_temp.V * U_l) * diagmat(1.0 / UdV_temp.d) * inv(U_r * UdV_temp.U);

	timing.stop("greenFromUdV");

	return green;
}


template<unsigned GC, typename V>
typename DetModelGC<GC,V>::MatV4 DetModelGC<GC,V>::greenFromUdV_timedisplaced(
		const UdVV& UdV_l, const UdVV& UdV_r) const {
	timing.start("greenFromUdV_timedisplaced");

	//Ul vs Vl to be compatible with labeling in the notes
	const MatV& Ul = UdV_l.V;   //!
	const VecV& dl = UdV_l.d;
	const MatV& Vl = UdV_l.U;   //!
	const MatV& Ur = UdV_r.U;
	const VecV& dr = UdV_r.d;
	const MatV& Vr = UdV_r.V;

	unsigned sz = Ul.n_rows;

	//submatrix view helpers for 2*N x 2*N matrices
	auto upleft = [sz](MatV& m) {
		return m.submat(0,0, sz-1,sz-1);
	};
	auto upright = [sz](MatV& m) {
		return m.submat(0,sz, sz-1,2*sz-1);
	};
	auto downleft = [sz](MatV& m) {
		return m.submat(sz,0, 2*sz-1,sz-1);
	};
	auto downright = [sz](MatV& m) {
		return m.submat(sz,sz, 2*sz-1,2*sz-1);
	};

	MatV temp(2*sz,2*sz);
	upleft(temp)    = arma::inv(Vr * Vl);
	upright(temp)   = arma::diagmat(dl);
	downleft(temp)  = arma::diagmat(-dr);
	downright(temp) = arma::inv(Ul * Ur);
	UdVV tempUdV = udvDecompose<V>(temp);

	MatV left(2*sz,2*sz);
	upleft(left) = arma::inv(Vr);
	upright(left).zeros();
	downleft(left).zeros();
	downright(left) = arma::inv(Ul);

	MatV right(2*sz,2*sz);
	upleft(right) = arma::inv(Vl);
	upright(right).zeros();
	downleft(right).zeros();
	downright(right) = arma::inv(Ur);

	MatV result = (left * arma::inv(tempUdV.V)) * arma::diagmat(1.0 / tempUdV.d)
					* (arma::inv(tempUdV.U) * right);

	timing.stop("greenFromUdV_timedisplaced");

	return MatV4(upleft(result), upright(result),
				   downleft(result), downright(result));
}



//compute the green function in timeslice s*(l-1) from scratch with the help
//of the B-matrices computed before in the last up-sweep
template<unsigned GC, typename V>
void DetModelGC<GC,V>::advanceDownGreen(unsigned l, unsigned greenComponent) {
	timing.start("advanceDownGreen");

	std::vector<UdVV>& storage = UdVStorage[greenComponent];

	MatV B_l = computeBmat[greenComponent](s*l, s*(l - 1));

	//U_l, d_l, V_l correspond to B(beta,l*s*dtau) [set in the last step]
	const MatV& U_l = storage[l].U;
	const VecV& d_l = storage[l].d;
	const MatV& V_l = storage[l].V;

	//UdV_L will correspond to B(beta,(l-1)*s*dtau)
	UdVV UdV_L = udvDecompose<V>(arma::diagmat(d_l) * (V_l * B_l));
	UdV_L.U = U_l * UdV_L.U;

	//UdV_R corresponds to B((l-1)*s*dtau,0) [set in last sweep]
	const UdVV& UdV_R = storage[l - 1];
	unsigned next = s * (l - 1);
	updateGreenFunctionUdV[greenComponent](next, UdV_L, UdV_R);
	storage[l - 1] = UdV_L;

	timing.stop("advanceDownGreen");
}

//compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
//also compute time-displaced Green functions
template<unsigned GC, typename V>
void DetModelGC<GC,V>::wrapDownGreen_timedisplaced(unsigned k, unsigned greenComponent) {
	timing.start("wrapDownGreen_timedisplaced");
	MatV B_k = computeBmat[greenComponent](k, k - 1);
	green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
	greenFwd[greenComponent].slice(k - 1) = arma::inv(B_k) * greenFwd[greenComponent].slice(k);
	greenBwd[greenComponent].slice(k - 1) = greenBwd[greenComponent](k) * B_k;
	timing.stop("wrapDownGreen_timedisplaced");
}

//compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
//only equal-time Green functions
template<unsigned GC, typename V>
void DetModelGC<GC,V>::wrapDownGreen(unsigned k, unsigned greenComponent) {
	timing.start("wrapDownGreen");
	MatV B_k = computeBmat[greenComponent](k, k - 1);
	green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
	timing.stop("wrapDownGreen");
}

//update the green function in timeslice s*(l+1) from scratch with the help
//of B-matrices computed before
template<unsigned GC, typename V>
void DetModelGC<GC,V>::advanceUpGreen(unsigned l, unsigned greenComponent) {
	timing.start("advanceUpGreen");

	std::vector<UdVV>& storage = UdVStorage[greenComponent];

	MatV B_lp1 = computeBmat[greenComponent](s*(l + 1), s*l);

	//The following is B(beta, (l+1)*s*dtau), valid from the last sweep
	const UdVV& UdV_lp1 = storage[l + 1];

	//from the last step the following are B(l*s*dtau, 0):
	const MatV& U_l = storage[l].U;
	const VecV& d_l = storage[l].d;
	const MatV& V_l = storage[l].V;

	//UdV_temp will be the new B((l+1)*s*dtau, 0):
	UdVV UdV_temp = udvDecompose<V>(((B_lp1 * U_l) * arma::diagmat(d_l)));
	UdV_temp.V *= V_l;

	unsigned next = s * (l + 1);
	updateGreenFunctionUdV[greenComponent](next, UdV_lp1, UdV_temp);

	//storage[l + 1] = UdV_temp;    //storage would be wrong after updateInSlice!

	timing.stop("advanceUpGreen");
}

//Given B(l*s*dtau, 0) from the last step in the storage, compute
//B((l+1)*s*dtau, 0) and put it into storage
template<unsigned GC, typename V>
void DetModelGC<GC,V>::advanceUpUpdateStorage(unsigned l, unsigned greenComponent) {
	timing.start("advanceUpUpdateStorage");

	std::vector<UdVV>& storage = UdVStorage[greenComponent];

	MatV B_lp1 = computeBmat[greenComponent](s*(l + 1), s*l);
	//from the last step the following are B(l*s*dtau, 0):
	const MatV& U_l = storage[l].U;
	const VecV& d_l = storage[l].d;
	const MatV& V_l = storage[l].V;
	//the new B((l+1)*s*dtau, 0):
	storage[l+1] = udvDecompose<V>(((B_lp1 * U_l) * arma::diagmat(d_l)));
	storage[l+1].V *= V_l;

	timing.stop("advanceUpUpdateStorage");
};

//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
//also handle the time-displaced Green functions
template<unsigned GC, typename V>
void DetModelGC<GC,V>::wrapUpGreen_timedisplaced(unsigned k, unsigned greenComponent) {
	timing.start("wrapUpGreen_timedisplaced");
	MatV B_kp1 = computeBmat[greenComponent](k + 1, k);
	green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
	greenFwd[greenComponent].slice(k + 1) = B_kp1 * greenFwd[greenComponent].slice(k);
	greenBwd[greenComponent].slice(k + 1) = greenBwd[greenComponent].slice(k) * arma::inv(B_kp1);
	timing.stop("wrapUpGreen_timedisplaced");
}

//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
//only compute equal-time Green functions
template<unsigned GC, typename V>
void DetModelGC<GC,V>::wrapUpGreen(unsigned k, unsigned greenComponent) {
	timing.start("wrapUpGreen");
	MatV B_kp1 = computeBmat[greenComponent](k + 1, k);
	green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
	timing.stop("wrapUpGreen");
}

template<unsigned GC, typename V>
template<class CallableUpdateInSlice>
void DetModelGC<GC,V>::sweepUp(CallableUpdateInSlice funcUpdateInSlice) {
	//We need to have computed the Green function for time slice k=0 so that the first
	//wrap-up step is correct.
	for (unsigned k = 1; k <= s-1; ++k) {
		for_each_gc( [this, k](unsigned gc) { wrapUp(k - 1, gc); } );
		funcUpdateInSlice(k);
	}
	//set storage at k=0 to unity for the upcoming sweep:
	for_each_gc( [this](unsigned gc) { UdVStorage[gc][0] = eye_UdV; } );
	for (unsigned l = 1; l < n; ++l) {
		for_each_gc( [this, l](unsigned gc) { this->advanceUpGreen(l-1, gc); });
		funcUpdateInSlice(l*s);
		for_each_gc( [this, l](unsigned gc) { this->advanceUpUpdateStorage(l-1, gc); });
		for (unsigned k = l*s + 1; k <= l*s + (s-1); ++k) {
			for_each_gc( [this, k](unsigned gc) { wrapUp(k - 1, gc); } );
			funcUpdateInSlice(k);
		}
	}
	funcUpdateInSlice(n*s);
	for_each_gc( [this](unsigned gc) { this->advanceUpUpdateStorage(n - 1, gc); } );
}

template<unsigned GC, typename V>
template <class CallableUpdateInSlice>
void DetModelGC<GC,V>::sweepDown(CallableUpdateInSlice funcUpdateInSlice) {
	//to compute green function for timeslice tau=beta:
	//we need VlDlUl = B(beta, beta) = I and UrDrVr = B(beta, 0).
	//The latter is given in storage slice m from the last sweep.
	for_each_gc( [this](unsigned gc) { updateGreenFunctionUdV[gc](m, eye_UdV, UdVStorage[gc][n]); } );
	for_each_gc( [this](unsigned gc) { UdVStorage[gc][n] = eye_UdV; } );
	for (unsigned l = n; l >= 1; --l) {
		funcUpdateInSlice(l*s);
		for (unsigned k = l*s - 1; k >= (l-1)*s + 1; --k) {
			for_each_gc( [this, k](unsigned gc) { wrapDown(k + 1, gc); } );
			funcUpdateInSlice(k);
		}
		//TODO: this will also compute the Green function at k=0, which technically is not necessary
		//but sensible for the following sweep up
		//TODO: alternatively just copy the k=m Green function to k=0  -- would that be up-to-date?
		for_each_gc( [this, l](unsigned gc) { this->advanceDownGreen(l, gc); } );
	}
}

template<unsigned GC, typename V>
void DetModelGC<GC,V>::sweep() {
	timing.start("sweep");

	if (lastSweepDir == SweepDirection::Up) {
		sweepDown([this](unsigned k){ this->updateInSlice(k); });
		lastSweepDir = SweepDirection::Down;
	} else if (lastSweepDir == SweepDirection::Down) {
		sweepUp([this](unsigned k){ this->updateInSlice(k); });
		lastSweepDir = SweepDirection::Up;
	}

	timing.stop("sweep");
}

template<unsigned GC, typename V>
void DetModelGC<GC,V>::sweepThermalization() {
	timing.start("sweep");

	if (lastSweepDir == SweepDirection::Up) {
		sweepDown([this](unsigned k){ this->updateInSliceThermalization(k); });
		lastSweepDir = SweepDirection::Down;
	} else if (lastSweepDir == SweepDirection::Down) {
		sweepUp([this](unsigned k){ this->updateInSliceThermalization(k); });
		lastSweepDir = SweepDirection::Up;
	}

	timing.stop("sweep");
}


//Special handling to allow passing either 'm' or 'dtau', but not both.
//Also check that 's' is set correctly.
ModelParams updateTemperatureParameters(ModelParams pars);


//compute e^{-scalar matrix}, matrix must be symmetric
MatNum computePropagator(num scalar, const MatNum& matrix);



template <typename Matrix> inline
void debugSaveMatrix(const Matrix& matrix, const std::string& basename) {
	matrix.save(basename + ".csv", arma::csv_ascii);
}

#endif /* DETMODEL_H_ */
