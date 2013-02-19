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

//base class for a model to be simulated by determinantal quantum Monte Carlo
//
//GreenComponents is the number of independent sectors of the Green's function,
//e.g. in the Hubbard model it is 2 for spin up and spin down


typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;
typedef arma::Mat<int> MatInt;
typedef arma::SpMat<num> SpMatNum;


template<unsigned GreenComponents>
class DetModel {
protected:
	DetModel(const ModelParams& pars);
public:
	virtual ~DetModel()
	{ }

	virtual unsigned getSystemN() const = 0;

	//Create a MetadataMap describing the parameters of the
	//simulated model
	virtual MetadataMap prepareModelMetadataMap() = 0;

	//perform measurements of all observables
    virtual void measure() = 0;

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
	virtual std::vector<ScalarObservable> getScalarObservables();
	virtual std::vector<VectorObservable> getVectorObservables();
	virtual std::vector<KeyValueObservable> getKeyValueObservables();

    //perform a sweep updating the auxiliary field with costly recomputations
    //of Green functions from scratch
    virtual void sweepSimple();

    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    virtual void sweep();

protected:
    typedef UdV<MatNum, VecNum> UdVnum;

    //update the auxiliary field and the green function in the single timeslice
	virtual void updateInSlice(unsigned timeslice) = 0;

    //functions to compute B-matrices for the different Green function sectors
    typedef std::function<MatNum(unsigned k2, unsigned k1)> FuncComputeBmat;
    std::array<FuncComputeBmat, GreenComponents> computeBmat;

    //functions to compute Green functions from scratch
    typedef std::function<MatNum(unsigned timeslice)> FuncComputeGreenFunction;
    std::array<FuncComputeGreenFunction, GreenComponents> computeGreenFunctionNaive;

    //functions that compute Green functions from UdV-decomposed matrices L/R
    //for a single timeslice
    //and update the member variables green -- and if desired -- greenFwd and greenBwd
    typedef std::function<void(unsigned targetSlice, const UdVnum& UdV_L, const UdVnum& UdV_R
    		                   )> FuncUpdateGreenFunctionUdV;
    std::array<FuncUpdateGreenFunctionUdV, GreenComponents> updateGreenFunctionUdV;
	//Given B(beta, tau) = V_l d_l U_l and B(tau, 0) = U_r d_r V_r
	//calculate a tuple of four NxN matrices (a,b,c,d) with
	// a = G(0), b = -(1-G(0))*B^(-1)(tau,0), c = B(tau,0)*G(0), d = G(tau)
	//b is the backward time-displaced Green function; c the forward time-
	//displaced Green function; d is the equal-time Green function
    //todo: get rid of this MatNum4 business
    typedef std::tuple<MatNum,MatNum,MatNum,MatNum> MatNum4;
    MatNum4 DetHubbard::greenFromUdV_timedisplaced(	const UdVnum& UdV_l, const UdVnum& UdV_r) const;
	//use a faster method that does not yield information about the time-displaced
	//Green functions
    MatNum greenFromUdV(const UdVnum& UdV_l, const UdVnum& UdV_r) const;

    //for each greenComponent call a function with the greenComponent as a parameter
    template<typename Callable>
    void for_each_gc(Callable func) {
    	for (unsigned gc = 0; gc < GreenComponents; ++gc) {
    		func(gc);
    	}
    }

    void setupUdVStorage();

    //helpers for sweep():
	void advanceDownGreen(unsigned l, unsigned greenComponent);
	void wrapDownGreen_timedisplaced(unsigned l, unsigned greenComponent);
	void wrapDownGreen(unsigned l, unsigned greenComponent);
	void advanceUpGreen(unsigned l, unsigned greenComponent);
	void advanceUpUpdateStorage(unsigned l, unsigned greenComponent);
	void wrapUpGreen_timedisplaced(unsigned l, unsigned greenComponent);
	void wrapUpGreen(unsigned l, unsigned greenComponent);
	void sweepUp();
	void sweepDown();
	typedef std::function<void(unsigned l, unsigned greenComponent)> FuncWrapGreen;
	FuncWrapGreen wrapDown;
	FuncWrapGreen wrapUp;


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
    std::array<CubeNum, GreenComponents> green;
    std::array<CubeNum, GreenComponents> greenFwd;
    std::array<CubeNum, GreenComponents> greenBwd;

    //The UdV-instances in UdVStorage will not move around much after setup, so storing
	//the (rather big) objects in the vector is fine
    virtual UdVnum give_eye_UdV() = 0;
    UdVnum eye_UdV;
	std::array<std::vector<UdVnum>, GreenComponents> UdVStorage;

    enum class SweepDirection: int {Up = 1, Down = -1};
	SweepDirection lastSweepDir;

	//observable handling -- these contain information about observables (such as their names)
	//as well as reference to their current value, which will be shared with simulation management
	//in a different class. The values reference there are to be updated here in the replica class.
	std::vector<ScalarObservable> obsScalar;
	std::vector<VectorObservable> obsVector;
	std::vector<KeyValueObservable> obsKeyValue;
};


template<unsigned GreenComponents>
DetModel::DetModel(const ModelParams& pars) :
	timedisplaced(pars.timedisplaced), beta(pars.beta), m(pars.m), s(pars.s), n(m / s), dtau(beta/m),
	obsScalar(), obsVector(), obsKeyValue()
{
	if (timedisplaced) {
		wrapUp = wrapUpGreen_timedisplaced;
		wrapDown = wrapDownGreen_timedisplaced;
	} else {
		wrapUp = wrapUpGreen;
		wrapDown = wrapDownGreen;
	}
}

template<unsigned GreenComponents>
void DetModel::setupUdVStorage() {
	eye_UdV = give_eye_UdV();

	auto setup = [this](unsigned gc) {
		std::vector<UdVnum>& storage = UdVStorage[gc];
		storage = std::vector<UdVnum>(n + 1);

		storage[0] = eye_UdV;
		storage[1] = udvNumDecompose(computeBmat[gc](s, 0));

		for (unsigned l = 1; l <= n - 1; ++l) {
			const MatNum& U_l = storage[l].U;
			const VecNum& d_l = storage[l].d;
			const MatNum& V_l = storage[l].V;
			MatNum B_lp1 = computeBmat[gc](s*(l + 1), s*l);
			UdVnum UdV_temp = udvNumDecompose((B_lp1 * U_l) * arma::diagmat(d_l));
			storage[l+1].U = UdV_temp.U;
			storage[l+1].d = UdV_temp.d;
			storage[l+1].V = UdV_temp.V * V_l;
		}
	};
	for_each_gc(setup);
	lastSweepDir = SweepDirection::Up;
}


template<unsigned GreenComponents>
void DetModel::sweepSimple() {
	for (unsigned timeslice = 1; timeslice <= m; ++timeslice) {
		for_each_gc( [this](unsigned gc) { computeGreenFunctionNaive[gc](timeslice); } );
		updateInSlice(timeslice);
	}
}


template<unsigned GreenComponents>
MatNum DetModel::greenFromUdV(const UdVnum& UdV_l, const UdVnum& UdV_r) const {
	//variable names changed according to labeling in notes
	const MatNum& V_l = UdV_l.U;   //!
	const VecNum& d_l = UdV_l.d;
	const MatNum& U_l = UdV_l.V;   //!
	const MatNum& U_r = UdV_r.U;
	const VecNum& d_r = UdV_r.d;
	const MatNum& V_r = UdV_r.V;

	using arma::inv; using arma::diagmat; using arma::eye;

	UdVnum UdV_temp = udvNumDecompose( inv(U_l * U_r) + diagmat(d_r) * (V_r * V_l) * diagmat(d_l) );

	MatNum green = inv(UdV_temp.V * U_l) * diagmat(1.0 / UdV_temp.d) * inv(U_r * UdV_temp.U);

	return green;
}


template<unsigned GreenComponents>
DetModel::MatNum4 DetHubbard::greenFromUdV_timedisplaced(
		const UdVnum& UdV_l, const UdVnum& UdV_r) const {
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
	UdVnum tempUdV = udvNumDecompose(temp);

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



//compute the green function in timeslice s*(l-1) from scratch with the help
//of the B-matrices computed before in the last up-sweep
template<unsigned GreenComponents>
void DetModel::advanceDownGreen(unsigned l, unsigned greenComponent) {
	std::vector<UdVnum>& storage = UdVStorage[greenComponent];

	MatNum B_l = computeBmat[greenComponent](s*l, s*(l - 1));

	//U_l, d_l, V_l correspond to B(beta,l*s*dtau) [set in the last step]
	const MatNum& U_l = storage[l].U;
	const VecNum& d_l = storage[l].d;
	const MatNum& V_l = storage[l].V;

	//UdV_L will correspond to B(beta,(l-1)*s*dtau)
	UdVnum UdV_L = udvNumDecompose(arma::diagmat(d_l) * (V_l * B_l));
	UdV_L.U = U_l * UdV_L.U;

	//UdV_R corresponds to B((l-1)*s*dtau,0) [set in last sweep]
	const UdVnum& UdV_R = storage[l - 1];
	unsigned next = s * (l - 1);
	updateGreenFunctionUdV[greenComponent](next, UdV_L, UdV_R);
	storage[l - 1] = UdV_L;
}


//compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
//also compute time-displaced Green functions
template<unsigned GreenComponents>
void DetModel::wrapDownGreen_timedisplaced(unsigned k, unsigned greenComponent) {
	MatNum B_k = computeBmat[greenComponent](k, k - 1);
	green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
	greenFwd[greenComponent].slice(k - 1) = arma::inv(B_k) * greenFwd[greenComponent].slice(k);
	greenBwd[greenComponent].slice(k - 1) = greenBwd[greenComponent](k) * B_k;
}

//compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
//only equal-time Green functions
template<unsigned GreenComponents>
void DetModel::wrapDownGreen(unsigned l, unsigned greenComponent) {
	MatNum B_k = computeBmat[greenComponent](k, k - 1);
	green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
}

//update the green function in timeslice s*(l+1) from scratch with the help
//of B-matrices computed before
template<unsigned GreenComponents>
void DetModel::advanceUpGreen(unsigned l, unsigned greenComponent) {
	std::vector<UdVnum>& storage = UdVStorage[greenComponent];

	MatNum B_lp1 = computeBmat[greenComponent](s*(l + 1), s*l);

	//The following is B(beta, (l+1)*s*dtau), valid from the last sweep
	const UdVnum& UdV_lp1 = storage[l + 1];

	//from the last step the following are B(l*s*dtau, 0):
	const MatNum& U_l = storage[l].U;
	const VecNum& d_l = storage[l].d;
	const MatNum& V_l = storage[l].V;

	//UdV_temp will be the new B((l+1)*s*dtau, 0):
	UdVnum UdV_temp = udvNumDecompose(((B_lp1 * U_l) * arma::diagmat(d_l)));
	UdV_temp.V *= V_l;

	unsigned next = s * (l + 1);
	updateGreenFunctionUdV[greenComponent](next, UdV_lp1, UdV_temp);

	//storage[l + 1] = UdV_temp;    //storage would be wrong after updateInSlice!
}

//Given B(l*s*dtau, 0) from the last step in the storage, compute
//B((l+1)*s*dtau, 0) and put it into storage
template<unsigned GreenComponents>
void DetModel::advanceUpUpdateStorage(unsigned l, unsigned greenComponent) {
	std::vector<UdVnum>& storage = UdVStorage[greenComponent];

	MatNum B_lp1 = computeBmat[greenComponent](s*(l + 1), s*l);
	//from the last step the following are B(l*s*dtau, 0):
	const MatNum& U_l = storage[l].U;
	const VecNum& d_l = storage[l].d;
	const MatNum& V_l = storage[l].V;
	//the new B((l+1)*s*dtau, 0):
	storage[l+1] = udvNumDecompose(((B_lp1 * U_l) * arma::diagmat(d_l)));
	storage[l+1].V *= V_l;
};

//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
//also handle the time-displaced Green functions
template<unsigned GreenComponents>
void DetModel::wrapUpGreen_timedisplaced(unsigned k, unsigned greenComponent) {
	MatNum B_kp1 = computeBmat[greenComponent](k + 1, k);
	green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
	greenFwd[greenComponent].slice(k + 1) = B_kp1 * greenFwd[greenComponent].slice(k);
	greenBwd[greenComponent].slice(k + 1) = greenBwd[greenComponent].slice(k) * arma::inv(B_kp1);
}

//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
//only compute equal-time Green functions
template<unsigned GreenComponents>
void DetModel::wrapUpGreen(unsigned k, unsigned greenComponent) {
	MatNum B_kp1 = computeBmat[greenComponent](k + 1, k);
	green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
}

template<unsigned GreenComponents>
void DetModel::sweepUp() {
	//We need to have computed the Green function for time slice k=0 so that the first
	//wrap-up step is correct.
	for (unsigned k = 1; k <= s-1; ++k) {
		for_each_gc( [this](unsigned gc) { wrapUp(k - 1, gc); } );
		updateInSlice(k);
	}
	//set storage at k=0 to unity for the upcoming sweep:
	for_each_gc( [this](unsigned gc) { UdVStorage[gc][0] = eye_UdV; } );
	for (unsigned l = 1; l < n; ++l) {
		for_each_gc( [this](unsigned gc) { advanceUpGreen(l-1, gc); });
		updateInSlice(l*s);
		for_each_gc( [this](unsigned gc) { advanceUpUpdateStorage(l-1, gc); });
		for (unsigned k = l*s + 1; k <= l*s + (s-1); ++k) {
			for_each_gc( [this](unsigned gc) { wrapUp(k - 1, gc); } );
			updateInSlice(k);
		}
	}
	updateInSlice(n*s);
	for_each_gc( [this](unsigned gc) { advanceUpUpdateStorage(n - 1, gc); } );
}

template<unsigned GreenComponents>
void DetModel::sweepDown() {
	//to compute green function for timeslice tau=beta:
	//we need VlDlUl = B(beta, beta) = I and UrDrVr = B(beta, 0).
	//The latter is given in storage slice m from the last sweep.
	for_each_gc( [this](unsigned gc) { updateGreenFunctionUdV[gc](eye_UdV, UdVStorage[gc][n]); } );
	for_each_gc( [this](unsigned gc) { UdVStorage[gc][n] = eye_UdV; } );
	for (unsigned l = n; l >= 1; --l) {
		updateInSlice(l*s);
		for (unsigned k = l*s - 1; k >= (l-1)*s + 1; --k) {
			for_each_gc( [this](unsigned gc) { wrapDown(k + 1, gc); } );
			updateInSlice(k);
		}
		//TODO: this will also compute the Green function at k=0, which technically is not necessary
		//but sensible for the following sweep up
		//TODO: alternatively just copy the k=m Green function to k=0  -- would that be up-to-date?
		for_each_gc( [this](unsigned gc) { advanceDownGreen(l, gc); } );
	}
}

template<unsigned GreenComponents>
void DetModel::sweep() {
	if (lastSweepDir == SweepDirection::Up) {
		sweepDown();
		lastSweepDir = SweepDirection::Down;
	} else if (lastSweepDir == SweepDirection::Down) {
		sweepUp();
		lastSweepDir = SweepDirection::Up;
	}
}





#endif /* DETMODEL_H_ */
