/*
 * dethubbard.h
 *
 *  Created on: Dec 3, 2012
 *      Author: gerlach
 */

#ifndef DETHUBBARD_H_
#define DETHUBBARD_H_

/*
 * Determinantal Quantum Monte Carlo (FTQMC) for a standard single-band
 * Hubbard model on a periodic hypercube
 *
 * H_t = -t \sum_{<i,j>,sigma} c^+_{i,sigma} c^_{j,sigma} + h.c. - mu \sum_j (n_{j,up} + n_{j,down})
 * H_U = U \sum_j (n_{j,up} - 0.5) * (n_{j,down} - 0.5)
 *
 * Parameters
 *   t -- hopping energy
 *   U -- potential energy
 *   mu -- chemical potential
 *   L -- linear spatial extent
 *   d -- spatial dimension
 *
 *   beta -- inverse temperature (in units of 1/t, kB=1)
 *
 *   m -- number of imaginary time discretization levels (beta = m*dtau)
 *
 */

#include <functional>
#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <tuple>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include <armadillo>
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"
#include "rngwrapper.h"
#include "parameters.h"
#include "metadata.h"
#include "observable.h"

typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;

typedef arma::Mat<int> MatInt;

typedef arma::SpMat<num> SpMatNum;

typedef arma::Mat<unsigned> tableSites;

struct ModelParams;			//definition in parameters.h
class DetHubbard;		//defined below in this file

//factory function to init DetHubbard from parameter struct
//
//(in the future: possibly create such factories for different models)
//will do parameter checking etc
std::unique_ptr<DetHubbard> createDetHubbard(RngWrapper& rng, ModelParams pars);

class DetHubbard {
private:
	//only initialize with the "factory" function ::createDetHubbard() declared above.
	//Give a reference to the RNG instance to be used
	DetHubbard(RngWrapper& rng, const ModelParams& pars);
public:
	friend std::unique_ptr<DetHubbard> createDetHubbard(RngWrapper& rng, ModelParams pars);
	virtual ~DetHubbard();

	unsigned getSystemN() const;

	//Create a MetadataMap describing the parameters of the
	//simulated model
	MetadataMap prepareModelMetadataMap();

	//perform measurements of all observables
    void measure();

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
	std::vector<ScalarObservable> getScalarObservables();
	std::vector<VectorObservable> getVectorObservables();
	std::vector<KeyValueObservable> getKeyValueObservables();

    //perform a sweep updating the auxiliary field with costly recomputations
    //of Green functions from scratch
    void sweepSimple();

    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    void sweep();

	enum class Spin: int {Up = +1, Down = -1};
protected:
	RngWrapper& rng;
	//parameters:
	const bool checkerboard;
	const num t;			//hopping energy scale
	const num U;			//interaction energy scale
	const num mu;			//chemical potential
	const unsigned L;	//linear lattice size
	const unsigned d;	//spatial dimension of lattice
	const unsigned z;   // lattice coordination number, 2 * d
	const unsigned N;   // L ** d
	const num beta;		//inverse temperature
	const unsigned m;	//number of imaginary time discretization steps (time slices) beta*m=dtau
	const unsigned s;	//interval between time slices where the Green-function is calculated from scratch
	const unsigned n;	//number of time slices where the Green-function is calculated from scratch n*s*dtau=beta
	const num dtau;     // beta / m
	const num alpha;    // cosh(alpha) = exp(dtau U / 2)

	//Neighbor table: columns index sites, rows index lattice directions
	//as in +x,-x,+y,-y,+z,-z, ...
	tableSites nearestNeigbors;
	enum class NeighDir : unsigned {
		XPLUS = 0, XMINUS = 1, YPLUS = 2, YMINUS = 3
	};

	std::function<MatNum(unsigned k2, unsigned k1, Spin spinz)> computeBmatFunc;


	//Matrix representing the kinetic energy part of the hamiltonian: H_t
	//for spin up or spin down -- tmat
	// related propagator e ** (-dtau * tmat):
	MatNum proptmat;

	//checker board decomposition matrices (2d square lattice only ATM)
	//used in computation of propagator
//	SpMatNum checkerX_a, checkerX_b;
//	SpMatNum checkerY_a, checkerY_b;

	//the following quantities vary during the course of the simulation


	//Auxiliary field represented by Ising spins, N entries of +/- 1.
	//There is one auxiliary field for each imaginary time slice. One time slice
	//corresponds to one column of the matrix. The time slices in auxfield are
	//indexed from 0 to m. So auxfield.col(n) refers to the timeslice dtau*n.
	//Most code, however, only uses timeslices >= 1 ! Don't rely on auxfield.col(0)
	//being valid.
	MatInt auxfield;

	//Equal imaginary time Green function
	//slices indexed k=0..m correspond to time slices at dtau*k,
	//which are then indexed by sites in row and column.
	//One cube for each value of spinz.
	//The Green functions for k=0 are conceptually equal to those for k=m.
	//Most code, however, only uses timeslices k >= 1 ! Don't rely on g*.slice(0)
	//being valid.
	CubeNum gUp, gDn;
	
	//Imaginary time displaced Green function
	// "forward" corresponds to G(tau, 0)
	// "backward" corresponds to G(0, tau)
	//the indexing works the same way as for the equal time case
	CubeNum gFwdUp, gFwdDn;
	CubeNum gBwdUp, gBwdDn;

	//matrices used in the computation of B-matrices decomposed into
	//(U,d,V) = (orthogonal matrix, diagonal matrix elements, row-normalized triangular matrix.
	//Initialize at beginning of simulation by member function setupUdVStorage()
	struct UdV {
		MatNum U;
		VecNum d;
		MatNum V;
		//default constructor: leaves everything empty
		UdV() : U(), d(), V() {}
		//specify matrix size: initialize to identity
		UdV(unsigned size) :
			U(arma::eye(size,size)), d(arma::ones(size)), V(arma::eye(size,size))
		{ }
	};
	UdV eye_UdV;	// U = d = V = 1
	static UdV udvDecompose(const MatNum& mat);				//wraps Armadillo functions
	//The UdV-instances in UdVStorage will not move around much after setup, so storing
	//the (rather big) objects in the vector is fine
	std::vector<UdV> UdVStorageUp;
	std::vector<UdV> UdVStorageDn;

	enum class SweepDirection: int {Up = 1, Down = -1};
	SweepDirection lastSweepDir;


	//observables, values for the current auxiliary field; averaged over aux. field
	num occUp;          //occupation spin up
	num occDn;          //occupation spin down
	num occTotal;	    //total occupation
	num eKinetic;		//energy -t \sum_<i,j>,\sigma c^+_i,\sigma c_j,\sigma
	                    //       -mu \sum_j,\sigma n_j\sigma
	num ePotential;		//energy U \sum_i (n_i,up - 0.5) (n_i,down - 0.5)
	num eTotal;			//total energy
	num occDouble;		//double occupation < nUp * nDown >
	num localMoment;	//local Moment: <m^2> = <(nUp - nDown)^2>
	num suscq0;			//q=0 susceptibility of z component of magnetization (nUp - nDown)

	//the same for vector observables, averaged over timeslices with the current auxiliary field
	VecNum zcorr;		//correlation function of magnetization density at site 0 with all other sites, z component

	//for d=2: the Fourier transform of the timedisplaced Green function for k=(pi/2, 2*pi/3)
	VecNum gf;			//values for different time-displacements
	VecNum gf_dt;		//time-displacements, where the values are evaluated



	//observable handling -- these contain information about observables (such as their names)
	//as well as reference to their current value, which will be shared with simulation management
	//in a different class. The values reference there are to be updated here in the replica class.
	std::vector<ScalarObservable> obsScalar;
	std::vector<VectorObservable> obsVector;
	std::vector<KeyValueObservable> obsKeyValue;

	void createNeighborTable();
	unsigned coordsToSite(const std::vector<unsigned>& coords) const;


	void setupRandomAuxfield();
	void setupPropTmat_direct();
	void setupPropTmat_checkerboard();
	void setupUdVStorage();


	//compute e^{-scalar matrix}, matrix must be symmetric
	MatNum computePropagator_direct(num scalar, const MatNum& matrix) const;

	//given the current auxiliary fields {s_n}, compute the matrix
	// B_{s_n}(tau_2, tau_1) = \prod_{n = n2}^{n = n1 + 1} e^V(s_n) e^{-dtau T}
	// n2 = tau_2 / m > tau_1 / m = n1
	//So n1, n2 run from 0 to m.
	// e^{-dtau T} = proptmat
	//Here the V(s_n) are computed for either the up or down Hubbard spins.
	//These functions naively multiply the matrices, which can be unstable.
	MatNum computeBmat_direct(unsigned k2, unsigned k1, Spin spinz) const;

	//compute the latter using a checker board decomposition with systematic
	//error of O[dtau^2]
	MatNum computeBmat_checkerBoard(unsigned k2, unsigned k1, Spin spinz) const;

	//Calculate (1 + B_s(tau, 0)*B_s(beta, tau))^(-1) from the given matrices
	//for the current aux field.
	//These functions perform a naive matrix product and inversion.
	MatNum computeGreenFunctionNaive(const MatNum& bTau0, const MatNum& bBetaTau) const;
	//Calculate the Green function from scratch for the given timeslice index
	//(tau = dtau * timeslice) for spin up or down
	MatNum computeGreenFunctionNaive(unsigned timeslice, Spin spinz) const;


	//calculate det[1 + B_after(beta, 0)] / det[1 + B_before(beta,0)]
	//by brute force and naively computed B-matrices
//	num weightRatioGenericNaive(const MatInt& auxfieldBefore,
//			const MatInt& auxfieldAfter) const;

	//ratio of weighting determinants if a single auxiliary field
	//spin (at site in timeslice) is flipped from the current configuration.
	//use pre-stored Green functions.
	//formula independent of system size or number of time slices, in this form
	//specific to the Hubbard model
	num weightRatioSingleFlip(unsigned site, unsigned timeslice) const;


	//Update the stored Green function matrices to reflect the state after
	//the auxiliary field spin at site in timeslice has been flipped. This
	//function expects this->auxfield to be in the state before the flip.
	void updateGreenFunctionWithFlip(unsigned site, unsigned timeslice);

	//update the HS auxiliary field and the green function in the single timeslice
	void updateInSlice(unsigned timeslice);

	//Given B(beta, tau) = V_l d_l U_l and B(tau, 0) = U_r d_r V_r
	//calculate a tuple of four NxN matrices (a,b,c,d) with
	// a = G(0), b = -(1-G(0))*B^(-1)(tau,0), c = B(tau,0)*G(0), d = G(tau)
	//b is the backward time-displaced Green function; c the forward time-
	//displaced Green function; d is the equal-time Green function
	typedef std::tuple<MatNum,MatNum,MatNum,MatNum> MatNum4;
	MatNum4 greenFromUdV_timedisplaced(const UdV& UdV_l, const UdV& UdV_r) const;

	//use a faster method that does not yield information about the time-displaced
	//Green functions
	MatNum greenFromUdV(const UdV& UdV_l, const UdV& UdV_r) const;

	void debugCheckBeforeSweepDown();
	void debugCheckBeforeSweepUp();
	void debugCheckGreenFunctions();
};

#endif /* DETHUBBARD_H_ */
