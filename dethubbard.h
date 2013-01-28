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

#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <tuple>
#include <armadillo>
#include "rngwrapper.h"
#include "parameters.h"
#include "metadata.h"

typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;

typedef arma::Mat<int> MatInt;

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

    //get values of observables normalized by system size:
    //obs corresponds to an observable additional to the energy
    //if values of obsIndex > 0 are supported, we measure multiple, different
    //observables (up to obsIndex == getNumberOfObservables() - 1 )
    unsigned getNumberOfObservables() const;
    num obsNormalized(unsigned obsIndex = 0) const;
    virtual std::string getObservableName(unsigned obsIndex = 0) const;
    virtual std::string getObservableShort(unsigned obsIndex = 0) const;

    //the same for vector observables
    unsigned getNumberOfVectorObservables() const;
    VecNum vecObsNormalized(unsigned obsIndex = 0) const;
    virtual std::string getVectorObservableName(unsigned obsIndex = 0) const;
    virtual std::string getVectorObservableShort(unsigned obsIndex = 0) const;


    //perform a sweep updating the auxiliary field with costly recomputations
    //of Green functions from scratch
    void sweepSimple();

    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time
    void sweep();

	enum class Spin: int {Up = +1, Down = -1};
protected:
	RngWrapper& rng;
	//parameters:
	num t;			//hopping energy scale
	num U;			//interaction energy scale
	num mu;			//chemical potential
	unsigned L;	//linear lattice size
	unsigned d;	//spatial dimension of lattice
	unsigned z;   // lattice coordination number, 2 * d
	unsigned N;   // L ** d
	num beta;		//inverse temperature
	unsigned m;	//number of imaginary time discretization steps (time slices) beta*m=dtau
	unsigned s;	//interval between time slices where the Green-function is calculated from scratch
	unsigned n;	//number of time slices where the Green-function is calculated from scratch n*s*dtau=beta
	num dtau;     // beta / m
	num alpha;    // cosh(alpha) = exp(dtau U / 2)

	//Neighbor table: columns index sites, rows index lattice directions
	//as in +x,-x,+y,-y,+z,-z, ...
	tableSites nearestNeigbors;

	//Matrix representing the kinetic energy part of the hamiltonian: H_t
	//for spin up or spin down -- tmat
	// related propagator e ** (-dtau * tmat):
	MatNum proptmat;


	//the following quantities vary during the course of the simulation


	//Auxiliary field represented by Ising spins, N entries of +/- 1.
	//There is one auxiliary field for each imaginary time slice. One time slice
	//corresponds to one column of the matrix. The time slices in auxfield are
	//indexed from 0 to m-1. So auxfield.col(n) refers to the timeslice dtau*(n+1).
	MatInt auxfield;

	//Equal imaginary time Green function
	//slices indexed n=0..m-1 correspond to time slices at dtau*(n+1),
	//which are then indexed by sites in row and column.
	//One cube for each value of spinz.
	//The Green functions for n=0 are equal to those for n=m.
	CubeNum gUp, gDn;
	
	//Imaginary time displaced Green function, tau > 0
	// "forward" corresponds to G(tau, 0)
	// "backward" corresponds to G(0, tau)
	//the indexing works the same way as for the equal time case
	CubeNum gFwdUp, gFwdDn;
	CubeNum gBwdUp, gBwdDn;

	//matrices used in the computation of B-matrices with singular value decomposition.
	//(U,d,V) = (orthogonal matrix, diagonal matrix elements, orthogonal matrix)
	//Initialize at beginning of simulation by member function setupUdVStorage()
	struct UdV {
		MatNum U;
		VecNum d;
		MatNum V;
	};
	UdV eye_UdV;	// U = d = V = 1
	static UdV svd(const MatNum& mat);				//wraps arma::svd()
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


	//multiple observable handling
	std::vector<std::string> obsNames;
    std::vector<std::string> obsShorts;
    //references to variables updated in measure():
    std::vector<std::reference_wrapper<const num>> obsValRefs;
    unsigned obsCount;

    std::vector<std::string> vecObsNames;
    std::vector<std::string> vecObsShorts;
    //references to variables updated in measure():
    std::vector<std::reference_wrapper<const VecNum>> vecObsValRefs;
    unsigned vecObsCount;

	void createNeighborTable();
	unsigned coordsToSite(const std::vector<unsigned>& coords) const;


	void setupRandomAuxfield();
	void setupPropTmat();
	void setupUdVStorage();


	//compute e^{-scalar matrix}, matrix must be symmetric
	MatNum computePropagator(num scalar, const MatNum& matrix) const;

	//given the current auxiliary fields {s_n}, compute the matrix
	// B_{s_n}(tau_2, tau_1) = \prod_{n = n2}^{n = n1 + 1} e^V(s_n) e^{-dtau T}
	// n2 = tau_2 / m > tau_1 / m = n1
	//So n1, n2 run from 0 to m.
	// e^{-dtau T} = proptmat
	//Here the V(s_n) are computed for either the up or down Hubbard spins.
	//These functions naively multiply the matrices, which can be unstable.
	MatNum computeBmatNaive(unsigned bn2, unsigned n1, Spin spinz) const;
	//calculate the B matrix for an arbitrary auxiliary field that need not match
	//the current one
	MatNum computeBmatNaive(unsigned n2, unsigned n1, Spin spinz,
			const MatInt& arbitraryAuxfield) const;

	//Calculate (1 + B_s(tau, 0)*B_s(beta, tau))^(-1) from the given matrices
	//for the current aux field.
	//These functions perform a naive matrix product and inversion.
	MatNum computeGreenFunctionNaive(const MatNum& bTau0, const MatNum& bBetaTau) const;
	//Calculate the Green function from scratch for the given timeslice index
	//(tau = dtau * timeslice) for spin up or down
	MatNum computeGreenFunctionNaive(unsigned timeslice, Spin spinz) const;


	//calculate det[1 + B_after(beta, 0)] / det[1 + B_before(beta,0)]
	//by brute force and naively computed B-matrices
	num weightRatioGenericNaive(const MatInt& auxfieldBefore,
			const MatInt& auxfieldAfter) const;

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
};

#endif /* DETHUBBARD_H_ */
