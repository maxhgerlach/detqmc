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

#include <vector>
#include <string>
#include <armadillo>

typedef double num;
typedef arma::Col<num> numvec;
typedef arma::Mat<num> nummat;
typedef arma::Cube<num> numcube;
typedef arma::Mat<int> intmat;
typedef arma::Mat<unsigned> tableSites;

class DetHubbard {
public:
	DetHubbard(num t, num U, num mu, unsigned L, unsigned d, num beta,
			unsigned m);
	virtual ~DetHubbard();

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


    void updateUnstabilized();

	enum class Spin: int {Up = +1, Down = -1};
protected:
	//parameters:
	num t;
	num U;
	num mu;
	unsigned L;
	unsigned d;
	unsigned latticeCoordination;      // 2 * d
	unsigned N;   // L ** d
	num beta;
	unsigned m;
	num dtau;     // beta / m
	num alpha;    // cosh(alpha) = exp(dtau U / 2)

	//Neighbor table: columns index sites, rows index lattice directions
	//as in +x,-x,+y,-y,+z,-z, ...
	tableSites nearestNeigbors;

	//Matrix representing the kinetic energy part of the hamiltonian: H_t
	//for spin up or spin down
	//TODO: no need to store tmat
	nummat tmat;
	// related propagator e ** (-dtau * tmat)
	nummat proptmat;


	//the following quantities vary during the course of the simulation


	//Auxiliary field represented by Ising spins, N entries of +/- 1.
	//There is one auxiliary field for each imaginary time slice. One time slice
	//corresponds to one column of the matrix. The time slices in auxfield are
	//indexed from 0 to m-1.
	intmat auxfield;

	//Equal imaginary time Green function
	//slices indexed n=0..m-1 correspond to time slices at dtau*n,
	//which are then indexed by sites in row and column.
	//One cube for each value of spinz.
	numcube gUp, gDn;


	//observables, values for the current auxiliary field; averaged over aux. field
	num occUp;          //occupation spin up
	num occDn;          //occupation spin down
	num occTotal;	    //total occupation
	num eKinetic;		//energy -t \sum_<i,j>,\sigma c^+_i,\sigma c_j,\sigma
	                    //       -mu \sum_j,\sigma n_j\sigma
	num ePotential;		//energy U \sum_i (n_i,up - 0.5) (n_i,down - 0.5)
	num eTotal;			//total energy


	//multiple observable handling
	std::vector<std::string> obsNames;
    std::vector<std::string> obsShorts;
    std::vector<num*> obsValPointers;			//pointers to variables updated in measure()
    unsigned obsCount;


	void createNeighborTable();
	unsigned coordsToSite(const std::vector<unsigned>& coords) const;


	void setupTmat();


	//compute e^{-scalar matrix}, matrix must be symmetric
	nummat computePropagator(num scalar, nummat matrix);

	//given the current auxiliary fields {s_n}, compute the matrix
	//  B_{s_n}(tau_2, tau_1) = \prod_{n = n2}^{n = n1 + 1} e^V(s_n) e^{-dtau T}
	// n2 = tau_2 / m > tau_1 / m = n1
	//So n1, n2 run from 0 to m.
	// e^{-dtau T} = proptmat
	//where the V(s_n) are computed for either the up or down Hubbard spins
	nummat computeBmat(unsigned n2, unsigned n1, Spin spinz);
	//calculate the B matrix for an arbitrary auxiliary field that need not match
	//the current one
	nummat computeBmat(unsigned n2, unsigned n1, Spin spinz,
			const nummat& arbitraryAuxfield);

	//Calculate (1 + B_s(tau, 0)*B_s(beta, tau))^(-1) from the given matrices
	//for the current aux field
	nummat computeGreenFunction(const nummat& bTau0, const nummat& bBetaTau);
	//Calculate the Green function from scratch for the given timeslice index
	//(tau = dtau * timeslice) for spin up or down
	nummat computeGreenFunction(unsigned timeslice, Spin spinz);

	//Calculate all entries of greenfctUp, greenfctDown from scratch for the
	//current aux field configuration.
	void computeAllGreenFunctions();

	//calculate det[1 + B_after(beta, 0)] / det[1 + B_before(beta,0)]
	//by brute force
	num weightRatioGeneric(const nummat& auxfieldBefore,
			const nummat& auxfieldAfter);

	//ratio of weighting determinants if a single auxiliary field
	//spin (at site in timeslice) is flipped from the current configuration.
	//use pre-stored Green functions.
	//formula independent of system size or number of time slices, in this form
	//specific to the Hubbard model
	num weightRatioSingleFlip(unsigned site, unsigned timeslice);


	//Update the stored Green function matrices to reflect the state after
	//the auxiliary field spin at site in timeslice has been flipped. This
	//function expects this->auxfield to be in the state before the flip.
	void updateGreenFunctionsAfterFlip(unsigned site, unsigned timeslice);
};

#endif /* DETHUBBARD_H_ */
