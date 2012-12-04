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
 *   U -- interaction energy
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
#include <armadillo>

typedef double num;
typedef arma::Col<num> numvec;
typedef arma::Mat<num> nummat;
typedef arma::Mat<unsigned> tableSites;

class DetHubbard {
public:
	DetHubbard(num t, num U, num mu, unsigned L, unsigned d, num beta,
			unsigned m);
	virtual ~DetHubbard();

	enum class Spin: int {Up = +1, Down = -1};
protected:
	//parameters:
	num t;
	num U;
	num mu;
	unsigned L;
	unsigned d;
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
	nummat tmat;
	// related propagator e ** (-dtau * tmat)
	nummat proptmat;


	//the following varies during the course of the simulation

	//Auxiliary field represented by Ising spins, N entries of +/- 1.
	//There is one auxiliary field for each imaginary time slice. One time slice
	//corresponds to one column of the matrix. The time slices in auxfield are
	//indexed from 0 to m-1.
	nummat auxfield;


	void createNeighborTable();
	unsigned coordsToSite(const std::vector<unsigned>& coords);


	void setupTmat();


	//compute e^{-scalar matrix}, matrix must be symmetric
	nummat computePropagator(num scalar, nummat matrix);

	//given the current auxiliary fields {s_n}, compute the matrix
	//  B_{s_n}(tau_2, tau_1) = \prod_{n = n2}^{n = n1 + 1} e^V(s_n) e^{-dtau T}
	// n2 = tau_2 / m > tau_1 / m = n1
	//So n1, n2 run from 0 to m.
	// e^{-dtau T} = proptmat
	//where the V(s_n) are computed for either the up or down Hubbard spins
	nummat computeBmat(unsigned n2, unsigned n1, Spin spinComponent);
};

#endif /* DETHUBBARD_H_ */
