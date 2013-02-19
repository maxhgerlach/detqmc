/*
 * neighbortable.h
 *
 *  Created on: Feb 19, 2013
 *      Author: gerlach
 */

#ifndef NEIGHBORTABLE_H_
#define NEIGHBORTABLE_H_

#include <cmath>
#include <cassert>
#include <vector>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include <armadillo>
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"

#include "tools.h"


//lattice directions are indexed as in +x,-x,+y,-y,+z,-z, ...
//enums provided for d <= 3, but code valid also for higher dimensions
enum class NeighDir : unsigned {
	XPLUS = 0, XMINUS = 1, YPLUS = 2, YMINUS = 3, ZPLUS = 4, ZMINUS = 5
};

class PeriodicCubicLatticeNearestNeighbors {
public:
	PeriodicCubicLatticeNearestNeighbors(unsigned d_, unsigned L_) :
		d(d_), L(L_), N(static_cast<unsigned>(uint_pow(L,d))), z(2*d),
		nearestNeigbors(z, N)
	{
		using std::floor;
		std::vector<unsigned> curCoords(d);     //holds the x, y, z coordinate components of the current site
		std::vector<unsigned> newCoords(d);     //newly calculated coords of the neighbor
		for (unsigned site = 0; site < N; ++site) {
			int reducedSite = site;
			for (int dim = d - 1; dim >= 0; --dim) {
				curCoords[dim] = unsigned(floor(reducedSite / uint_pow(L, dim)));
				reducedSite -= curCoords[dim] * uint_pow(L, dim);
			}
			assert(reducedSite == 0);
			for (unsigned dim = 0; dim < d; ++dim) {
				//neighbor in + direction, periodic
				newCoords = curCoords;
				newCoords[dim] = (newCoords[dim] + 1) % L;
				nearestNeigbors(dim * 2, site) = coordsToSite(newCoords);
				//neighbor in - direction, periodic
				newCoords = curCoords;
				newCoords[dim] = (newCoords[dim] - 1 + L) % L;
				nearestNeigbors(dim * 2 + 1, site) = coordsToSite(newCoords);
			}
		}
	}

	//get site index of nearest neighbor of site in NeighDir latticeDirection
	unsigned operator()(unsigned latticeDirection, unsigned site) const {
		assert(latticeDirection < z);
		assert(site < N);
		return nearestNeighbors(latticeDirection, site);
	}

	//iterators over the nearest neighbors of a site:
	auto beginNeighbors(unsigned site) -> decltype(const nearestNeighbors.begin_col(site)) const {
		return nearestNeighbors.begin_col(site);
	}

	auto endNeighbors(unsigned site) -> decltype(const nearestNeighbors.end_col(site)) const {
		return nearestNeighbors.end_col(site);
	}


	unsigned coordsToSite(const std::vector<unsigned>& coords) const {
		assert(coords.size() == d);
		unsigned site = 0;
		for (unsigned dim = 0; dim < d; ++dim) {
			site += coords[dim] * uint_pow(L, dim);
		}
		return site;
	}
private:
	unsigned d; 	//spatial dimension
	unsigned L;		//linear extent
	unsigned N;		//number of sites
	unsigned z;		//lattice coordination number
	//Neighbor table: columns index sites, rows index lattice directions
	typedef arma::Mat<unsigned> tableSites;
	tableSites nearestNeigbors;
};



#endif /* NEIGHBORTABLE_H_ */
