/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

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
#include <utility>
#include <armadillo>
#include <string>

#include "tools.h"

#include "boost_serialize_armadillo.h"


//lattice directions are indexed as in +x,-x,+y,-y,+z,-z, ...
//enums provided for d <= 3, but code valid also for higher dimensions
enum NeighDir {
    XPLUS = 0, XMINUS = 1, YPLUS = 2, YMINUS = 3, ZPLUS = 4, ZMINUS = 5
};

inline std::string neighToString(NeighDir neigh) {
	switch (neigh) {
	case XPLUS: return "XPLUS";
	case XMINUS: return "XMINUS";
	case YPLUS: return "YPLUS";
	case YMINUS: return "YMINUS";
	case ZPLUS: return "ZPLUS";
	case ZMINUS: return "ZMINUS";
	default: return "#NODIR";
	}
}

typedef arma::Mat<uint32_t> tableSites;


class PeriodicCubicLatticeNearestNeighbors {
public:
    PeriodicCubicLatticeNearestNeighbors(uint32_t d_, uint32_t L_) :
        d(d_), L(L_), N(static_cast<uint32_t>(uint_pow(L,d))), z(2*d),
        nearestNeighbors(z, N)
    {
        using std::floor;
        std::vector<uint32_t> curCoords(d);     //holds the x, y, z coordinate components of the current site
        std::vector<uint32_t> newCoords(d);     //newly calculated coords of the neighbor
        for (uint32_t site = 0; site < N; ++site) {
            uint32_t reducedSite = site;
            for (int idim = int(d) - 1; idim >= 0; --idim) {
                uint32_t dim = uint32_t(idim);
                curCoords[dim] = uint32_t(floor(reducedSite / uint_pow(L, dim)));
                reducedSite -= curCoords[dim] * uint_pow(L, dim);
            }
            assert(reducedSite == 0);
            for (uint32_t dim = 0; dim < d; ++dim) {
                //neighbor in + direction, periodic
                newCoords = curCoords;
                newCoords[dim] = (newCoords[dim] + 1) % L;
                nearestNeighbors(dim * 2, site) = coordsToSite(newCoords);
                //neighbor in - direction, periodic
                newCoords = curCoords;
                newCoords[dim] = (newCoords[dim] - 1 + L) % L;
                nearestNeighbors(dim * 2 + 1, site) = coordsToSite(newCoords);
            }
        }
    }
    virtual ~PeriodicCubicLatticeNearestNeighbors() {
    }

    //get site index of nearest neighbor of site in NeighDir latticeDirection
    uint32_t operator()(uint32_t latticeDirection, uint32_t site) const {
        assert(latticeDirection < z);
        assert(site < N);
        return nearestNeighbors(latticeDirection, site);
    }
    uint32_t operator()(NeighDir latticeDirection, uint32_t site) const {
        return operator()((uint32_t) latticeDirection, site);
    }

    //iterators over the nearest neighbors of a site:
    auto beginNeighbors(uint32_t site) const -> tableSites::const_col_iterator {
        return nearestNeighbors.begin_col(site);
    }

    auto endNeighbors(uint32_t site) const -> tableSites::const_col_iterator {
        return nearestNeighbors.end_col(site);
    }


    uint32_t coordsToSite(const std::vector<uint32_t>& coords) const {
        assert(coords.size() == d);
        uint32_t site = 0;
        for (uint32_t dim = 0; dim < d; ++dim) {
            site += coords[dim] * uint_pow(L, dim);
        }
        return site;
    }
protected:
    uint32_t d;     //spatial dimension
    uint32_t L;     //linear extent
    uint32_t N;     //number of sites
    uint32_t z;     //lattice coordination number
    //Neighbor table: columns index sites, rows index lattice directions
    tableSites nearestNeighbors;

protected:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar & d & L & N & z;
        ar & nearestNeighbors;
    }
};


class PeriodicSquareLatticeNearestNeighbors : public PeriodicCubicLatticeNearestNeighbors {
public:
    PeriodicSquareLatticeNearestNeighbors(uint32_t L_)
        : PeriodicCubicLatticeNearestNeighbors(2, L_)
    {
    }
    void save(const std::string& file = "neighbors.csv") {
    	nearestNeighbors.save(file, arma::csv_ascii);
    }
};


enum class ChainDir : uint32_t {
    PLUS = 0, MINUS = 1
};


template <uint32_t startWith=0>     //start the indexing with 0 or maybe 1 or something else...
class PeriodicChainNearestNeighbors : public PeriodicCubicLatticeNearestNeighbors {
public:
    PeriodicChainNearestNeighbors(uint32_t L_)
        : PeriodicCubicLatticeNearestNeighbors(1, L_)
    {
        //nearestNeighbors.save("timeneighbors.csv", arma::csv_ascii);
    }

    uint32_t operator()(uint32_t latticeDirection, uint32_t site) const {
        return PeriodicCubicLatticeNearestNeighbors::operator()(latticeDirection, site - startWith)
               + startWith;
    }
    uint32_t operator()(NeighDir latticeDirection, uint32_t site) const {
        return PeriodicChainNearestNeighbors::operator()((uint32_t) latticeDirection, site);
    }
    uint32_t operator()(ChainDir latticeDirection, uint32_t site) const {
        return PeriodicChainNearestNeighbors::operator()((uint32_t) latticeDirection, site);
    }
private:
    //hide these as they are not easy to make work correctly and are not important

    typedef PeriodicCubicLatticeNearestNeighbors Base;
    //iterators over the nearest neighbors of a site:
    auto beginNeighbors(uint32_t site) const -> tableSites::const_col_iterator {
    	return Base::beginNeighbors(site);
    }
    auto endNeighbors(uint32_t site) const -> tableSites::const_col_iterator {
    	return Base::endNeighbors(site);
    }
};



#endif /* NEIGHBORTABLE_H_ */
