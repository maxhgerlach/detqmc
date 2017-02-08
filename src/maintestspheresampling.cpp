/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

// tests for sampling points in spherical coordinates; in the end we
// settled for picking a new point from a box centered around the
// previous point

#if defined (MAX_DEBUG) && ! defined(DUMA_NO_DUMA)
#include "dumapp.h"
#endif

/*
 * maintestspheresampling.cpp
 *
 *  Created on: Feb 14, 2014
 *      Author: gerlach
 */

#include <armadillo>
#include <cmath>
#include <array>
#include <tuple>
#include <cassert>
#include <string>
#include "rngwrapper.h"
#include "normaldistribution.h"

typedef double num;

RngWrapper rng;


// make a 2d histogram (pa, pb: series of orthogonal coordinate values), use nbins bins
// in each dimension, from -extend to +extend
inline arma::mat hist2d(const arma::vec& pa, const arma::vec& pb, uint32_t nbins, num extend) {
	assert(pa.n_elem == pb.n_elem);
	assert(extend > 0.);

	const num min = -extend;
	const num max = +extend;
	const num size = 2.*extend;
	auto getBin = [nbins,min,max,size](num val) -> uint32_t {
		return uint32_t(std::floor(((val - min) / size) * num(nbins)));
	};

	arma::mat out_hist = arma::zeros(nbins, nbins);
	for (uint32_t counter = 0; counter < pa.n_elem; ++counter) {
		num a = pa[counter];
		num b = pb[counter];
		uint32_t bina = getBin(a);
		uint32_t binb = getBin(b);
		out_hist(bina,binb) += 1.;
	}
	out_hist /= num(pa.n_elem);
	out_hist /= (size*size);

	return out_hist;
}

// TODO: Histogram, for instance only for z-values in the range corresponding to one bin -- this should appear uniform

inline arma::mat hist2d_onezbin(const arma::vec& px, const arma::vec& py, const arma::vec& pz, num target_zvalue, uint32_t nbins, num extend) {
	assert(px.n_elem == py.n_elem);
	assert(px.n_elem == pz.n_elem);
	assert(extend > 0.);

	const num min = -extend;
	const num max = +extend;
	const num size = 2.*extend;
	auto getBin = [nbins,min,max,size](num val) -> uint32_t {
		return uint32_t(std::floor(((val - min) / size) * num(nbins)));
	};

	const num target_bin = getBin(target_zvalue);

	arma::mat out_hist = arma::zeros(nbins, nbins);
	uint32_t total_count = 0;
	for (uint32_t counter = 0; counter < px.n_elem; ++counter) {
		uint32_t binx = getBin(px[counter]);
		uint32_t biny = getBin(py[counter]);
		uint32_t binz = getBin(pz[counter]);

		if (binz == target_bin) {
			out_hist(binx,biny) += 1.;
			total_count += 1;
		}
	}
	out_hist /= num(total_count);
	out_hist /= (size*size);

	return out_hist;
}


typedef std::array<num,3> Point;


// choose an independent random point from the 3d spherical volume of Rmax
inline Point drawRandomPointSphericalCoordinates(num Rmax) {
	using namespace std;
	num phi = rng.randRange(0., 2.*M_PI);
	num costheta = rng.randRange(-1., 1.);		// theta in (0, pi)
	num sintheta = sqrt(1. - pow(costheta, 2.));
	num rcubed = rng.randRange(0., pow(Rmax, 3.));
	num r = pow(rcubed, 1./3.);

	num x = r * cos(phi) * sintheta;
	num y = r * sin(phi) * sintheta;
	num z = r * costheta;

	return {{x, y, z}};
}

// choose an independent random point from the 3d spherical volume of Rmax,
// but choose the coordinate r from a uniform distribution, not r^3 -- should
// give botched results
inline Point drawRandomPointSphericalCoordinates_botched(num Rmax) {
	using namespace std;
	num phi = rng.randRange(0., 2.*M_PI);
	num costheta = rng.randRange(-1., 1.);		// theta in (0, pi)
	num sintheta = sqrt(1. - pow(costheta, 2.));
	num r = rng.randRange(0., Rmax);

	num x = r * cos(phi) * sintheta;
	num y = r * sin(phi) * sintheta;
	num z = r * costheta;

	return {{x, y, z}};
}


// choose new random point from a Markov chain -- importance sampling for the spherical
// volume, keeps a memory of the last point sampled
// 1) use the "box" modification technique"
class randomPointsFromMarkovChain_box {
	Point last_point;
	const num delta;
public:
	randomPointsFromMarkovChain_box()
		: last_point{{1.,1.,1.}}, delta(0.4)
	{ }

	Point operator()(num Rmax) {
		using namespace std;
		Point new_point;
		for (uint32_t i = 0; i < 3; ++i) {
			num d = rng.randRange(-delta, +delta);
			new_point[i] = last_point[i] + d;
		}
		num new_r = pow(pow(new_point[0], 2.) +
						pow(new_point[1], 2.) +
						pow(new_point[2], 2.),    1./2.);
		if (new_r > Rmax) {
			// do not change the current point, resample the last point
			return last_point;
		}
		else {
			// we just sample the uniform distribution 1.0/N,
			// Metropolis tells us to accept with prob = min(1., 1.) = 1
			last_point = new_point;
			return new_point;
		}
	}
};



// choose new random point from a Markov chain -- importance sampling for the spherical
// volume, keeps a memory of the last point sampled
// 1) rotate and then scale
class randomPointsFromMarkovChain_rotate_scale {
	NormalDistribution normal_distribution;
	Point last_point;
	const num scaleDelta;
	const num angleDelta;
public:
	randomPointsFromMarkovChain_rotate_scale()
		: normal_distribution(rng), last_point{{1.,1.,1.}},
		  scaleDelta(1.), angleDelta(0.8)
	{ }

	Point operator()(num Rmax) {
		using namespace std;

		//old orientation
		num x = last_point[0];
		num y = last_point[1];
		num z = last_point[2];
		//squares:
		num x2 = pow(x, 2);
		num y2 = pow(y, 2);
		num z2 = pow(z, 2);
		//squared length
		num r2 = x2 + y2 + z2;
		//length
		num r = sqrt(r2);

		//cubed length
		num r3 = pow(r, 3.0);
		num new_r3 = normal_distribution.get(scaleDelta, r3);
		if (new_r3 <= 0 or new_r3 > pow(Rmax, 3.0)) {
			// The gaussian-distributed new r^3 might be negative or zero or too big,
			// in that case the proposed new spin must
			// be rejected -- we sample r only from (0, Rmax).  In this case we just return the original spin again
			return last_point;
		} else {
			// otherwise we scale the spin appropriately and also change its orientation

			//new angular coordinates:
			num cosTheta = rng.rand01() * (1.0 - angleDelta) + angleDelta;     // \in [angleDelta, 1.0] since rand() \in [0, 1.0]
			num phi = rng.rand01() * 2.0 * M_PI;
			num sinTheta = sqrt(1.0 - pow(cosTheta, 2.0));
			num cosPhi = cos(phi);
			num sinPhi = sin(phi);

			// To find the new orientation, we first consider the normalized old spin
			num x2n = x2 / r2;
			num y2n = y2 / r2;
			num xn = x / r;
			num yn = y / r;
			num zn = z / r;
			// new spin (rotated so that cone from which the new spin is chosen has its center axis precisely aligned with the old spin);
			// this gives a normalized vector
			num newx = (sinTheta / (x2n+y2n)) * ((x2n*zn + y2n)*cosPhi + (zn-1)*xn*yn*sinPhi) + xn*cosTheta;
			num newy = (sinTheta / (x2n+y2n)) * ((zn-1)*xn*yn*cosPhi + (x2n + y2n*zn)*sinPhi) + yn*cosTheta;
			num newz = -sinTheta * (xn*cosPhi + yn*sinPhi) + zn*cosTheta;

			// Then we set the length of the new spin appropriately
			num new_r = pow(new_r3, (1.0 / 3.0));
			newx *= new_r;
			newy *= new_r;
			newz *= new_r;

			// Accept with prob = min(1., 1.)
			Point new_point {{newx, newy, newz}};
			last_point = new_point;
			return new_point;
		}
	}
};



template<class Function>
void createHistograms2d(Function getPoint,
		num Rmax, int samples, int nbins, const std::string& fileprefix) {
	arma::vec px(samples), py(samples), pz(samples);

	for (int i = 0; i < samples; ++i) {
		Point p = getPoint(Rmax);
		px[i] = p[0];
		py[i] = p[1];
		pz[i] = p[2];
	}

	arma::mat hist_xy = hist2d(px, py, nbins, Rmax);
	arma::mat hist_zx = hist2d(pz, px, nbins, Rmax);

	hist_xy.save(fileprefix + "_xy.csv", arma::csv_ascii);
	hist_zx.save(fileprefix + "_zx.csv", arma::csv_ascii);

	arma::mat hist_xy_z0 = hist2d_onezbin(px, py, pz, 0.0, nbins, Rmax);
	hist_xy_z0.save(fileprefix + "_z0_xy.csv", arma::csv_ascii);
}




int main() {
	const num Rmax = 4.0;
	const int samples = 1000*1000*10;
	const int nbins = 100;

	createHistograms2d(drawRandomPointSphericalCoordinates, Rmax, samples, nbins, "simple");
	createHistograms2d(drawRandomPointSphericalCoordinates_botched, Rmax, samples, nbins, "simple_botched");
	createHistograms2d(randomPointsFromMarkovChain_box(), Rmax, samples, nbins, "box");
	createHistograms2d(randomPointsFromMarkovChain_rotate_scale(), Rmax, samples, nbins, "rotate_scale");
}

