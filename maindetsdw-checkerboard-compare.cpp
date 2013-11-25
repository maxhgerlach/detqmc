/*
 * maindetsdw-checkerboard-compare.cpp
 *
 *  Created on: Nov 20, 2013
 *      Author: gerlach
 */


#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <vector>
#include "git-revision.h"
#include "metadata.h"
#include "detsdw.h"
#include "parameters.h"
#include "rngwrapper.h"


// class derived from DetSDW that only adds a read-only interface to the Green's
// function
// // this makes thinks horribly hard
//template <bool CheckerBoard>
//class DetSDW_getGreen : public DetSDW<false, CheckerBoard> {
//public:
//	typedef DetSDW<false, CheckerBoard> Base;
//	DetSDW_getGreen(const Base& base) : Base(base) {
//	}
//	const CubeCpx& get_g() {
//		return g;
//	}
//private:
//	using Base::g;
//};



ModelParams setupParameters(bool checkerboard_, uint32_t m_, num dtau_, uint32_t s_, std::string bc_) {
	ModelParams params;

#define SET(par, val) { \
		params.par = val; \
		params.specified.insert(#par); \
    }

	SET(model, "sdw");
	SET(timedisplaced, false);
	SET(checkerboard, checkerboard_);
	SET(r, 1.0);
	SET(txhor, -1.0);
	SET(txver, -0.5);
	SET(tyhor, 0.5);
	SET(tyver, 1.0);
	SET(mu, 0.5);
	SET(L, 4);
	SET(d, 2);
	SET(m, m_);
	SET(dtau, dtau_);
	SET(s, s_);
	SET(accRatio, 0.5);
	SET(bc, bc_);
	SET(rescale, false);
#undef SET

	return params;
}


int main() {
	std::string bc = "pbc";
	std::vector<num> dtau_values = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5};
	num beta = 10;
	uint32_t s = 10;
	std::vector<uint32_t> m_values;
	m_values.reserve(dtau_values.size());
	for (num dtau : dtau_values) {
		m_values.push_back(std::ceil(beta / dtau));
	}

	for (uint32_t i = 0; i < dtau_values.size(); ++i) {
		num dtau = dtau_values[i];
		uint32_t m = m_values[i];

		// create one instance with checkerboard decomposition and one without
		// same parameters and same rng seed

		ModelParams pars_cb = setupParameters(true, m, dtau, s, bc);
		RngWrapper rng_cb(1111);
		typedef DetSDW<false, true> Cb;
		std::unique_ptr<DetModel> cb_model_ptr = createDetSDW(rng_cb, pars_cb);
		std::unique_ptr<Cb> sdw_cb = std::unique_ptr<Cb>(dynamic_cast<Cb*>(cb_model_ptr.release()));

		ModelParams pars_reg = setupParameters(false, m, dtau, s, bc);
		RngWrapper rng_reg(1111);
		typedef DetSDW<false, false> Reg;
		std::unique_ptr<DetModel> reg_model_ptr = createDetSDW(rng_reg, pars_reg);
		std::unique_ptr<Reg> sdw_reg = std::unique_ptr<Reg>(dynamic_cast<Reg*>(reg_model_ptr.release()));

		// do two sweeps (one up, one down) for each instance

		sdw_cb->sweepThermalization();
		sdw_cb->sweepThermalization();
		sdw_reg->sweepThermalization();
		sdw_reg->sweepThermalization();

		// compute differences of the Green's functions

		CubeNum g_abs_diff = arma::abs(sdw_reg->get_green() - sdw_cb->get_green());
		//CubeNum g_abs_diff = arma::abs(sdw_reg->get_green());
		//CubeNum g_abs_diff = arma::real(sdw_reg->get_green());
		num max_abs_diff = g_abs_diff.max();
		num mean_abs_diff = 0.0;
		for (uint32_t k = 1; k <= m; ++k) {
			num mean_k = 0.0;
			uint32_t elems = g_abs_diff.slice(k).size();
			for (uint32_t i = 0; i < elems; ++i) {
				mean_k += g_abs_diff.slice(k)[i];
			}
			mean_k /= num(elems);
			mean_abs_diff += mean_k;
		}
		mean_abs_diff /= num(m);

		std::cout << dtau << '\t' << max_abs_diff << '\t' << mean_abs_diff << std::endl;
	}

	return 0;
}


