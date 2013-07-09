/*
 * dethubbard.cpp
 *
 *  Created on: Dec 3, 2012
 *      Author: gerlach
 */

#include <functional>
#include <vector>
#include <cmath>
#include <complex>
#include <cassert>
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#include "tools.h"
#include "exceptions.h"
#include "rngwrapper.h"
#include "dethubbard.h"
#include "parameters.h"

//using std::acosh;    //Intel compiler chokes with std::acosh



std::unique_ptr<DetModel> createDetHubbard(RngWrapper& rng, ModelParams pars) {
    pars = updateTemperatureParameters(pars);

    //check parameters: passed all that are necessary
    using namespace boost::assign;
    std::vector<std::string> neededModelPars;
    neededModelPars += "t", "U", "mu", "L", "d", "checkerboard";
    for (auto p = neededModelPars.cbegin(); p != neededModelPars.cend(); ++p) {
        if (pars.specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

    if (pars.checkerboard and pars.L % 2 != 0) {
        throw ParameterWrong("Checker board decomposition only supported for even linear lattice sizes");
    }
    if (pars.checkerboard and pars.d != 2) {
        throw ParameterWrong("Checker board decomposition only supported for 2d lattices");
    }
    if (pars.bc != "pbc") {
        throw ParameterWrong("Boundary conditions " + pars.bc + " not supported for Hubbard model (only pbc)");
    }

    //check that only positive values are passed for certain parameters
#define IF_NOT_POSITIVE(x) if (pars.specified.count(#x) > 0 and pars.x <= 0)
#define CHECK_POSITIVE(x)   {                                           \
                                IF_NOT_POSITIVE(x) {                    \
                                    throw ParameterWrong(#x, pars.x);   \
                                }                                       \
                            }
    CHECK_POSITIVE(L);
    CHECK_POSITIVE(d);
#undef CHECK_POSITIVE
#undef IF_NOT_POSITIVE

    //since pars is not a constant expression, we need this stupid if:
    if (pars.timedisplaced == true and pars.checkerboard == true) {
        return std::unique_ptr<DetModel>(new DetHubbard<true,true>(rng, pars));
    } else
    if (pars.timedisplaced == true and pars.checkerboard == false) {
        return std::unique_ptr<DetModel>(new DetHubbard<true,false>(rng, pars));
    } else
    if (pars.timedisplaced == false and pars.checkerboard == true) {
        return std::unique_ptr<DetModel>(new DetHubbard<false,true>(rng, pars));
    } else
    if (pars.timedisplaced == false and pars.checkerboard == false) {
        return std::unique_ptr<DetModel>(new DetHubbard<false,false>(rng, pars));
    } else {
        //this can't be reached
        return 0;
    }
}


template <bool TD, bool CB>
DetHubbard<TD,CB>::DetHubbard(RngWrapper& rng_, const ModelParams& pars) :
        DetModelGC<2,num,TD>(pars, static_cast<uint32_t>(uint_pow(pars.L,pars.d))),
        rng(rng_),
        checkerboard(CB),
        timedisplaced(TD),
        t(pars.t), U(pars.U), mu(pars.mu), L(pars.L), d(pars.d),
        z(2*d), //coordination number: 2*d
        N(static_cast<uint32_t>(uint_pow(L,d))),
//      beta(pars.beta), m(pars.m), s(pars.s), n(m / s), dtau(beta/m),
        alpha(acosh(std::exp(dtau * U * 0.5))),
        neigh(d, L),
        proptmat(N,N),
        auxfield(N, m+1),                //m+1 columns of N rows
        gUp(green[GreenCompSpinUp]), gDn(green[GreenCompSpinDown]),
        gFwdUp(greenFwd[GreenCompSpinUp]), gFwdDn(greenFwd[GreenCompSpinDown]),
        gBwdUp(greenBwd[GreenCompSpinUp]), gBwdDn(greenBwd[GreenCompSpinDown]),
        UdVStorageUp(UdVStorage[GreenCompSpinUp]), UdVStorageDn(UdVStorage[GreenCompSpinDown]),
        occUp(), occDn(), occTotal(), eKinetic(), ePotential(), eTotal(),
        occDouble(), localMoment(), suscq0(), zcorr(), gf(m), gf_dt(m)
{
    gUp = CubeNum(N,N,m+1);
    gDn = CubeNum(N,N,m+1);
    gFwdUp = CubeNum(N,N,m+1);
    gFwdDn = CubeNum(N,N,m+1);
    gBwdUp = CubeNum(N,N,m+1);
    gBwdDn = CubeNum(N,N,m+1);

    setupRandomAuxfield();
    if (CB) {
        setupPropTmat_checkerboard();
    } else {
        setupPropTmat_direct();
    }
    setupUdVStorage_skeleton(hubbardComputeBmat(this));

    lastSweepDir = Base::SweepDirection::Up;        //first sweep will be downwards

    using namespace boost::assign;         // bring operator+=() into scope
    using std::cref;

    obsScalar += ScalarObservable(cref(occUp), "occupationUp", "nUp"),
                 ScalarObservable(cref(occDn), "occupationDown", "nDown"),
                 ScalarObservable(cref(occTotal), "totalOccupation", "n"),
                 ScalarObservable(cref(occDouble), "doubleOccupation", "n2"),
                 ScalarObservable(cref(localMoment), "localMoment", "m^2"),
                 ScalarObservable(cref(eKinetic), "kineticEnergy", "e_t"),
                 ScalarObservable(cref(ePotential), "potentialEnergy", "e_U"),
                 ScalarObservable(cref(eTotal), "totalEnergy", "e");

    zcorr.zeros(N);
    obsVector += VectorObservable(cref(zcorr), N, "spinzCorrelationFunction", "zcorr");

    if (timedisplaced) {
        obsScalar += ScalarObservable(cref(suscq0), "susceptibilityQ0", "chi_q0");
        if (d == 2) {
            for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
                gf_dt[timeslice - 1] = dtau * timeslice;
            }
            gf_dt.reshape(gf_dt.n_elem, 1);     //column vector
            gf.zeros(gf_dt.n_elem);

            obsKeyValue += KeyValueObservable(cref(gf), gf_dt, "dt", "greenFourier", "gf");
        }
    }
}

template <bool TD, bool CB>
DetHubbard<TD,CB>::~DetHubbard() {
}


template <bool TD, bool CB>
MetadataMap DetHubbard<TD,CB>::prepareModelMetadataMap() const {
    MetadataMap meta;
    meta["model"] = "hubbard";
    meta["checkerboard"] = (checkerboard ? "true" : "false");
    meta["timedisplaced"] = (timedisplaced ? "true" : "false");
#define META_INSERT(VAR) {meta[#VAR] = numToString(VAR);}
    META_INSERT(t);
    META_INSERT(U);
    META_INSERT(mu);
    META_INSERT(L);
    META_INSERT(d);
    META_INSERT(N);
    META_INSERT(beta);
    META_INSERT(m);
    META_INSERT(dtau);
    META_INSERT(s);
    META_INSERT(alpha);
#undef META_INSERT
    return meta;
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::updateInSlice(uint32_t timeslice) {
    // picking sites linearly: system seemed to alternate between two configurations
    // sweep after sweep
//  for (uint32_t site = 0; site < N; ++site) {
//  std::cout << timeslice << " ";  //DEBUG
    for (uint32_t count = 0; count < N; ++count) {
        uint32_t site = rng.randInt(0, N-1);

        num ratio =  weightRatioSingleFlip(site, timeslice);

//      //DEBUG: comparison of weight ratio calculation
//      MatInt newAuxfield = auxfield;
//      newAuxfield(site, timeslice) *= -1;
//      num refRatio = weightRatioGenericNaive(auxfield, newAuxfield);
//      std::cout << ratio << " vs. " << refRatio << std::endl;

        assert(ratio > 0.0);

//      std::cout << ratio << '\n';

        //Metropolis
        if (ratio > 1.0 or rng.rand01() < ratio) {
            //          if (refRatio > 1 or rng.rand01() < refRatio) {
            //              std::cout << "acc" << '\n';
            //              auxfield = newAuxfield;

            updateGreenFunctionWithFlip(site, timeslice);
            auxfield(site, timeslice) *= -1;
        }
    }
}

//
//void DetHubbard::sweepSimple() {
//  for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
//      gUp.slice(timeslice) = computeGreenFunctionNaive(timeslice, Spin::Up);
//      gDn.slice(timeslice) = computeGreenFunctionNaive(timeslice, Spin::Down);
//      updateInSlice(timeslice);
//  }
//}

template <bool TD, bool CB>
uint32_t DetHubbard<TD,CB>::getSystemN() const {
    return N;
}

//DetHubbard::MatNum4 DetHubbard::greenFromUdV_timedisplaced(
//      const UdVnum& UdV_l, const UdVnum& UdV_r) const {
//  //Ul vs Vl to be compatible with labeling in the notes
//  const MatNum& Ul = UdV_l.V;   //!
//  const VecNum& dl = UdV_l.d;
//  const MatNum& Vl = UdV_l.U;   //!
//  const MatNum& Ur = UdV_r.U;
//  const VecNum& dr = UdV_r.d;
//  const MatNum& Vr = UdV_r.V;
//
//  //submatrix view helpers for 2*N x 2*N matrices
//  //#define upleft(m) m.submat(0,0, N-1,N-1)
//  //#define upright(m) m.submat(0,N, N-1,2*N-1)
//  //#define downleft(m) m.submat(N,0, 2*N-1,N-1)
//  //#define downright(m) m.submat(N,N, 2*N-1,2*N-1)
//  auto upleft = [N](MatNum& m) {
//      return m.submat(0,0, N-1,N-1);
//  };
//  auto upright = [N](MatNum& m) {
//      return m.submat(0,N, N-1,2*N-1);
//  };
//  auto downleft = [N](MatNum& m) {
//      return m.submat(N,0, 2*N-1,N-1);
//  };
//  auto downright = [N](MatNum& m) {
//      return m.submat(N,N, 2*N-1,2*N-1);
//  };
//
//  MatNum temp(2*N,2*N);
//  upleft(temp)    = arma::inv(Vr * Vl);
//  upright(temp)   = arma::diagmat(dl);
//  downleft(temp)  = arma::diagmat(-dr);
//  downright(temp) = arma::inv(Ul * Ur);
//  UdVnum tempUdV = udvNumDecompose(temp);
//
//  MatNum left(2*N,2*N);
//  upleft(left) = arma::inv(Vr);
//  upright(left).zeros();
//  downleft(left).zeros();
//  downright(left) = arma::inv(Ul);
//
//  MatNum right(2*N,2*N);
//  upleft(right) = arma::inv(Vl);
//  upright(right).zeros();
//  downleft(right).zeros();
//  downright(right) = arma::inv(Ur);
//
//  MatNum result = (left * arma::inv(tempUdV.V)) * arma::diagmat(1.0 / tempUdV.d)
//                  * (arma::inv(tempUdV.U) * right);
//  return MatNum4(upleft(result), upright(result),
//                 downleft(result), downright(result));
//}

//MatNum DetHubbard::greenFromUdV(const UdVnum& UdV_l, const UdVnum& UdV_r) const {
//  //variable names changed according to labeling in notes
//  const MatNum& V_l = UdV_l.U;   //!
//  const VecNum& d_l = UdV_l.d;
//  const MatNum& U_l = UdV_l.V;   //!
//  const MatNum& U_r = UdV_r.U;
//  const VecNum& d_r = UdV_r.d;
//  const MatNum& V_r = UdV_r.V;
//
//  using arma::inv; using arma::diagmat; using arma::eye;
//
//  UdVnum UdV_temp = udvNumDecompose( inv(U_l * U_r) + diagmat(d_r) * (V_r * V_l) * diagmat(d_l) );
//
//  MatNum green = inv(UdV_temp.V * U_l) * diagmat(1.0 / UdV_temp.d) * inv(U_r * UdV_temp.U);
//
//  return green;
//}


//void DetHubbard::setupUdVStorage() {
//  eye_UdV.U = arma::eye(N,N);
//  eye_UdV.d = arma::ones(N);
//  eye_UdV.V = arma::eye(N,N);
//
//  auto setup = [this](std::vector<UdVnum>& storage, Spin spinz) {
//      storage = std::vector<UdVnum>(n + 1);
//
//      storage[0] = eye_UdV;
//      storage[1] = udvNumDecompose(computeBmatFunc(s, 0, spinz));
//
//      for (uint32_t l = 1; l <= n - 1; ++l) {
//          const MatNum& U_l = storage[l].U;
//          const VecNum& d_l = storage[l].d;
//          const MatNum& V_l = storage[l].V;
//          MatNum B_lp1 = computeBmatFunc(s*(l + 1), s*l, spinz);
//          UdVnum UdV_temp = udvNumDecompose((B_lp1 * U_l) * arma::diagmat(d_l));
//          storage[l+1].U = UdV_temp.U;
//          storage[l+1].d = UdV_temp.d;
//          storage[l+1].V = UdV_temp.V * V_l;
//      }
//  };
//
//  setup(UdVStorageUp, Spin::Up);
//  setup(UdVStorageDn, Spin::Down);
//
//  lastSweepDir = SweepDirection::Up;
//}


//void DetHubbard::debugCheckBeforeSweepDown() {
//  std::cout << "Before sweep down:\n";
//  std::cout << "up: ";
//  for (uint32_t l = 0; l <= n; ++l) {
//      UdVnum udv = UdVStorageUp[l];
//      MatNum diff = computeBmatFunc(l*s, 0, Spin::Up) - udv.U * arma::diagmat(udv.d) * udv.V;
//      std::cout << diff.max() << " ";
//  }
//  std::cout << "\n";
//  std::cout << "down: ";
//  for (uint32_t l = 0; l <= n; ++l) {
//      UdVnum udv = UdVStorageDn[l];
//      MatNum diff = computeBmatFunc(l*s, 0, Spin::Down) - udv.U * arma::diagmat(udv.d) * udv.V;
//      std::cout << diff.max() << " ";
//  }
//  std::cout << "\n\n";
//}
//
//void DetHubbard::debugCheckBeforeSweepUp() {
//  std::cout << "Before sweep up:\n";
//  std::cout << "up: ";
//  for (uint32_t l = 0; l <= n; ++l) {
//      UdVnum udv = UdVStorageUp[l];
//      MatNum diff = computeBmatFunc(m, l*s, Spin::Up) - udv.U * arma::diagmat(udv.d) * udv.V;
//      std::cout << diff.max() << " ";
//  }
//  std::cout << "\n";
//  std::cout << "down: ";
//  for (uint32_t l = 0; l <= n; ++l) {
//      UdVnum udv = UdVStorageDn[l];
//      MatNum diff = computeBmatFunc(m, l*s, Spin::Down) - udv.U * arma::diagmat(udv.d) * udv.V;
//      std::cout << diff.max() << " ";
//  }
//  std::cout << "\n\n";
//}
//
//
//
//void DetHubbard::debugCheckGreenFunctions() {
//  std::cout << "debugCheckGreenFunctions:\n";
//  std::cout << "up: ";
//  for (uint32_t k = 1; k <= m; ++k) {
//      MatNum green = gUp.slice(k);
//      MatNum reldiff = (computeGreenFunctionNaive(k, Spin::Up) - green) / green;
//      std::cout << reldiff.max() << " ";
//  }
//  std::cout << "\n";
//  std::cout << "down: ";
//  for (uint32_t k = 1; k <= m; ++k) {
//      MatNum green = gDn.slice(k);
//      MatNum reldiff = (computeGreenFunctionNaive(k, Spin::Down) - green) / green;
//      std::cout << reldiff.max() << " ";
//  }
//  std::cout << "\n\n";
//}



//void DetHubbard::sweep() {
//  using std::tie; using std::ignore; using std::get;
//
//  //compute the green function in timeslice s*(l-1) from scratch with the help
//  //of the B-matrices computed before in the last up-sweep
//  auto advanceDownGreen = [this](uint32_t l, std::vector<UdVnum>& storage,
//          CubeNum& green, CubeNum& greenFwd,
//          CubeNum& greenBwd, Spin spinz) -> void {
//      MatNum B_l = computeBmatFunc(s*l, s*(l - 1), spinz);
//
//      //U_l, d_l, V_l correspond to B(beta,l*s*dtau) [set in the last step]
//      const MatNum& U_l = storage[l].U;
//      const VecNum& d_l = storage[l].d;
//      const MatNum& V_l = storage[l].V;
//
//      //UdV_L will correspond to B(beta,(l-1)*s*dtau)
//      UdVnum UdV_L = udvNumDecompose(arma::diagmat(d_l) * (V_l * B_l));
//      UdV_L.U = U_l * UdV_L.U;
//
//      //UdV_R corresponds to B((l-1)*s*dtau,0) [set in last sweep]
//      const UdVnum& UdV_R = storage[l - 1];
//
//      uint32_t next = s * (l - 1);
//      tie(ignore, greenBwd.slice(next), greenFwd.slice(next), green.slice(next)) =
//              greenFromUdV_timedisplaced(UdV_L, UdV_R);
//      storage[l - 1] = UdV_L;
//  };
//
//  //compute the green function at k-1 by wrapping the one at k (accumulates rounding errors)
//  auto wrapDownGreen = [this](uint32_t k, CubeNum& green, CubeNum& greenFwd,
//          CubeNum& greenBwd, Spin spinz) -> void {
//      MatNum B_k = computeBmatFunc(k, k - 1, spinz);
//      green.slice(k - 1) = arma::inv(B_k) * green.slice(k) * B_k;
//      greenFwd.slice(k - 1) = arma::inv(B_k) * greenFwd.slice(k);
//      greenBwd.slice(k - 1) = greenBwd(k) * B_k;
//  };
//
//  //update the green function in timeslice s*(l+1) from scratch with the help
//  //of B-matrices computed before
//  auto advanceUpGreen = [this](uint32_t l, const std::vector<UdVnum>& storage,
//          CubeNum& green, CubeNum& greenFwd,
//          CubeNum& greenBwd, Spin spinz) -> void {
//      MatNum B_lp1 = computeBmatFunc(s*(l + 1), s*l, spinz);
//
//      //The following is B(beta, (l+1)*s*dtau), valid from the last sweep
//      const UdVnum& UdV_lp1 = storage[l + 1];
//
//      //from the last step the following are B(l*s*dtau, 0):
//      const MatNum& U_l = storage[l].U;
//      const VecNum& d_l = storage[l].d;
//      const MatNum& V_l = storage[l].V;
//
//      //UdV_temp will be the new B((l+1)*s*dtau, 0):
//      UdVnum UdV_temp = udvNumDecompose(((B_lp1 * U_l) * arma::diagmat(d_l)));
//      UdV_temp.V *= V_l;
//
//      uint32_t next = s * (l + 1);
//      tie(ignore, greenBwd.slice(next), greenFwd.slice(next), green.slice(next)) =
//              greenFromUdV_timedisplaced(UdV_lp1, UdV_temp);
//
//      //storage[l + 1] = UdV_temp;    //storage would be wrong after updateInSlice!
//  };
//
//  //Given B(l*s*dtau, 0) from the last step in the storage, compute
//  //B((l+1)*s*dtau, 0) and put it into storage
//  auto advanceUpUpdateStorage = [this](uint32_t l, std::vector<UdVnum>& storage,
//          Spin spinz) -> void {
//      MatNum B_lp1 = computeBmatFunc(s*(l + 1), s*l, spinz);
//      //from the last step the following are B(l*s*dtau, 0):
//      const MatNum& U_l = storage[l].U;
//      const VecNum& d_l = storage[l].d;
//      const MatNum& V_l = storage[l].V;
//      //the new B((l+1)*s*dtau, 0):
//      storage[l+1] = udvNumDecompose(((B_lp1 * U_l) * arma::diagmat(d_l)));
//      storage[l+1].V *= V_l;
//  };
//
//  //compute the green function at k+1 by wrapping the one at k (accumulates rounding errors)
//  auto wrapUpGreen = [this](uint32_t k, CubeNum& green, CubeNum& greenFwd,
//          CubeNum& greenBwd, Spin spinz) -> void {
//      MatNum B_kp1 = computeBmatFunc(k + 1, k, spinz);
//      green.slice(k + 1) = B_kp1 * green.slice(k) * arma::inv(B_kp1);
//      greenFwd.slice(k + 1) = B_kp1 * greenFwd.slice(k);
//      greenBwd.slice(k + 1) = greenBwd.slice(k) * arma::inv(B_kp1);
//  };
//
//  if (lastSweepDir == SweepDirection::Up) {
////        debugCheckBeforeSweepDown();
////        debugCheckGreenFunctions();
//      //to compute green function for timeslice tau=beta:
//      //we need VlDlUl = B(beta, beta) = I and UrDrVr = B(beta, 0).
//      //The latter is given in storage slice m from the last sweep.
//      tie(ignore, gBwdUp.slice(m), gFwdUp.slice(m), gUp.slice(m)) =
//              greenFromUdV_timedisplaced(eye_UdV, UdVStorageUp[n]);
//      tie(ignore, gBwdDn.slice(m), gFwdDn.slice(m), gDn.slice(m)) =
//              greenFromUdV_timedisplaced(eye_UdV, UdVStorageDn[n]);
//      UdVStorageUp[n] = eye_UdV;
//      UdVStorageDn[n] = eye_UdV;
//      for (uint32_t l = n; l >= 1; --l) {
//          updateInSlice(l*s);
//          for (uint32_t k = l*s - 1; k >= (l-1)*s + 1; --k) {
//              wrapDownGreen(k + 1, gUp, gFwdUp, gBwdUp, Spin::Up);
//              wrapDownGreen(k + 1, gDn, gFwdDn, gBwdDn, Spin::Down);
//              updateInSlice(k);
//          }
//          //TODO: this will also compute the Green function at k=0, which technically is not necessary
//          //but sensible for the following sweep up
//          //TODO: alternatively just copy the k=m Green function to k=0  -- would that be up-to-date?
//          advanceDownGreen(l, UdVStorageUp, gUp, gFwdUp, gBwdUp, Spin::Up);
//          advanceDownGreen(l, UdVStorageDn, gDn, gFwdDn, gBwdDn, Spin::Down);
//      }
//      lastSweepDir = SweepDirection::Down;
//  } else if (lastSweepDir == SweepDirection::Down) {
////        debugCheckGreenFunctions();
////        debugCheckBeforeSweepUp();
//      //We need to have computed the Green function for time slice k=0 so that the first
//      //wrap-up step is correct.
//      for (uint32_t k = 1; k <= s-1; ++k) {
//          wrapUpGreen(k - 1, gUp, gFwdUp, gBwdUp, Spin::Up);
//          wrapUpGreen(k - 1, gDn, gFwdDn, gBwdDn, Spin::Down);
//          updateInSlice(k);
//      }
//      //set storage at k=0 to unity for the upcoming sweep:
//      UdVStorageUp[0] = eye_UdV;
//      UdVStorageDn[0] = eye_UdV;
//      for (uint32_t l = 1; l < n; ++l) {
//          advanceUpGreen(l-1, UdVStorageUp, gUp, gFwdUp, gBwdUp, Spin::Up);
//          advanceUpGreen(l-1, UdVStorageDn, gDn, gFwdDn, gBwdDn, Spin::Down);
//          updateInSlice(l*s);
//          advanceUpUpdateStorage(l - 1, UdVStorageUp, Spin::Up);
//          advanceUpUpdateStorage(l - 1, UdVStorageDn, Spin::Down);
//          for (uint32_t k = l*s + 1; k <= l*s + (s-1); ++k) {
//              wrapUpGreen(k - 1, gUp, gFwdUp, gBwdUp, Spin::Up);
//              wrapUpGreen(k - 1, gDn, gFwdDn, gBwdDn, Spin::Down);
//              updateInSlice(k);
//          }
//      }
//      updateInSlice(n*s);
//      advanceUpUpdateStorage(n - 1, UdVStorageUp, Spin::Up);
//      advanceUpUpdateStorage(n - 1, UdVStorageDn, Spin::Down);
//      lastSweepDir = SweepDirection::Up;
//  }
////    std::cout << std::endl;     //DEBUG
//}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::measure() {
    //used to measure occupation:
    num sum_GiiUp = 0;
    num sum_GiiDn = 0;
    //used to measure kinetic energy:
    num sum_GneighUp = 0;
    num sum_GneighDn = 0;
//  //used to measure double occupancy / potential energy:
    num sum_GiiUpDn = 0;
    //FORMULA-TEST -- made no difference
//  num sum_doubleoccupancy = 0;
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for (uint32_t site = 0; site < N; ++site) {
            //use diagonal elements of Green functions:
            sum_GiiUp += gUp(site, site, timeslice);
            sum_GiiDn += gDn(site, site, timeslice);
            sum_GiiUpDn += gUp(site, site, timeslice) * gDn(site, site, timeslice);
            //use nearest neighbor elements of Green functions:
//          for (uint32_t neighIndex = 0; neighIndex < z; ++neighIndex) {
            for (auto p = neigh.beginNeighbors(site); p != neigh.endNeighbors(site); ++p) {
                uint32_t neigh = *p;
                sum_GneighUp += gUp(site, neigh, timeslice);
                sum_GneighDn += gDn(site, neigh, timeslice);
            }
            //FORMULA-TEST -- made no difference
//          sum_doubleoccupancy += (1 - gUp(site,site, timeslice))
//                               * (1 - gDn(site,site, timeslice));
        }
    }
    occUp = 1.0 - (1.0 / (N*m)) * sum_GiiUp;
    occDn = 1.0 - (1.0 / (N*m)) * sum_GiiDn;
    occTotal = occUp + occDn;

    //FORMULA-TEST -- made no difference
//  std::cout << (1.0 / (N*m)) * sum_doubleoccupancy;
    occDouble = 1.0 + (1.0 / (N*m)) * (sum_GiiUpDn - sum_GiiUp - sum_GiiDn);
//  std::cout << " vs. " << occDouble << std::endl;

    localMoment = occTotal - 2*occDouble;

//  ePotential = (U / (N*m)) * (sum_GiiUpDn + 0.5 * sum_GiiUp + 0.5 * sum_GiiDn);
//  ePotential = U * ( 0.25 + (1.0 / (N*m)) * (sum_GiiUpDn - 0.5 * (sum_GiiUp + sum_GiiDn)) );
    ePotential = U * occDouble;

    //Note: chemical potential term included in kinetic energy:
    eKinetic   = (t / (N*m)) * (sum_GneighUp + sum_GneighDn) - mu * occTotal;
    eTotal = eKinetic + ePotential;

    //susceptibility
    if (timedisplaced) {
        uint32_t mm = m;        // I don't understand why m can't be captured for the lambda without this line
        auto sumTrace = [this, mm](const CubeNum& green) -> num {
            num sum = 0;
            for (uint32_t timeslice = 1; timeslice <= mm; ++timeslice) {
                sum += arma::trace(green.slice(timeslice));
            }
            return sum;
        };
        num sumTrGreenUp = sumTrace(gUp);
        num sumTrGreenDn = sumTrace(gDn);
        auto sumProdTrace = [this, mm](const CubeNum& green1, const CubeNum& green2) -> num {
            num sum = 0;
            for (uint32_t timeslice = 1; timeslice <= mm; ++timeslice) {
                sum += arma::trace(green1.slice(timeslice) * green2.slice(timeslice));
            }
            return sum;
        };
        num sumTrGreenDisplacedUp = sumProdTrace(gBwdUp, gFwdUp);
        num sumTrGreenDisplacedDn = sumProdTrace(gBwdDn, gFwdDn);
        num trGreenUp_0 = arma::trace(gUp.slice(m));            //g(beta) = g(0)
        num trGreenDn_0 = arma::trace(gDn.slice(m));
        suscq0 = (1.0 / num(N)) * dtau * ( (trGreenUp_0 - trGreenDn_0) * (sumTrGreenUp - sumTrGreenDn)
                - (sumTrGreenDisplacedUp + sumTrGreenDisplacedDn)
        );
    }

    // vector observables
    zcorr.zeros();
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        num gUp_00 = gUp(0,0, timeslice);
        num gDn_00 = gDn(0,0, timeslice);
        zcorr[0] += -2.0 * gUp_00 * gDn_00 + gUp_00 + gDn_00;
        for (uint32_t siteJ = 1; siteJ < N; ++siteJ) {
            num gUp_0j = gUp(0,siteJ,     timeslice);
            num gDn_0j = gDn(0,siteJ,     timeslice);
            num gUp_jj = gUp(siteJ,siteJ, timeslice);
            num gDn_jj = gDn(siteJ,siteJ, timeslice);
            using std::pow;
            zcorr[siteJ] += gUp_00 * gUp_jj - gUp_00 * gDn_jj + gDn_00 * gDn_jj - gDn_00 * gUp_jj
                    - pow(gUp_0j, 2) - pow(gDn_0j, 2);
        }
    }
    zcorr /= num(m);

    if (timedisplaced and d == 2) {
        //compute the Fourier transform of the imaginary time-displaced forward Green function for
        //k = (pi/3, 2*pi/3) at some preset values of tau
        const num kx = M_PI / 3.0;
        const num ky = 2.0 * M_PI / 3.0;
        for (uint32_t timeslice = 1; timeslice < m; ++timeslice) {
//          std::complex<num> sum(0,0);
            num sum = 0.0;
            for (uint32_t siteI = 0; siteI < N; ++siteI) {
                const num yI = siteI / L;
                const num xI = siteI % L;
                for (uint32_t siteJ = 0; siteJ < N; ++siteJ) {
                    const num yJ = siteJ / L;
                    const num xJ = siteJ % L;
//                  const std::complex<num> expFactor =
//                          std::exp(std::complex<num>(0, kx * (xJ - xI) + ky * (yJ - yI)));
                    const num cosFactor = std::cos(kx * (xJ - xI) + ky * (yJ - yI));

//                  sum += expFactor * gFwdUp(siteI, siteJ, timeslice);
//                  sum += expFactor * gFwdDn(siteI, siteJ, timeslice);
                    sum += cosFactor * gFwdUp(siteI, siteJ, timeslice);
                    sum += cosFactor * gFwdDn(siteI, siteJ, timeslice);
                }
            }
            sum /= num(N);
            gf[timeslice - 1] = sum;
        }
    }
}


template <bool TD, bool CB>
void DetHubbard<TD,CB>::setupRandomAuxfield() {
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for (uint32_t site = 0; site < N; ++site) {
            if (rng.rand01() <= 0.5) {
                auxfield(site, timeslice) = +1;
            } else {
                auxfield(site, timeslice) = -1;
            }
        }
    }
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::setupPropTmat_direct() {
    MatNum tmat = -mu * arma::eye(N, N);

    for (uint32_t site = 0; site < N; ++site) {
        //hopping between nearest neighbors
        for (auto p = neigh.beginNeighbors(site);
                 p != neigh.endNeighbors(site); ++p) {
            tmat(*p, site) -= t;
        }
    }

    proptmat = computePropagator(dtau, tmat);
}


template <bool TD, bool CB>
void DetHubbard<TD,CB>::setupPropTmat_checkerboard() {
    //Checkerboard break up as in dos Santos 2003

    assert(d == 2);
    auto xyToSite = [this](uint32_t x, uint32_t y) {
        return y * L + x;
    };
    SpMatNum kxa = SpMatNum(N, N);
    SpMatNum kxb = SpMatNum(N, N);
    SpMatNum kya = SpMatNum(N, N);
    SpMatNum kyb = SpMatNum(N, N);
    for (uint32_t y = 0; y < L; ++y) {
        for (uint32_t x = 0; x < L; x += 2) {
            //sub board a
            uint32_t siteA = xyToSite(x, y);
            uint32_t neighA = neigh(NeighDir::XPLUS, siteA);
            kxa(siteA, neighA) = 1.0;
            kxa(neighA, siteA) = 1.0;
            //sub board b
            uint32_t siteB = neighA;
            uint32_t neighB = neigh(NeighDir::XPLUS, siteB);
            kxb(siteB, neighB) = 1.0;
            kxb(neighB, siteB) = 1.0;
        }
    }
    for (uint32_t x = 0; x < L; ++x) {
        for (uint32_t y = 0; y < L; y += 2) {
            //sub board a
            uint32_t siteA = xyToSite(x, y);
            uint32_t neighA = neigh(NeighDir::YPLUS, siteA);
            kya(siteA, neighA) = 1.0;
            kya(neighA, siteA) = 1.0;
            //sub board b
            uint32_t siteB = neighA;
            uint32_t neighB = neigh(NeighDir::YPLUS, siteB);
            kyb(siteB, neighB) = 1.0;
            kyb(neighB, siteB) = 1.0;
        }
    }
    using std::cosh; using std::sinh; using std::pow;
    const num ch = cosh(dtau * t);
    const num sh = sinh(dtau * t);
    using arma::eye; using arma::conv_to;
    //see notes / Mathematica experimentation:
    //[Actually, by performing the calculations in this way, there is no numerical
    //advantage whatsoever -- more sensible implementation to be done for the SDW model]x
    proptmat = pow(ch, 4) * eye(N,N)
             + pow(ch, 3)*sh * (kxa + kxb + kya + kyb)
             + pow(ch, 2)*pow(sh, 2 ) * (kxa*kxb + kxa*kya + kxb*kya +
                                          kxa*kyb + kxb*kyb + kya*kyb)
             + ch*pow(sh, 3) * (kxa*kxb*kya + kxa*kxb*kyb + kxa*kya*kyb +
                                kxb*kya*kyb)
             + pow(sh, 4) * kxa*kxb*kya*kyb;
}

template <bool TD, bool CB>
inline MatNum DetHubbard<TD,CB>::computeBmat(uint32_t k2, uint32_t k1, Spin spinz) const {
    using namespace arma;

    if (k2 == k1) {
        return eye(N, N);
    }

    assert(k2 > k1);
    assert(k2 <= m);
    //assert(n1 >= 0);

    num sign = num(int(spinz));

    //Propagator using the HS-field potential for the given timeslice
    auto singleTimeslicePropagator = [this, sign](uint32_t timeslice) -> MatNum {
        //the cast with conv_to is necessary here, else everything would result in integers!
        //-- an Armadillo bug IMHO
        return diagmat(exp(sign * alpha *
                conv_to<VecNum>::from(auxfield.col(timeslice)))) * proptmat;
    };

    MatNum B = singleTimeslicePropagator(k2);

    for (uint32_t k = k2 - 1; k >= k1 + 1; --k) {
        B *= singleTimeslicePropagator(k);
    }

    return B;
}

//template <bool TD, bool CB>
//MatNum DetHubbard<TD,CB>::computeBmat_checkerBoard(uint32_t k2, uint32_t k1,
//      Spin spinz) const {
//  //for now: in the checkerboard decomposition we generate the same type
//  //of proptmat as with the direct calculation
//  return computeBmat_direct(k2, k1, spinz);
//}
//
//
//inline MatNum DetHubbard::computeGreenFunctionNaive(
//      const MatNum& bTau0, const MatNum& bBetaTau) const {
//  return arma::inv(arma::eye(N,N) + bTau0 * bBetaTau);
//}
//
//inline MatNum DetHubbard::computeGreenFunctionNaive(uint32_t timeslice,
//      Spin spinz) const {
//  //TODO: should use stored B-matrices, for the timeslices that have not changed
//  return computeGreenFunctionNaive(computeBmat_direct(timeslice, 0, spinz),
//                              computeBmat_direct(m, timeslice, spinz));
//}


//num DetHubbard::weightRatioGenericNaive(const MatInt& auxfieldBefore,
//      const MatInt& auxfieldAfter) const {
//  using namespace arma;
//
//  num weightAfterUp  = det(eye(N,N) + computeBmat_direct(m, 0, Spin::Up, auxfieldAfter));
//  num weightBeforeUp = det(eye(N,N) + computeBmat_direct(m, 0, Spin::Up, auxfieldBefore));
//  num ratioUp = weightAfterUp / weightBeforeUp;
//
//  num weightAfterDown  = det(eye(N,N) + computeBmat_direct(m, 0, Spin::Down, auxfieldAfter));
//  num weightBeforeDown = det(eye(N,N) + computeBmat_direct(m, 0, Spin::Down, auxfieldBefore));
//  num ratioDown = weightAfterDown / weightBeforeDown;
//
////    std::cout << weightafterup << " " << weightbeforeup << " " << weightafterdown << " " << weightbeforedown << "\n";
//
//  return ratioUp * ratioDown;
//}

template <bool TD, bool CB>
inline num DetHubbard<TD,CB>::weightRatioSingleFlip(uint32_t site, uint32_t timeslice) const {
    using std::exp;
    //TODO: possibly precompute the exponential factors (auxfield is either +/- 1), would require an if though.

    //exponential factors
    //results again do not seem to for the location of the -sign
    num expUp   = exp(-2.0 * alpha * num(auxfield(site, timeslice)));
    num expDown = exp( 2.0 * alpha * num(auxfield(site, timeslice)));

    num ratioUp   = 1.0 + (expUp   - 1.0) * (1.0 - gUp(site,site, timeslice));
    num ratioDown = 1.0 + (expDown - 1.0) * (1.0 - gDn(site,site, timeslice));

//  std::cout << ratioUp << " " << ratioDown << std::endl;

    return ratioUp * ratioDown;
}

template <bool TD, bool CB>
inline void DetHubbard<TD,CB>::updateGreenFunctionWithFlip(uint32_t site, uint32_t timeslice) {
    auto update = [this, site](MatNum& green, num deltaSite) {
        const MatNum& greenOld = green;     //reference
        MatNum greenNew = green;            //copy
        const MatNum oneMinusGreenOld = arma::eye(N,N) - greenOld;
        num divisor = 1.0 + deltaSite * oneMinusGreenOld(site, site);
        num greenFactor = deltaSite / divisor;
        for (uint32_t y = 0; y < N; ++y) {
            for (uint32_t x = 0; x < N; ++x) {
                greenNew(x, y) -=  greenOld(x, site) * greenFactor *
                        oneMinusGreenOld(site, y);

                //experimental index swap below -- this seemed to be the choice in Santos 2003
//              greenNew(x, y) -=  greenOld(site, y) * greenFactor *
//                      oneMinusGreenOld(x, site);
            }
        }
        green = greenNew;
    };

    using std::exp;
    update(gUp.slice(timeslice), exp(-2.0 * alpha * num(auxfield(site, timeslice))) - 1.0);
    update(gDn.slice(timeslice), exp(+2.0 * alpha * num(auxfield(site, timeslice))) - 1.0);
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::sweepSimple() {
    sweepSimple_skeleton(hubbardComputeBmat(this));
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::sweepSimpleThermalization() {
    sweepSimpleThermalization_skeleton(hubbardComputeBmat(this));
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::sweep() {
    sweep_skeleton(hubbardLeftMultiplyBmat(this), hubbardRightMultiplyBmat(this),
                   hubbardLeftMultiplyBmatInv(this), hubbardRightMultiplyBmatInv(this));
}

template <bool TD, bool CB>
void DetHubbard<TD,CB>::sweepThermalization() {
    sweepThermalization_skeleton(hubbardLeftMultiplyBmat(this), hubbardRightMultiplyBmat(this),
                                 hubbardLeftMultiplyBmatInv(this), hubbardRightMultiplyBmatInv(this));
}

