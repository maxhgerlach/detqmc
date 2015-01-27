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

#include "exceptions.h"
#include "detmodel.h"
#include "neighbortable.h"
#include "rngwrapper.h"
#include "dethubbardparams.h"
#include "metadata.h"
#include "observable.h"
#include "udv.h"

// factory function to init DetHubbard from parameter struct
//
// will do parameter checking etc
void createReplica(std::unique_ptr<DetHubbard>& replica_out, RngWrapper& rng, ModelParams<DetHubbard> pars,
                   DetModelLoggingParams logpar = DetModelLoggingParams() /*ignore this one*/);



class DetHubbard : public DetModelGC<2, num, false> {
private:
    //only initialize with the "factory" function ::createReplica() declared above.
    //Give a reference to the RNG instance to be used
    DetHubbard(RngWrapper& rng, const ModelParams<DetHubbard>& pars);
public:
    friend void createReplica(std::unique_ptr<DetHubbard>& replica_out, RngWrapper& rng, ModelParams<DetHubbard> pars,
                              DetModelLoggingParams /*ignored*/);
    virtual ~DetHubbard();

    virtual uint32_t getSystemN() const;

    //Create a MetadataMap describing the parameters of the
    //simulated model
    virtual MetadataMap prepareModelMetadataMap() const;

    //    //perform measurements of all observables
    //    virtual void measure();

    virtual void sweep(bool takeMeasurements);
    virtual void sweepThermalization();
    virtual void sweepSimple(bool takeMeasurements);
    virtual void sweepSimpleThermalization();

    // These currently are just dummy functions. DetHubbard does not
    // support writing system configurations to disk. [and should fail if
    // that program option is invoked]
    void saveConfigurationStreamText(const std::string& directory = ".") {
        (void)directory;
        throw GeneralError("DetHubbard::saveConfigurationStreamText not implemented");
    }
    void saveConfigurationStreamBinary(const std::string& directory = ".") {
        (void)directory;
        throw GeneralError("DetHubbard::saveConfigurationStreamBinary not implemented");
    }
    void saveConfigurationStreamTextHeader(const std::string& simInfoHeaderText,
                                           const std::string& directory = ".") {
        (void)simInfoHeaderText; (void)directory;
        throw GeneralError("DetHubbard::saveConfigurationStreamTextHeader not implemented");
    }
    void saveConfigurationStreamBinaryHeaderfile(const std::string& simInfoHeaderText,
                                                 const std::string& directory = ".") {
        (void)simInfoHeaderText; (void)directory;
        throw GeneralError("DetHubbard::saveConfigurationStreamBinaryHeaderfile not implemented");
    }
protected:
    typedef DetModelGC<2, num, false> Base;
    // stupid C++ weirdness forces us to explicitly "import" these protected base
    // class member variables:
    // (see: http://stackoverflow.com/questions/11405/gcc-problem-using-a-member-of-a-base-class-that-depends-on-a-template-argument )
    using Base::dtau;
    using Base::m;
    using Base::green;
    // using Base::greenFwd;
    // using Base::greenBwd;
    using Base::UdVStorage;
    using Base::lastSweepDir;
    using Base::obsScalar;
    using Base::obsVector;
    using Base::obsKeyValue;
    using Base::beta;
    using Base::s;

    enum class Spin: int {Up = +1, Down = -1};
    enum {GreenCompSpinUp = 0, GreenCompSpinDown = 1};
    RngWrapper& rng;
    //parameters:
    const bool checkerboard;
    const bool timedisplaced;
    const num t;            //hopping energy scale
    const num U;            //interaction energy scale
    const num mu;           //chemical potential
    const uint32_t L;   //linear lattice size
    const uint32_t d;   //spatial dimension of lattice
    const uint32_t z;   // lattice coordination number, 2 * d
    const uint32_t N;   // L ** d
//  const num beta;     //inverse temperature
//  const uint32_t m;   //number of imaginary time discretization steps (time slices) beta*m=dtau
//  const uint32_t s;   //interval between time slices where the Green-function is calculated from scratch
//  const uint32_t n;   //number of time slices where the Green-function is calculated from scratch n*s*dtau=beta
//  const num dtau;     // beta / m
    const num alpha;    // cosh(alpha) = exp(dtau U / 2)

    PeriodicCubicLatticeNearestNeighbors neigh;

    //std::function<MatNum(uint32_t k2, uint32_t k1, Spin spinz)> computeBmatFunc;

    //Matrix representing the kinetic energy part of the hamiltonian: H_t
    //for spin up or spin down -- tmat
    // related propagator e ** (-dtau * tmat):
    MatNum proptmat;


    //the following quantities vary during the course of the simulation


    //Auxiliary field represented by Ising spins, N entries of +/- 1.
    //There is one auxiliary field for each imaginary time slice. One time slice
    //corresponds to one column of the matrix. The time slices in auxfield are
    //indexed from 0 to m. So auxfield.col(n) refers to the timeslice dtau*n.
    //Most code, however, only uses timeslices >= 1 ! Don't rely on auxfield.col(0)
    //being valid.
    MatInt auxfield;

    //Equal imaginary time Green function
    //One cube for each value of spinz.
    MatNum& gUp;
    MatNum& gDn;
    
    //Imaginary time displaced Green function
    // "forward" corresponds to G(tau, 0)
    // "backward" corresponds to G(0, tau)
    //the indexing works the same way as for the equal time case
    // MatNum& gFwdUp;
    // MatNum& gFwdDn;
    // MatNum& gBwdUp;
    // MatNum& gBwdDn;

//  UdVnum eye_UdV; // U = d = V = 1
    std::vector<UdVnum>& UdVStorageUp;
    std::vector<UdVnum>& UdVStorageDn;

    //helper variables to compute observables
    //used to measure occupation:
    num sum_GiiUp;
    num sum_GiiDn;
    //used to measure kinetic energy:
    num sum_GneighUp;
    num sum_GneighDn;
    //used to measure double occupancy / potential energy:
    num sum_GiiUpDn;

    //observables, values for the current auxiliary field; averaged over aux. field
    num occUp;          //occupation spin up
    num occDn;          //occupation spin down
    num occTotal;       //total occupation
    num eKinetic;       //energy -t \sum_<i,j>,\sigma c^+_i,\sigma c_j,\sigma
                        //       -mu \sum_j,\sigma n_j\sigma
    num ePotential;     //energy U \sum_i (n_i,up - 0.5) (n_i,down - 0.5)
    num eTotal;         //total energy
    num occDouble;      //double occupation < nUp * nDown >
    num localMoment;    //local Moment: <m^2> = <(nUp - nDown)^2>
    num suscq0;         //q=0 susceptibility of z component of magnetization (nUp - nDown)

    //the same for vector observables, averaged over timeslices with the current auxiliary field
    VecNum zcorr;       //correlation function of magnetization density at site 0 with all other sites, z component

    //for d=2: the Fourier transform of the timedisplaced Green function for k=(pi/2, 2*pi/3)
    VecNum gf;          //values for different time-displacements
    VecNum gf_dt;       //time-displacements, where the values are evaluated

    void setupRandomAuxfield();
    void setupPropTmat_direct();
    void setupPropTmat_checkerboard();

    //given the current auxiliary fields {s_n}, compute the matrix
    // B_{s_n}(tau_2, tau_1) = \prod_{n = n2}^{n = n1 + 1} e^V(s_n) e^{-dtau T}
    // n2 = tau_2 / m > tau_1 / m = n1
    //So n1, n2 run from 0 to m.
    // e^{-dtau T} = proptmat
    //Here the V(s_n) are computed for either the up or down Hubbard spins.
    //These functions naively multiply the matrices, which can be unstable.
    MatNum computeBmat(uint32_t k2, uint32_t k1, Spin spinz) const;

    //compute the latter using a checker board decomposition with systematic
    //error of O[dtau^2]
    // -- for now this is just the same code
//  MatNum computeBmat_checkerBoard(uint32_t k2, uint32_t k1, Spin spinz) const;

//  //Calculate (1 + B_s(tau, 0)*B_s(beta, tau))^(-1) from the given matrices
//  //for the current aux field.
//  //These functions perform a naive matrix product and inversion.
//  MatNum computeGreenFunctionNaive(const MatNum& bTau0, const MatNum& bBetaTau) const;
//  //Calculate the Green function from scratch for the given timeslice index
//  //(tau = dtau * timeslice) for spin up or down
//  MatNum computeGreenFunctionNaive(uint32_t timeslice, Spin spinz) const;


    //calculate det[1 + B_after(beta, 0)] / det[1 + B_before(beta,0)]
    //by brute force and naively computed B-matrices
//  num weightRatioGenericNaive(const MatInt& auxfieldBefore,
//          const MatInt& auxfieldAfter) const;

    //ratio of weighting determinants if a single auxiliary field
    //spin (at site in timeslice) is flipped from the current configuration.
    //use pre-stored Green functions.
    //formula independent of system size or number of time slices, in this form
    //specific to the Hubbard model
    num weightRatioSingleFlip(uint32_t site, uint32_t timeslice) const;


    //Update the stored Green function matrices to reflect the state after
    //the auxiliary field spin at site in timeslice has been flipped. This
    //function expects this->auxfield to be in the state before the flip.
    void updateGreenFunctionWithFlip(uint32_t site, uint32_t timeslice);

    //update the HS auxiliary field and the green function in the single timeslice
    void updateInSlice(uint32_t timeslice);

    //measuring observables
    void initMeasurements();				//reset stored observable values (beginning of a sweep)
    void measure(uint32_t timeslice);		//measure observables for one timeslice
    void finishMeasurements();				//finalize stored observable values (end of a sweep)

    virtual void consistencyCheck();

    
    //Given B(beta, tau) = V_l d_l U_l and B(tau, 0) = U_r d_r V_r
    //calculate a tuple of four NxN matrices (a,b,c,d) with
    // a = G(0), b = -(1-G(0))*B^(-1)(tau,0), c = B(tau,0)*G(0), d = G(tau)
    //b is the backward time-displaced Green function; c the forward time-
    //displaced Green function; d is the equal-time Green function
//  typedef std::tuple<MatNum,MatNum,MatNum,MatNum> MatNum4;
//  MatNum4 greenFromUdV_timedisplaced(const UdVnum& UdV_l, const UdVnum& UdV_r) const;

    //use a faster method that does not yield information about the time-displaced
    //Green functions
//  MatNum greenFromUdV(const UdVnum& UdV_l, const UdVnum& UdV_r) const;

//  void debugCheckBeforeSweepDown();
//  void debugCheckBeforeSweepUp();
//  void debugCheckGreenFunctions();

    //wrappers to use to instantiate template functions of the base class
    struct hubbardComputeBmat {
        DetHubbard* parent;
        hubbardComputeBmat(DetHubbard* parent_) :
            parent(parent_)
        { }
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
        MatNum operator()(uint32_t gc, uint32_t k2, uint32_t k1) {
            assert(gc == GreenCompSpinUp or gc == GreenCompSpinDown);
            if (gc == GreenCompSpinUp) {
                return parent->computeBmat(k2, k1, Spin::Up);
            } else if (gc == GreenCompSpinDown) {
                return parent->computeBmat(k2, k1, Spin::Down);
            }
        }
#pragma GCC diagnostic pop
    };

    struct hubbardLeftMultiplyBmat {
        DetHubbard* parent;
        hubbardLeftMultiplyBmat(DetHubbard* parent_) :
            parent(parent_)
        { }
        MatNum operator()(uint32_t gc, const MatNum& mat, uint32_t k2, uint32_t k1) {
            return hubbardComputeBmat(parent)(gc, k2, k1) * mat;
        }
    };

    struct hubbardRightMultiplyBmat {
        DetHubbard* parent;
        hubbardRightMultiplyBmat(DetHubbard* parent_) :
            parent(parent_)
        { }
        MatNum operator()(uint32_t gc, const MatNum& mat, uint32_t k2, uint32_t k1) {
            return mat * hubbardComputeBmat(parent)(gc, k2, k1);
        }
    };

    struct hubbardLeftMultiplyBmatInv {
        DetHubbard* parent;
        hubbardLeftMultiplyBmatInv(DetHubbard* parent_) :
            parent(parent_)
        { }
        MatNum operator()(uint32_t gc, const MatNum& mat, uint32_t k2, uint32_t k1) {
            return arma::inv(hubbardComputeBmat(parent)(gc, k2, k1)) * mat;
        }
    };

    struct hubbardRightMultiplyBmatInv {
        DetHubbard* parent;
        hubbardRightMultiplyBmatInv(DetHubbard* parent_) :
            parent(parent_)
        { }
        MatNum operator()(uint32_t gc, const MatNum& mat, uint32_t k2, uint32_t k1) {
            return mat * arma::inv(hubbardComputeBmat(parent)(gc, k2, k1));
        }
    };

public:
    // serialization by selected DetQMC methods
    template<class Archive>
    void saveContents(Archive &ar) {
        Base::saveContents(ar);      //base class
        serializeContentsCommon(ar);
    }

    //after loadContents() a sweep must be performed before any measurements are taken:
    //else the green function would not be in a valid state
    template<class Archive>
    void loadContents(Archive &ar) {
        Base::loadContents(ar);      //base class
        serializeContentsCommon(ar);
        //the fields now have a valid state, update UdV-storage to start
        //sweeping again
        //setupUdVStorage_and_calculateGreen_skeleton(hubbardComputeBmat(this));
        setupUdVStorage_and_calculateGreen_skeleton(hubbardLeftMultiplyBmat(this));      
        //now: lastSweepDir == SweepDirection::Up --> the next sweep will be downwards
    }

    template<class Archive>
    void serializeContentsCommon(Archive &ar) {
        ar & auxfield;
        ar & occUp & occDn & occTotal & eKinetic & ePotential & eTotal
            & occDouble & localMoment & suscq0;
        ar & zcorr & gf & gf_dt;
    }
    
};

#endif /* DETHUBBARD_H_ */
