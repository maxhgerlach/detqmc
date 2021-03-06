/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

/*
 * detsdwopdim.h
 *
 *  Created on: Aug 04, 2014
 *      Author: gerlach
 */

#ifndef DETSDWOPDIM_H_
#define DETSDWOPDIM_H_

#include <tuple>
#include <vector>
#include <complex>
#include <type_traits>          // std::conditional
#include <string>
#include "rngwrapper.h"
#include "detmodel.h"
#include "detsdwparams.h"
#include "neighbortable.h"
#include "RunningAverage.h"
#include "checkarray.h"
#include "normaldistribution.h"
#include "symmat.h"
#include "detsdwsystemconfig.h"
#include "detsdwsystemconfigfilehandle.h"

typedef std::complex<num> cpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<cpx> SpMatCpx;
typedef arma::Col<cpx> VecCpx;
typedef arma::Cube<cpx> CubeCpx;

std::string cbmToString(CheckerboardMethod cbm);


template<CheckerboardMethod CBM, int OPDIM> class DetSDW;

template<CheckerboardMethod CBM, int OPDIM>
void createReplica(std::unique_ptr<DetSDW<CBM, OPDIM>>& replica_out,
                   RngWrapper& rng, ModelParamsDetSDW pars,
                   DetModelLoggingParams loggingPars = DetModelLoggingParams(),
                   // optionally allow passing a directory, where this
                   // replica is allowed to write output into a logfile
                   const std::string& logfiledir = "");


// template parameters:
// Checkerboard:
//   do a checker-board decomposition of which kind?
// OPDIM:
//   - values of 1, 2, 3 to be supported
//   - dimension of the order parameter field phi
template<CheckerboardMethod Checkerboard, int OPDIM>
// class DetSDW: public DetModelGC<1,
//                                 // the basic data type is complex for
//                                 // O(3) or O(2) order parameters, real
//                                 // for O(1):
//                                 typename std::conditional<OPDIM==1, num, cpx>::type > {
class DetSDW: public DetModelGC<1,
                                // to support an extra magnetic field, we use
                                // a complex basic data type for the O(1), O(2),
                                // and O(3) models. ~~ no good way around this, unfortunately.
                                cpx> {
public:
    static_assert(OPDIM==1 or OPDIM==2 or OPDIM==3,
                  "Supported order parameter dimensions: 1, 2, or 3");
    typedef ModelParamsDetSDW ModelParams;
private:
    DetSDW(RngWrapper& rng, const ModelParams& pars,
           const DetModelLoggingParams& loggingPars = DetModelLoggingParams(),
           const std::string& logfiledir = "");
public:
    template<CheckerboardMethod CBM, int OrderParameterDimension>
    friend void createReplica(std::unique_ptr<DetSDW<CBM, OrderParameterDimension>>& replica_out,
                              RngWrapper& rng, ModelParams pars,
                              DetModelLoggingParams loggingPars,
                              const std::string& logfiledir);

    virtual ~DetSDW();
    virtual uint32_t getSystemN() const;
    virtual MetadataMap prepareModelMetadataMap() const; // also contains dynamic data

    virtual void thermalizationOver();
    //the following overloaded version is for use in an MPI
    //context. It writes to stdout and prepends the output
    //by the processIndex
    virtual void thermalizationOver(int processIndex);

    virtual void sweep(bool takeMeasurements);
    virtual void sweepThermalization();
    virtual void sweepSimple(bool takeMeasurements);
    virtual void sweepSimpleThermalization();

    //Write out current system configuration samples to disk: ASCII or
    //binary. Proper filenames detected set up automatically.
    //----------------------------------------------------------------
    void saveConfigurationStreamText(const std::string& directory = ".");
    void saveConfigurationStreamBinary(const std::string& directory = ".");
    // If the file does not already exist, write an informative human
    // readable header. For binary: write it to a separate text file.
    // Only write this file if it does not exist already.
    void saveConfigurationStreamTextHeader(const std::string& simInfoHeaderText,
                                           const std::string& directory = ".");
    void saveConfigurationStreamBinaryHeaderfile(const std::string& simInfoHeaderText,
                                                 const std::string& directory = ".");    

    //Methods used to implement writing out system configurations in a
    //replica exchange simulation (needed by DetQMCPT).
    //----------------------------------------------------------------
    typedef DetSDW_SystemConfig SystemConfig;
    typedef DetSDW_SystemConfig_FileHandle SystemConfig_FileHandle;
    DetSDW_SystemConfig getCurrentSystemConfiguration();
    // the file handle created by the following method will ultimately
    // be used by other replicas simulated on other processors if this
    // is used by DetQMCPT
    DetSDW_SystemConfig_FileHandle prepareSystemConfigurationStreamFileHandle(
        bool binaryStream, bool textStream,
        const std::string& directory = ".");

    
    //Methods to implement a replica-exchange / parallel tempering scheme
    //-------------------------------------------------------------------
    // get / set pars.r
    num  get_exchange_parameter_value() const;
    void set_exchange_parameter_value(num r);
    const char* get_exchange_parameter_name() const;
    // return 1/2 \int_0^\beta d\tau \sum_i [\vec{\phi}_i(\tau)]^2
    num get_exchange_action_contribution() const;
    // get/set control parameter [r] specific data -- here:
    // MC stepsize adjustment done during thermalization.  This is for
    // Boost MPI
    //// and therefore uses this plain pointer to a buffer for
    //// simplicity
    // constexpr uint32_t get_control_data_buffer_size() const {
    //     return ad.bufsize;
    // }
    // void get_control_data(double* buffer) const {
    //     ad.set_buffer_to_data(buffer);
    // }
    // void set_control_data(const double* buffer) {
    //     ad.get_data_from_buffer(buffer);
    // }
    void get_control_data(std::string& buffer) const;
    void set_control_data(const std::string& buffer);
protected:
/*
  The Green's function etc. are 4*N x 4*N matrices for the O(3) model,
  2*N x 2*N for the O(2) and O(1) models.  In the case of turned off
  fermions we do not really use those matrices in the base class
  DetModelGC and just set the matrix size to 1.
*/    
    static constexpr uint32_t MatrixSizeFactor = (OPDIM == 3 ? 4 : 2);

/*
    Originally the basic data type here used to be complex for O(3) or
    O(2) order parameters, real for O(1).  However, this is not compatible
    to having a magnetic field => need to use complex variables in all cases.
*/
    // typedef typename std::conditional<OPDIM==1, num, cpx>::type DataType;
    typedef cpx DataType;
    typedef arma::Mat<DataType> MatData;
    typedef typename arma::Mat<DataType>::template fixed<MatrixSizeFactor,MatrixSizeFactor> MatSmall;
    typedef arma::Col<DataType> VecData;
    typedef arma::Cube<DataType> CubeData;
    
    static num dataReal(const cpx& value) { return value.real(); } // could just use std::real
    static num dataReal(const num& value) { return value; }        //      -- " --
    static void setReal(cpx& value, num realPart) { value.real(realPart); }
    static void setReal(num& value, num realPart) { value = realPart; }
    static void setImag(cpx& value, num imagPart) { value.imag(imagPart); }
    static void setImag(num& value, num imagPart) { value = imagPart; }

    //need these to have the same code handle real/complex subviews
    //subviews do not have members set_real or set_imag.
    //For real subviews: just discard the imagPart.
    template<class Matrix1, class Matrix2>
    static void setRealImag(arma::subview<num> subv,
                            const Matrix1& realPart, const Matrix2& imagPart) {
        (void) imagPart;
        subv = realPart;
    }
    template<class Matrix1, class Matrix2>
    static void setRealImag(arma::subview<cpx> subv,
                            const Matrix1& realPart, const Matrix2& imagPart) {
        subv = MatCpx(realPart, imagPart);
    }
    
   
    
/*

    Things referenced from the base class

*/
    typedef DetModelGC<1, DataType> Base;
    // stupid C++ weirdness (visibility in different look-up phases)
    // forces us to explicitly "import" these protected base class
    // member variables:
    // (see: http://stackoverflow.com/questions/11405/gcc-problem-using-a-member-of-a-base-class-that-depends-on-a-template-argument )
    using Base::dtau;
    using Base::m;
    using Base::n;
    using Base::green;
    using Base::green_inv_sv;
    // using Base::greenFwd;
    // using Base::greenBwd;
    using Base::UdVStorage;
    using Base::eye_gc;
    using Base::lastSweepDir;
    using Base::obsScalar;
    using Base::obsVector;
    using Base::obsKeyValue;
    using Base::beta;
    using Base::s;
    using Base::currentTimeslice;
    using Base::loggingParams;
    // also for typedefs:
    typedef typename Base::UdVV UdVV;
    typedef typename Base::SweepDirection SweepDirection;
    // also for methods (g++ 4.7+)
    using Base::sweep_skeleton;
    using Base::sweepThermalization_skeleton;
    using Base::sweepSimple_skeleton;
    using Base::sweepSimpleThermalization_skeleton;
    using Base::setupUdVStorage_and_calculateGreen_skeleton;
    using Base::setupUdVStorage_and_calculateGreen_forTimeslice_skeleton;    
/*

    More type definitions

*/
    typedef VecNum::fixed<OPDIM> Phi;       //value of the (1,2,3)-component field at a single site and timeslice
    
    // small identity matrix (2x2 or 4x4)
    const MatSmall smalleye;
/*

    Random numbers

*/
    RngWrapper& rng;
    NormalDistribution normal_distribution;

/*

    Parameters and some necessary enums

*/
    // this contains the simulation parameters relevant to our model.
    // There is some duplication: the base class also has "free
    // standing" member variables dtau, m, beta, s
    ModelParams pars;
    // also in use: loggingParams from base class (see "using"
    // statement above)
    
    enum Band {XBAND = 0, YBAND = 1};
    enum Spin {SPINUP = 0, SPINDOWN = 1};
    enum BandSpin {XUP = 0, YDOWN = 1, XDOWN = 2, YUP = 3};
    static inline BandSpin getBandSpin(Band b, Spin s) {
        if      (b == XBAND and s == SPINUP)   return XUP;
        else if (b == YBAND and s == SPINDOWN) return YDOWN;
        else if (b == XBAND and s == SPINDOWN) return XDOWN;
        else /* (b == YBAND and s == SPINUP)*/ return YUP;
    }
    typedef ModelParamsDetSDW::BC_Type BC_Type;
    typedef ModelParamsDetSDW::UpdateMethod_Type UpdateMethod_Type;
    typedef ModelParamsDetSDW::SpinProposalMethod_Type SpinProposalMethod_Type;

/*

    Statistics about [global] updates
    
    This needs to be synchronized during parallel tempering

*/
    struct UpdateStatistics {
        uint32_t acceptedGlobalShifts;
        uint32_t attemptedGlobalShifts;
        uint32_t acceptedWolffClusterUpdates;
        uint32_t attemptedWolffClusterUpdates;
        uint32_t acceptedWolffClusterShiftUpdates;
        uint32_t attemptedWolffClusterShiftUpdates;
        num addedWolffClusterSize;

        UpdateStatistics() :
            acceptedGlobalShifts(0), attemptedGlobalShifts(0),
            acceptedWolffClusterUpdates(0), attemptedWolffClusterUpdates(0),
            acceptedWolffClusterShiftUpdates(0), attemptedWolffClusterShiftUpdates(0),
            addedWolffClusterSize(0.0)
        { }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int /*version*/) {
            ar & acceptedGlobalShifts;
            ar & attemptedGlobalShifts;
            ar & acceptedWolffClusterUpdates;
            ar & attemptedWolffClusterUpdates;
            ar & acceptedWolffClusterShiftUpdates;
            ar & attemptedWolffClusterShiftUpdates;
            ar & addedWolffClusterSize;
        }        
    } us;
    

/*
    
    Hopping constant representations

    These are useful for the checkerboard break-up (4-site matrices)
    in the case without a magnetic field.

*/  
    //hopping constants for XBAND and YBAND
    //these just contain the same values as t{x|y}{hor|ver} for historical reasons
    checkarray<num,2> hopHor;
    checkarray<num,2> hopVer;
    // sinh|cosh(-dtau * hop..)
    checkarray<num,2> sinhHopHor;
    checkarray<num,2> sinhHopVer;
    checkarray<num,2> coshHopHor;
    checkarray<num,2> coshHopVer;
    // sinh|cosh(-dtau * 0.5*hop..)
    checkarray<num,2> sinhHopHorHalf;
    checkarray<num,2> sinhHopVerHalf;
    checkarray<num,2> coshHopHorHalf;
    checkarray<num,2> coshHopVerHalf;


/*
    
  Precalculated 4-site hopping matrix exponentials

  These are used in the case of a perpendicular magnetic field.

  We consider the following plaquettes:
  
    subgroup == 0: plaquettes A = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
      i = (2m, 2n), with m,n integer, 2m < L, 2n < L
      j = i + XPLUS
      k = i + YPLUS
      l = k + XPLUS
    subgroup == 1: plaquettes B = [i j k l] = bonds (<ij>,<ik>,<kl>,<jl>)
      i = (2m+1, 2n+1), with m,n integer, 2m+1 < L, 2n+1 < L
      j = i + XPLUS
      k = i + YPLUS
      l = k + XPLUS

  To precalculate these, call precalc_4site_hopping_exponentials 
  
*/
    typedef MatData::fixed<4,4> Mat4Site;
    typedef checkarray<std::map<uint32_t, Mat4Site>, 2> ExpHop4SiteStorage; // indexed first by subgroup, then by site i

    // these are first indexed by band index XBAND or YBAND:
    checkarray<ExpHop4SiteStorage, 2> expHop4Site_minus;       // contains terms to compute exp(-dtau K_band)
    checkarray<ExpHop4SiteStorage, 2> expHop4Site_minusHalf;   // contains terms to compute exp(-0.5 dtau K_band)
    checkarray<ExpHop4SiteStorage, 2> expHop4Site_plus;        // contains terms to compute exp(+dtau K_band)
    checkarray<ExpHop4SiteStorage, 2> expHop4Site_plusHalf;    // contains terms to compute exp(+0.5 dtau K_band)

    
/*

    Per-orbital chemical potential

*/  
    checkarray<num, 2> mu;


/*

    Per-spinor-component magnetic field in z-direction (orthogonal to
    lattice plane), given in units of the flux quantum Phi_0

    use one field strength for (XUP, YDOWN) and its negative for
    (XDOWN, YUP)
  
*/
    checkarray<num, 4> zmag;    // to be indexed by enum BandSpin
    

/*    

    Neighbor tables
      
*/
    PeriodicSquareLatticeNearestNeighbors spaceNeigh;
    PeriodicChainNearestNeighbors<1> timeNeigh;



/*

    propK and variations

    These are used for non-checkerboard routines (non-performance-
    sensitive) when we need to involve the kinetic part of the
    Hamiltonian.

    They are complex to allow for a magnetic field.

*/
    checkarray<MatCpx, 2> propK; // e^{-dtau K[band]}
    MatCpx& propKx;
    MatCpx& propKy;

    //for shifting green functions to obtain equivalency of symmetric Trotter decomposition
    //[checker board decomposition could be applied alternatively]
    checkarray<MatCpx, 2> propK_half;       //factor of -dtau/2 in exponential
    MatCpx& propKx_half;
    MatCpx& propKy_half;
    checkarray<MatCpx, 2> propK_half_inv;   //factor of +dtau/2 in exponential
    MatCpx& propKx_half_inv;
    MatCpx& propKy_half_inv;



/*

    Green's function

    Note for OPDIM < 3: In the O(1) and O(2) model we only need to
    store and consider the upper left 2Nx2N block of the 4Nx4N
    matrices, corresponding to the XUP,YDOWN subspace.  The lower
    right block representing the XDOWN,YUP subspace is the complex
    conjugate of the upper block.  Transition probabilities are then
    given by the squared absolute value of the determinant of g.
*/    
    MatData& g;

    // Singular values of g^{-1}, valid after advance-steps /
    // setupudvstorage
    VecNum& g_inv_sv;


/*

    Auxiliary / bosonic fields

*/    
    // //three component sdw-order parameter,
    // //column indexes timeslice, row indexes site
    // MatNum phi0;
    // MatNum phi1;
    // MatNum phi2;
    
    // OPDIM-dimensional sdw-order parameter
    
    // slice indexes timeslice, column indexes order parameter
    // dimension, row indexes site    
    // [Index order: row, col, slice]
    // [Data layout; slice after slice, matrices then column major] 
    CubeNum phi;
    
    //discrete field for cdwU: l_i(\tau_k) == cdwl(i,k).
    //row indexes site, column indexes timeslice
    MatInt cdwl;
    //evaluation of element-wise functions of phi and cdwl:
    MatNum coshTermPhi;		// cosh(\lambda \dtau |\phi|)
    MatNum sinhTermPhi; 	// sinh(\lambda \dtau |\phi|) / |\phi|
    MatNum coshTermCDWl;	// cosh(\sqrt(\dtau) cdwU \eta)
    MatNum sinhTermCDWl;	// sinh(\sqrt(\dtau) cdwU \eta)


    
/*

    Adjustment of MC step size

    This needs to be synchronized during parallel tempering

*/   
    struct AdjustmentData {
        //adjustment of phiDelta, angleDelta, scaleDelta:
        // [these const's must also be defined in a source file, detsdwopdim.cpp;
        //  but the values are set here in the struct definition]
#if __GNUC__ == 4 && __GNUC_MINOR < 7
//work around g++ bug
#define const
#endif// __GNUC__ == 4 && __GNUC_MINOR < 7
        constexpr static const num InitialPhiDelta = 0.5;
        constexpr static const num InitialAngleDelta = 0.0;
        constexpr static const num InitialScaleDelta = 0.1;
        constexpr static const num MinScaleDelta = 0.0;
        constexpr static const num MaxScaleDelta = 1.0;
        constexpr static const num MinAngleDelta = -1.0;
        constexpr static const num MaxAngleDelta = 1.0;
        constexpr static const uint32_t AccRatioAdjustmentSamples = 100;
        constexpr static const num phiDeltaGrowFactor = 1.05;
        constexpr static const num phiDeltaShrinkFactor = 0.95;
#if __GNUC__ == 4 && __GNUC_MINOR < 7
#undef const
#endif// __GNUC__ == 4 && __GNUC_MINOR < 7

        num phiDelta;   //MC step size for field components (box update)
        num angleDelta; //  for rotation update (this is the minimal cos\theta)
        num scaleDelta; //  for scaling update (the standard deviation of the gaussian radius update)
        //used to adjust phiDelta/angleDelta/scaleDelta; acceptance ratios for local field updates
        num targetAccRatioLocal_phi;
        num lastAccRatioLocal_phi;
        RunningAverage accRatioLocal_box_RA, accRatioLocal_rotate_RA, accRatioLocal_scale_RA;

        //initialized to the respective constants above in the constructor
        num curminAngleDelta;
        num curmaxAngleDelta;
        num curminScaleDelta;
        num curmaxScaleDelta;

        AdjustmentData(const ModelParams& pars) :
            phiDelta(InitialPhiDelta), angleDelta(InitialAngleDelta),
            scaleDelta(InitialScaleDelta),
            targetAccRatioLocal_phi(pars.accRatio),
            lastAccRatioLocal_phi(0),
            accRatioLocal_box_RA(AccRatioAdjustmentSamples),
            accRatioLocal_rotate_RA(AccRatioAdjustmentSamples),
            accRatioLocal_scale_RA(AccRatioAdjustmentSamples),
            curminAngleDelta(MinAngleDelta), curmaxAngleDelta(MaxAngleDelta),
            curminScaleDelta(MinScaleDelta), curmaxScaleDelta(MaxScaleDelta)
        { }

        template<class Archive>
        void serialize(Archive & ar, const unsigned int /*version*/) {
            ar & phiDelta;   
            ar & angleDelta; 
            ar & scaleDelta; 
            ar & targetAccRatioLocal_phi;
            ar & lastAccRatioLocal_phi;
            ar & accRatioLocal_box_RA
               & accRatioLocal_rotate_RA
               & accRatioLocal_scale_RA;
            ar & curminAngleDelta;
            ar & curmaxAngleDelta;
            ar & curminScaleDelta;
            ar & curmaxScaleDelta;
        }

        
        
        // static constexpr uint32_t bufsize = numelemens * sizeof(double);

        // void get_data_from_buffer(const uint8_t* buffer) {
        //     const double* double_buffer = reinterpret_cast<const double*>(buffer);
            
        //     phiDelta = *buffer;
        //     angleDelta = *(buffer + 1);
        //     scaleDelta = *(buffer + 2);
        //     curminAngleDelta = *(buffer + 3);
        //     curmaxAngleDelta = *(buffer + 4);
        //     curminScaleDelta = *(buffer + 5);
        //     curmaxScaleDelta = *(buffer + 6);
        //     assert(bufsize - 1 == 6);
        //     // reset acc ratio running averagers
        //     accRatioLocal_box_RA = RunningAverage(AccRatioAdjustmentSamples);
        //     accRatioLocal_rotate_RA = RunningAverage(AccRatioAdjustmentSamples);
        //     accRatioLocal_scale_RA = RunningAverage(AccRatioAdjustmentSamples);
            
        // }

        // void set_buffer_to_data(double* buffer) const {
        //     *buffer = phiDelta;
        //     *(buffer + 1) = angleDelta;
        //     *(buffer + 2) = scaleDelta;
        //     *(buffer + 3) = curminAngleDelta;
        //     *(buffer + 4) = curmaxAngleDelta;
        //     *(buffer + 5) = curminScaleDelta;
        //     *(buffer + 6) = curmaxScaleDelta;
        //     assert(bufsize - 1 == 6);
        // }
    } ad;


/*
    
    internal counter of performed sweeps. This should be serialized

*/
    uint32_t performedSweeps;


/*
    
    Observables:

*/
    Phi meanPhi;		//averaged field [~ magnetization]
    num normMeanPhi;            //norm of averaged field

    // for O(2) model without fermions: parts of the bosonic spin stiffness
    // rhoS = (beta / L**2) * ( <Gc> + <Gs>**2 - <Gs**2> )
    num phiRhoS_Gs;
    num phiRhoS_Gc;

    num associatedEnergy;       // to enable reweighting measure: 1/2 1/(N*m) sum phi^2

    // // diagonal blocks of the (equal time) Green's function in
    // // momentum space
    // VecNum kgreenXUP;
    // VecNum kgreenYDOWN;
    // // the following two are equal to the upper to in the cases of
    // // OPDIM == 1 or 2
    // VecNum kgreenXDOWN;
    // VecNum kgreenYUP;
    // helpers for the above: real space, summed and averaged over
    // imaginary time during a sweep
    MatData greenXUPXUP_summed;
    MatData greenYDOWNYDOWN_summed;
    MatData greenXDOWNXDOWN_summed;
    MatData greenYUPYUP_summed;
    // also include some of the off-diagonal blocks:
    MatData greenXUPYDOWN_summed;
    MatData greenYDOWNXUP_summed;
    
    // greenK0 is the sum over all entries of G
    num greenK0;
    // greenLocal is the trace of G divided by the size of G
    num greenLocal;
        
    
    checkarray<VecNum, 2> kOcc;     //Fermion occupation number in momentum space for x/y-band; site-index: k-vectors
    VecNum& kOccX;
    VecNum& kOccY;
//    checkarray<VecNum, 2> kOccImag;
//    VecNum& kOccXimag;
//    VecNum& kOccYimag;

    // checkarray<VecNum, 2> occ;     //Fermion occupation number in Real space for x/y-band; indexed by site
    // VecNum& occX;
    // VecNum& occY;
//    checkarray<VecNum, 2> occImag;
//    VecNum& occXimag;
//    VecNum& occYimag;

    //pairing correlations near maximum range x_max = (L/2, L/2)
    num pairPlusMax;
    num pairMinusMax;
//    num pairPlusMaximag;
//    num pairMinusMaximag;
    //and between site 0 and any other site
    VecNum pairPlus;
    VecNum pairMinus;
//    VecNum pairPlusimag;
//    VecNum pairMinusimag;

    //fermion energy
    num fermionEkinetic;        //kinetic
//    num fermionEkinetic_imag;
    num fermionEcouple;         //coupling
//    num fermionEcouple_imag;

    // band occupation / charge correlations
    //
    // Real space for ab = xx, xy, yx, yy band combinations
    // rows: site i, cols: site j
    // <n_{a,i} n_{b,j}>
    // Here: Do not use symmat because of site resolution
    // GenMat<MatNum, 2> occCorr;  // index by [XBAND,YBAND] etc.
    // // charge-charge correlations
    // // <n_i n_j> = <n_{x,i}n_{x,i}> + <n_{x,i}n_{y,i}> + <n_{y,i}n_{x,i}> + <n_{y,i}n_{y,i}>
    // MatNum chargeCorr;
    // Fourier transforms, index by wave-vector k = (kx,ky)
    // represented as ksite (like kOcc)
    //
    // We only need to write out the measurements in the form of the
    // fourier-transformed version
    // GenMat<VecNum, 2> occCorrFT;
    // VecNum chargeCorrFT;

    // This measures the local occupation difference of the x and y
    // orbitals: 1/N \sum_i <(n_{xi} - n_{yi})^2>
    num occDiffSq;


/*

    Static helper member functions, defined in header

*/
    static inline std::string bandstr(Band b);
    static inline std::string spinstr(Spin s);
    static inline std::string bandspinstr(BandSpin bs);
    static inline std::string updateMethodstr(UpdateMethod_Type um);
    static inline std::string spinProposalMethodstr(SpinProposalMethod_Type sp);


    //lookup values for discrete field:
    static inline num cdwl_gamma(int32_t l);
    static inline num cdwl_eta(int32_t l);


/*
  
    Helpers for functors, defined in header
    
*/
    //these for_each functions don't really work well with the class template
    template<typename Callable>
    void for_each_band(Callable func);
    template<typename Callable> inline
    void for_each_site(Callable func);
    template<typename Callable> inline
    void for_each_timeslice(Callable func);
    template<typename CallableSiteTimeslice, typename V> inline
    V sumWholeSystem(CallableSiteTimeslice f, V init);
    template<typename CallableSiteTimeslice, typename V> inline
    V averageWholeSystem(CallableSiteTimeslice f, V init);

/*
  
    Helpers, defined in header
    
*/
    uint32_t coordsToSite(uint32_t x, uint32_t y) const {
        return y*pars.L + x;
    }
    Phi getPhi(uint32_t site, uint32_t timeslice) const {
//      assert(site >= 0);    //fulfilled always for uint
        assert(site <  pars.N);
        assert(timeslice >= 1);
        assert(timeslice <= pars.m);
        Phi result;
        if (OPDIM == 3 or OPDIM == 2 or OPDIM == 1) {
            result[0] = phi(site, 0, timeslice);
        }
        if (OPDIM == 3 or OPDIM == 2) {
            result[1] = phi(site, 1, timeslice);
        }
        if (OPDIM == 3) {
            result[2] = phi(site, 2, timeslice);
        }
        return result;
    }
    void setPhi(uint32_t site, uint32_t timeslice, Phi new_phi) {
        for (uint32_t dim = 0; dim < OPDIM; ++dim) {
            phi(site, dim, timeslice) = new_phi[dim];
        }
    }


/*
  
    Functions related to auxiliary / bosonic fields, setting up propK,
    setting up the UdV storage and setting up the Green's function
    
*/
    void setupRandomField();
    void setupConstantField();    
    // return (coshTerm*, sinhTerm*) for these field values
    std::tuple<num,num> getCoshSinhTermPhi(Phi phi);
    std::tuple<num,num> getCoshSinhTermCDWl(int32_t cdwl);

    //compute new entries of coshTerm* and sinhTerm* from current phi, cdwl
    void updateCoshSinhTermsPhi();
    void updateCoshSinhTermsCDWl();
    void updateCoshSinhTerms();
    //the same, for a single site
    void updateCoshSinhTermsPhi(uint32_t site, uint32_t timeslice);
    void updateCoshSinhTermsCDWl(uint32_t site, uint32_t timeslice);
    void updateCoshSinhTerms(uint32_t site, uint32_t timeslice);
    void setupPropK();          //compute e^(-dtau*K..) matrices by diagonalization
    void setupUdVStorage_and_calculateGreen();
    void setupUdVStorage_and_calculateGreen_forTimeslice(uint32_t timeslice);

    
/*
  
    Checkerboard multiplication routines
    
*/
    //the following are template functions to allow applying them
    //to submatrices as well
    //"cb" = "checkerboard"
    //Reminder: These do not include chemical potential
    //if invertedCbOrder = false: use the following checkerbaord decomposition:
    //    e^{+-K} = e^{+-K_b} e^{+-K_a}
    //if invertedCbOrder = true: use the following checkerbaord decomposition:
    //    e^{+-K} = e^{+-K_a} e^{+-K_b}
    //for the symmetric checkerboard break-up (CB_ASSAAD_BERG) this is ignored
    // with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to E^(sign * dtau * K_band) * A
    template <class Matrix>
    MatData cbLMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    // with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to A * E^(sign * dtau * K_band)
    template <class Matrix>
    MatData cbRMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder = false);

    //cbLMultHoppingExp and cbRMultHoppingExp need separate implementations for each CheckerboardMethod,
    //this cannot be realized by a direct partial template specialization, but we need to have a proxy
    //with function overloading
    //compare first solution in winning answer at:
    //http://stackoverflow.com/questions/1501357/template-specialization-of-particular-members
    template <class Matrix>
    MatData cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
                                   const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatData cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
                                   const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatData cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
                                   const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatData cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
                                   const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    //functions called by the above if no magnetic field is applied:
    template<class Matrix>
    void cb_assaad_applyBondFactorsLeft(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver);
    template<class Matrix>
    void cb_assaad_applyBondFactorsRight(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver);

    // in the case with a non-zero magnetic field in z-direction, we
    // need to compute the matrix exponentials of the 4-site (plaquette)
    // hopping matrices numerically in advance:
    void precalc_4site_hopping_exponentials();
    template<class Matrix>    
    void cb_assaad_applyBondFactorsLeft_precalcedMatrices(Matrix& result, uint32_t subgroup,
                                                          const ExpHop4SiteStorage& expHop4SiteMatrices);
    template<class Matrix>    
    void cb_assaad_applyBondFactorsRight_precalcedMatrices(Matrix& result, uint32_t subgroup,
                                                           const ExpHop4SiteStorage& expHop4SiteMatrices);
    
    

    
    // the following function prepares numerically exponentiated 4x4 

    //the following take a MatrixSizeFactor*N x MatrixSizeFactor*N
    //matrix A and effectively multiply B(k2,k1) or its inverse to the
    //left or right of it and return the result
    MatData checkerboardLeftMultiplyBmat(const MatData& A, uint32_t k2, uint32_t k1);
    MatData checkerboardRightMultiplyBmat(const MatData& A, uint32_t k2, uint32_t k1);
    MatData checkerboardLeftMultiplyBmatInv(const MatData& A, uint32_t k2, uint32_t k1);
    MatData checkerboardRightMultiplyBmatInv(const MatData& A, uint32_t k2, uint32_t k1);

    //helpers for the checkerboardMultiplyFunctions
    MatData rightMultiplyBk(const MatData& orig, uint32_t k);     //multiply B(k,k-1) from right to orig, return result
    MatData rightMultiplyBkInv(const MatData& orig, uint32_t k);  //multiply B(k,k-1)^-1 from right to orig, return result
    MatData leftMultiplyBk(const MatData& orig, uint32_t k);      //multiply B(k,k-1) from left to orig, return result
    MatData leftMultiplyBkInv(const MatData& orig, uint32_t k);   //multiply B(k,k-1)^-1 from left to orig, return result

/*
  
    Checkerboard or dense computation of B-matrix, some helpers for both cases
    
*/
    //compute B-matrix using dense matrix products or checkerboardLeftMultiplyBmat
    MatData computeBmatSDW(uint32_t k2, uint32_t k1);            

    // compute sqrt(dtau) * cdwU * eta_{cdwl[i]}
    template<class Vec>
    VecNum compute_d_for_cdwl(const Vec& cdwl);
    // compute sqrt(dtau) * cdwU * eta_{cdwl} for cdwl
    num compute_d_for_cdwl_site(int32_t cdwl);


    
/*

    Routines for consistency checks / debugging
    
*/
    // Compute e^( sign * dtau * V(phi-configuration) ) as a dense matrix.
    // for one timeslice. This may be useful for checks.
    MatData computePotentialExponential(int sign, checkarray<VecNum, OPDIM> phi, VecInt cdwl);

    //several checks, many are normally commented out
    virtual void consistencyCheck();

    // will be used in DetModelGC::advance{Up,Down} functions

    // Compares the Green's function as computed from scratch with the
    // wrapped version.
    void greenConsistencyCheck(const MatData& g1, const MatData& g2, SweepDirection cur_sweep_dir);


    // reference computation of the Green's function determinant
    // ratio.  Compare this during the estimation of the fermionic
    // transition probability for going from the current to a new
    // phi-spin configuration
    num computeGreenDetRatioFromScratch(uint32_t timeslice, const CubeNum& newPhi);
    // helper wrapping the above for a single spin update
    num computeGreenDetRatioFromScratch(uint32_t site, uint32_t timeslice, Phi singleNewPhi);


    // reference computation of the new Green's function after
    // switching to the new phi-spin configuration
    MatData computeGreenFromScratch(uint32_t timeslice, const CubeNum& newPhi);
    // helper wrapping the above for a single spin update
    MatData computeGreenFromScratch(uint32_t site, uint32_t timeslice, Phi singleNewPhi);    
/*
  
    Routines related to measurements
    
*/
    //Perform the correction of the Green's function to ensure an
    //effectively symmetric Trotter decomposition.  This returns the
    //shifted matrix for the current timeslice.  This should be done
    //before measurements.
    virtual MatData shiftGreenSymmetric();
    
    //the upper function shiftGreenSymmetric() calls the following
    //helper with functors RightMultiply, LeftMultiply depending on
    //the CheckerboardMethod
    template<class RightMultiply, class LeftMultiply>
    MatData shiftGreenSymmetric_impl(RightMultiply, LeftMultiply);

    //measuring observables
    void initMeasurements();				//reset stored observable values (beginning of a sweep)
    void measure(uint32_t timeslice);                   //measure observables for one timeslice
    void finishMeasurements();				//finalize stored observable values (end of a sweep)
    std::set<uint32_t> timeslices_included_in_measurement; 	//for a consistency check -- sweep includes correct #timeslices
    // compute the structure factor from a matrix of real space correlations
    void computeStructureFactor(VecNum& out_k, const MatNum& in_r);
    void computeStructureFactor(VecNum& out_k, const MatCpx& in_r); // this computes the real part of the Fourier transform of in_r

    
/*

    Monte Carlo updates, some related data structures
    
*/
    //update the auxiliary field and the green function in the single timeslice,
    //using Metropolis updates
    void updateInSlice(uint32_t timeslice);
    //this template member function calls the right one of those below
    //return acceptance ratio for that sweep
    template<class CallableProposeSpin>
    num callUpdateInSlice_for_updateMethod(uint32_t timeslice, CallableProposeSpin proposeSpin) {
    	switch(pars.updateMethod) {
    	case UpdateMethod_Type::ITERATIVE:
    		return updateInSlice_iterative(timeslice, proposeSpin);
    	case UpdateMethod_Type::WOODBURY:
    		return updateInSlice_woodbury(timeslice, proposeSpin);
    	case UpdateMethod_Type::DELAYED:
    		return updateInSlice_delayed(timeslice, proposeSpin);
    	default:
    		return 0;
    	}
    };
    //specific implementations of updateInSlice, the callable parameter is for one of the field proposal routines below,
    //return the acceptance ratio for that sweep
    //CallableProposeLocalUpdate should have the signature tuple<Changed,Phi,uint32_t>(uint32_t site, uint32_t timeslice).
    //This either changes phi or cdwl locally and returns the new values for both.
    //The returned enum Changed is defined a bit below and contains what has changed or whether the update
    //is rejected (e.g. negative magnitude of phi)
    template<class CallableProposeLocalUpdate>
    num updateInSlice_iterative(uint32_t timeslice, CallableProposeLocalUpdate proposeLocalUpdate);
    template<class CallableProposeLocalUpdate>
    num updateInSlice_woodbury(uint32_t timeslice, CallableProposeLocalUpdate proposeLocalUpdate);
   	//this one uses a nested struct because for some status
    template<class CallableProposeLocalUpdate>
    num updateInSlice_delayed(uint32_t timeslice, CallableProposeLocalUpdate proposeLocalUpdate);
    struct DelayedUpdatesData {		//some helper data for updatesInSlice_delayed that should not be realloced all the time
    	MatData X;
    	MatData Y;
    	MatData Rj;
        MatSmall Sj;
    	MatData Cj;
    	MatSmall tempBlock;
	MatSmall Mj;
    	DelayedUpdatesData(uint32_t N, uint32_t delaySteps)
            : X(MatrixSizeFactor*N, MatrixSizeFactor*delaySteps),
              Y(MatrixSizeFactor*delaySteps, MatrixSizeFactor*N),
              Rj(MatrixSizeFactor, MatrixSizeFactor*N),
              Cj(MatrixSizeFactor*N, MatrixSizeFactor) {
            X.zeros();
            Y.zeros();
    	}
    } dud;

    //this one does some adjusting of the box size from which new fields are chosen:
    void updateInSliceThermalization(uint32_t timeslice);

    //functions used by updateInSlice:
    //--------------------------------
    enum Changed {
    	NONE = 0, 	// reject update
    	PHI = 1,	// new phi, consider bosonic action & fermion determinant
    	CDWL = 2	// new discrete field, consider only fermion determinant
    };
    typedef std::tuple<Changed,Phi,int32_t> changedPhiInt;
    //generate a new phi:
    changedPhiInt proposeNewPhiBox(uint32_t site, uint32_t timeslice);		//from a box
    changedPhiInt proposeRotatedPhi(uint32_t site, uint32_t timeslice);		//same length, new angles
    changedPhiInt proposeScaledPhi(uint32_t site, uint32_t timeslice);		//new length, same angles
    changedPhiInt proposeRotatedScaledPhi(uint32_t site, uint32_t timeslice);	//new length, new angles
    //generate a new cdwl:
    changedPhiInt proposeNewCDWl(uint32_t site, uint32_t timeslice);		//choose Metropolis-randomly from +-1, +-2
    //more generic helpers
    num deltaSPhi(uint32_t site, uint32_t timeslice, Phi newphi);
    MatSmall get_delta_forsite(
        Phi newphi, int32_t new_cdwl, uint32_t timeslice, uint32_t site);

    // over-relaxation update
    //-----------------------
    // this is only for simulations without fermions.  This variation
    // picks out random sites in a single timeslice.  This makes less
    // sense than just going linearly through the entire space-time
    // lattice.
    void updateInSlice_overRelaxation(uint32_t timeslice);


    // Entire space-time sweeps of local updates
    // -----------------------------------------
    // This is only for simulations without fermions.  Linear sweep over
    // entire space-time lattice.
    void overRelaxationSweep();

    //global updates
    //--------------
    // These are not restricted to operation on a single timeslice.
    // True global updates: cluster / global shift moves.  For
    // non-fermion simulations: overRelaxationSweep is also to be
    // called in globalMove()

    void globalMove();
    //Try a global move, where the fields on all sites and timeslices are shifted
    //by the same constant amount.  Only do this after a certain number of sweeps.
    void attemptGlobalShiftMove();
    //Try a Wolff single-cluster update. The cluster is constructed using
    //probabilities determined from the scalar action [proposal probability].
    //Then it is accepted by a Metropolis rule according to the fermion
    //determinants.
    void attemptWolffClusterUpdate();
    // attempt a combined update
    void attemptWolffClusterShiftUpdate(); 
    struct GlobalMoveData {		//some helper data that should not be reallocated all the time
        // quantities to be backed up:
        CubeNum phi;
    	MatNum coshTermPhi;
    	MatNum sinhTermPhi;
    	MatData g;
        VecNum  g_inv_sv;
    	std::unique_ptr<checkarray<std::vector<UdVV>, 1>> UdVStorage;

    	//for the cluster update: mark visited sites by 1, else 0
    	MatUint visited;		//as usual: row is spatial index, col is time index
    	typedef std::tuple<uint32_t, uint32_t> SpaceTimeIndex;
    	//for the cluster update: keep book about which sites to try to add next
    	std::stack<SpaceTimeIndex> next_sites;
    	GlobalMoveData(uint32_t N, uint32_t m, bool noFermions = false) {
            phi.resize(N, OPDIM, m+1);
            visited.resize(N, m+1);
            next_sites = std::stack<SpaceTimeIndex>();            

            if (not noFermions) {
                coshTermPhi.resize(N, m+1);
                sinhTermPhi.resize(N, m+1);
                g.resize(MatrixSizeFactor*N, MatrixSizeFactor*N);
                g_inv_sv.resize(MatrixSizeFactor*N);
                UdVStorage = std::unique_ptr<checkarray<std::vector<UdVV>, 1>>(
                    new checkarray<std::vector<UdVV>, 1>);
            }                
        }
    } gmd;
    //helper functions for global updates:
    void addGlobalRandomDisplacement(); // works directly on phi0,phi1,phi2
    uint32_t buildAndFlipCluster(bool updateCoshSinh = true); // returns size of cluster
    void globalMoveStoreBackups();
    void globalMoveRestoreBackups();
    

    //compute the total value of the action associated with the field phi
    num phiAction();



    
/*    

      wrappers to use to instantiate template functions of the base class

*/
    struct sdwRightMultiplyBmatInv;
    struct sdwLeftMultiplyBmatInv;
    struct sdwRightMultiplyBmat;
    struct sdwLeftMultiplyBmat;
    struct sdwComputeBmat;


/*

     For debugging / consistency purposes:
  
 */
    // wrapping/advancing green consistency checks (optional)
    struct GreenConsistencyLogger {
        std::string logfiledir;
        // greenConsistency results during sweep up / down:
        std::ofstream up_log, down_log;
        
        GreenConsistencyLogger(const std::string& logfiledir_, bool enabled);
    } greenConsistencyLogger;

    // used for further logging / consistency checks
    std::unique_ptr<DoubleVectorWriterSuccessive> detRatioLogging, greenLogging;
public:
/*

    Serialization

*/
    // serialization by DetQMC methods
    template<class Archive>
    void saveContents(Archive &ar) {
        Base::saveContents(ar);            //base class
        serializeContentsCommon(ar);
    }

    //after loadContents() a sweep must be performed before any measurements are taken:
    //else the green function would not be in a valid state
    template<class Archive>
    void loadContents(Archive &ar) {
        Base::loadContents(ar);            //base class
        serializeContentsCommon(ar);
        //the fields now have a valid state, update UdV-storage to start
        //sweeping again
        setupUdVStorage_and_calculateGreen();
        //now: lastSweepDir == SweepDirection::Up --> the next sweep will be downwards
    }

    template<class Archive>
    void serializeContentsCommon(Archive& ar) {
    	// ar & acceptedGlobalShifts;
    	// ar & attemptedGlobalShifts;
    	// ar & acceptedWolffClusterUpdates;
    	// ar & attemptedWolffClusterUpdates; 
        // ar & acceptedWolffClusterShiftUpdates;  
        // ar & attemptedWolffClusterShiftUpdates;         
    	// ar & addedWolffClusterSize;
        ar & us;
        ar & phi;
        ar & cdwl;
        ar & coshTermPhi & sinhTermPhi;
        ar & coshTermCDWl & sinhTermCDWl;
        // ar & ad.phiDelta & ad.angleDelta & ad.scaleDelta;
        // ar & ad.targetAccRatioLocal_phi & ad.lastAccRatioLocal_phi;
        // ar & ad.accRatioLocal_box_RA & ad.accRatioLocal_rotate_RA & ad.accRatioLocal_scale_RA;
        // ar & ad.curminAngleDelta & ad.curmaxAngleDelta;
        // ar & ad.curminScaleDelta & ad.curmaxScaleDelta;
        ar & ad;
        ar & performedSweeps;
    }
};





/*

	Static helper member functions

*/
template <CheckerboardMethod CBM, int OPDIM>
inline std::string DetSDW<CBM, OPDIM>::bandstr(Band b) {
	return (b == XBAND) ? "x" : (b == YBAND ? "y" : "N");
}
template <CheckerboardMethod CBM, int OPDIM>
inline std::string DetSDW<CBM, OPDIM>::spinstr(Spin s) {
	return (s == SPINUP) ? "up" : (s == SPINDOWN ? "dn" : "N");
}
template <CheckerboardMethod CBM, int OPDIM>
inline std::string DetSDW<CBM, OPDIM>::bandspinstr(BandSpin bs) {
    switch (bs) {
    case XUP:   return "XUp";
    case XDOWN: return "XDown";
    case YUP:   return "YUp";
    case YDOWN: return "YDown";
    default:    return "N";
    }
}
template <CheckerboardMethod CBM, int OPDIM>
inline std::string DetSDW<CBM, OPDIM>::updateMethodstr(UpdateMethod_Type um) {
	switch (um) {
	case UpdateMethod_Type::ITERATIVE:
		return "iterative";
	case UpdateMethod_Type::WOODBURY:
		return "woodbury";
	case UpdateMethod_Type::DELAYED:
		return "delayed";
	default:
		return "invalid";
	}
}
template <CheckerboardMethod CBM, int OPDIM>
inline std::string DetSDW<CBM, OPDIM>::spinProposalMethodstr(SpinProposalMethod_Type sp) {
	switch (sp) {
	case SpinProposalMethod_Type::BOX:
		return "box";
	case SpinProposalMethod_Type::ROTATE_THEN_SCALE:
		return "rotate_then_scale";
	case SpinProposalMethod_Type::ROTATE_AND_SCALE:
		return "rotate_and_scale";
	default:
		return "invalid";
	}
}

    

//lookup values for discrete field:
template <CheckerboardMethod CBM, int OPDIM>
inline num DetSDW<CBM, OPDIM>::cdwl_gamma(int32_t l) {
    switch (l) {
    case +1:
    case -1:
        return (3. + std::sqrt(6.));
    case +2:
    case -2:
        return (3. - std::sqrt(6.));
    default:
        return 0;
    }
}
template <CheckerboardMethod CBM, int OPDIM>
inline num DetSDW<CBM, OPDIM>::cdwl_eta(int32_t l) {
    switch (l) {
    case +1:
        return  std::sqrt(2. * (3. - std::sqrt(6.)));
    case -1:
        return -std::sqrt(2. * (3. - std::sqrt(6.)));
    case +2:
        return  std::sqrt(2. * (3. + std::sqrt(6.)));
    case -2:
        return -std::sqrt(2. * (3. + std::sqrt(6.)));
    default:
        return 0;
    }
}


/*

    Helpers for functors

*/
template<CheckerboardMethod CBM, int OPDIM>
template<typename Callable>
void DetSDW<CBM, OPDIM>::for_each_band(Callable func) {
    func(XBAND);
    func(YBAND);
}
template<CheckerboardMethod CBM, int OPDIM>
template<typename Callable> inline
void DetSDW<CBM, OPDIM>::for_each_site(Callable func) {
    for (uint32_t site = 0; site < pars.N; ++site) {
        func(site);
    }
}
template<CheckerboardMethod CBM, int OPDIM>
template<typename Callable> inline
void DetSDW<CBM, OPDIM>::for_each_timeslice(Callable func) {
    for (uint32_t k = 1; k <= m; ++k) {
        func(k);
    }
}
// template<typename CallableSiteTimeslice, typename V> inline
// V sumWholeSystem(CallableSiteTimeslice f, V init);

template<CheckerboardMethod CBM, int OPDIM>
template<typename CallableSiteTimeslice, typename V> inline
V DetSDW<CBM, OPDIM>::sumWholeSystem(CallableSiteTimeslice f, V init) {
    for (uint32_t timeslice = 1; timeslice <= pars.m; ++timeslice) {
        for (uint32_t site = 0; site < pars.N; ++site) {
            init += f(site, timeslice);
        }
    }
    return init;
}
template<CheckerboardMethod CBM, int OPDIM>
template<typename CallableSiteTimeslice, typename V> inline
V DetSDW<CBM, OPDIM>::averageWholeSystem(CallableSiteTimeslice f, V init) {
    V sum = sumWholeSystem(f, init);
    return sum / num(pars.m * pars.N);
}




/*    

      wrappers to use to instantiate template functions of the base class

*/
template<CheckerboardMethod CBM, int OPDIM>
struct DetSDW<CBM, OPDIM>::sdwComputeBmat {
    DetSDW<CBM, OPDIM>* parent;
    sdwComputeBmat(DetSDW<CBM, OPDIM>* parent_) :
        parent(parent_)
        { }
    MatData operator()(uint32_t gc, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        return parent->computeBmatSDW(k2, k1);
    }
};

template<CheckerboardMethod CBM, int OPDIM>
struct DetSDW<CBM, OPDIM>::sdwLeftMultiplyBmat {
    DetSDW<CBM, OPDIM>* parent;
    sdwLeftMultiplyBmat(DetSDW<CBM, OPDIM>* parent_) :
        parent(parent_)
    { }
    MatData operator()(uint32_t gc, const MatData& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (CBM != CB_NONE) {
            return parent->checkerboardLeftMultiplyBmat(mat, k2, k1);
        } else {
            return parent->computeBmatSDW(k2, k1) * mat;
        }
    }
};

template<CheckerboardMethod CBM, int OPDIM>
struct DetSDW<CBM, OPDIM>::sdwRightMultiplyBmat {
    DetSDW<CBM, OPDIM>* parent;
    sdwRightMultiplyBmat(DetSDW<CBM, OPDIM>* parent_) :
        parent(parent_)
        { }
    MatData operator()(uint32_t gc, const MatData& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (CBM != CB_NONE) {
            return parent->checkerboardRightMultiplyBmat(mat, k2, k1);
        } else {
            return mat * parent->computeBmatSDW(k2, k1);
        }
    }
};

template<CheckerboardMethod CBM, int OPDIM>
struct DetSDW<CBM, OPDIM>::sdwLeftMultiplyBmatInv {
    DetSDW<CBM, OPDIM>* parent;
    sdwLeftMultiplyBmatInv(DetSDW<CBM, OPDIM>* parent_) :
        parent(parent_)
        { }
    MatData operator()(uint32_t gc, const MatData& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (CBM != CB_NONE) {
            return parent->checkerboardLeftMultiplyBmatInv(mat, k2, k1);
        } else {
            return arma::inv(parent->computeBmatSDW(k2, k1)) * mat;
        }
    }
};

template<CheckerboardMethod CBM, int OPDIM>
struct DetSDW<CBM, OPDIM>::sdwRightMultiplyBmatInv {
    DetSDW<CBM, OPDIM>* parent;
    sdwRightMultiplyBmatInv(DetSDW<CBM, OPDIM>* parent_) :
        parent(parent_)
        { }
    MatData operator()(uint32_t gc, const MatData& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (CBM != CB_NONE) {
            return parent->checkerboardRightMultiplyBmatInv(mat, k2, k1);
        } else {
            return mat * arma::inv(parent->computeBmatSDW(k2, k1));
        }
    }
};




// return min{1, e^{-\Delta}}
template<>
num get_replica_exchange_probability<DetSDW<CB_NONE, 1>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_ASSAAD_BERG, 1>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_NONE, 2>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_ASSAAD_BERG, 2>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_NONE, 3>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_ASSAAD_BERG, 3>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
 



#endif /* DETSDWOPDIM_H_ */
