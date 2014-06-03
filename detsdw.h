/*
 * detsdw.h
 *
 *  Created on: Feb 21, 2013
 *      Author: gerlach
 */

#ifndef DETSDW_H_
#define DETSDW_H_

#include <tuple>
#include <vector>
#include <complex>
#include "rngwrapper.h"
#include "detmodel.h"
#include "detsdwparams.h"
#include "neighbortable.h"
#include "RunningAverage.h"
#include "checkarray.h"
#include "normaldistribution.h"
#include "symmat.h"

typedef std::complex<num> cpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<cpx> SpMatCpx;
typedef arma::Col<cpx> VecCpx;
typedef arma::Cube<cpx> CubeCpx;

class SerializeContentsKey;

std::string cbmToString(CheckerboardMethod cbm);


template<CheckerboardMethod CBM> class DetSDW;

template<>
std::unique_ptr<DetSDW<CB_NONE>>
createReplica<DetSDW<CB_NONE>,ModelParamsDetSDW>(RngWrapper& rng, ModelParamsDetSDW pars);
template<>
std::unique_ptr<DetSDW<CB_ASSAAD_BERG>>
createReplica<DetSDW<CB_ASSAAD_BERG>,ModelParamsDetSDW>(RngWrapper& rng, ModelParamsDetSDW pars);


// return min{1, e^{-\Delta}}
template<>
num get_replica_exchange_probability<DetSDW<CB_NONE>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);
template<>
num get_replica_exchange_probability<DetSDW<CB_ASSAAD_BERG>>(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2);


// template parameters: do a checker-board decomposition of
// which kind?
template<CheckerboardMethod Checkerboard>
class DetSDW: public DetModelGC<1, cpx> {
public:
    typedef ModelParamsDetSDW ModelParams;
private:
    DetSDW(RngWrapper& rng, const ModelParams& pars);
public:
    friend std::unique_ptr<DetSDW<Checkerboard>>
    createReplica<DetSDW<Checkerboard>,ModelParamsDetSDW>(RngWrapper& rng, ModelParams pars);

    virtual ~DetSDW();
    virtual uint32_t getSystemN() const;
    virtual MetadataMap prepareModelMetadataMap() const;

    virtual void thermalizationOver();

    virtual void sweep(bool takeMeasurements);
    virtual void sweepThermalization();
    virtual void sweepSimple(bool takeMeasurements);
    virtual void sweepSimpleThermalization();

    //Methods to implement a replica-exchange / parallel tempering scheme
    //-------------------------------------------------------------------
    // get / set pars.r
    num  get_exchange_parameter_value() const;
    void set_exchange_parameter_value(num r);
    constexpr std::string get_exchange_parameter_name() const;
    // return 1/2 \int_0^\beta d\tau \sum_i [\vec{\phi}_i(\tau)]^2
    num get_exchange_action_contribution() const;

    //Perform the correction of the Green's function to ensure an
    //effectively symmetric Trotter decomposition.  This returns the
    //shifted matrix for the current timeslice.  This should be done
    //before measurements.
    virtual MatCpx shiftGreenSymmetric();
protected:
/*

    Things referenced from the base class

*/
    typedef DetModelGC<1, cpx> Base;
    // stupid C++ weirdness forces us to explicitly "import" these protected base
    // class member variables:
    // (see: http://stackoverflow.com/questions/11405/gcc-problem-using-a-member-of-a-base-class-that-depends-on-a-template-argument )
    using Base::dtau;
    using Base::m;
    using Base::n;
    using Base::green;
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
    typedef typename Base::UdVV UdVV;

/*

    More type definitions

*/
    typedef VecNum::fixed<3> Phi;       //value of the three-component field at a single site and timeslice
    
    //complex 4x4 identity matrix
    const MatCpx::fixed<4,4> eye4cpx;

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
    
    enum Band {XBAND = 0, YBAND = 1};
    enum Spin {SPINUP = 0, SPINDOWN = 1};
    enum BandSpin {XUP = 0, XDOWN = 1, YUP = 2, YDOWN = 3};
    typedef ModelParamsDetSDW::BC_Type BC_Type;
    typedef ModelParamsDetSDW::UpdateMethod_Type UpdateMethod_Type;
    typedef ModelParamsDetSDW::SpinProposalMethod_Type SpinProposalMethod_Type;

/*

    Statistics about updates

*/
    uint32_t acceptedGlobalShifts;
    uint32_t attemptedGlobalShifts;
    uint32_t acceptedWolffClusterUpdates;
    uint32_t attemptedWolffClusterUpdates;
    uint32_t acceptedWolffClusterShiftUpdates;
    uint32_t attemptedWolffClusterShiftUpdates;
    num addedWolffClusterSize;

    

/*
    
    Hopping constant representations

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

    Neighbor tables
      
*/
    PeriodicSquareLatticeNearestNeighbors spaceNeigh;
    PeriodicChainNearestNeighbors<1> timeNeigh;



/*

    propK and variations

*/
    checkarray<MatNum, 2> propK;
    MatNum& propKx;
    MatNum& propKy;

    //for shifting green functions to obtain equivalency of symmetric Trotter decomposition
    //[checker board decomposition could be applied alternatively]
    checkarray<MatNum, 2> propK_half;       //factor of -dtau/2 in exponential
    MatNum& propKx_half;
    MatNum& propKy_half;
    checkarray<MatNum, 2> propK_half_inv;   //factor of +dtau/2 in exponential
    MatNum& propKx_half_inv;
    MatNum& propKy_half_inv;



/*

    Green's function

*/    
    MatCpx& g;
    // MatCpx& gFwd;
    // MatCpx& gBwd;


/*

    Auxiliary / bosonic fields

*/    
    //three component sdw-order parameter,
    //column indexes timeslice, row indexes site
    MatNum phi0;
    MatNum phi1;
    MatNum phi2;
    //discrete field for cdwU: l_i(\tau_k) == cdwl(i,k)
    MatInt cdwl;
    //evaluation of element-wise functions of phi and cdwl:
    MatNum coshTermPhi;		// cosh(\lambda \dtau |\phi|)
    MatNum sinhTermPhi; 	// sinh(\lambda \dtau |\phi|) / |\phi|
    MatNum coshTermCDWl;	// cosh(\sqrt(\dtau) cdwU \eta)
    MatNum sinhTermCDWl;	// sinh(\sqrt(\dtau) cdwU \eta)


    
/*

    Adjustment of MC step size

*/   
    num phiDelta;   //MC step size for field components (box update)
    num angleDelta; //  for rotation update (this is the minimal cos\theta)
    num scaleDelta; //  for scaling update (the standard deviation of the gaussian radius update)
    //used to adjust phiDelta/angleDelta/scaleDelta; acceptance ratios for local field updates
    num targetAccRatioLocal_phi;
    num lastAccRatioLocal_phi;
    RunningAverage accRatioLocal_box_RA, accRatioLocal_rotate_RA, accRatioLocal_scale_RA;

    //adjustment of phiDelta, angleDelta, scaleDelta:
    static constexpr num InitialPhiDelta = 0.5;
    static constexpr num InitialAngleDelta = 0.0;
    static constexpr num InitialScaleDelta = 0.1;
    static constexpr num MinScaleDelta = 0.0;
    static constexpr num MaxScaleDelta = 1.0;
    static constexpr num MinAngleDelta = -1.0;
    static constexpr num MaxAngleDelta = 1.0;
    static constexpr uint32_t AccRatioAdjustmentSamples = 100;
    static constexpr num phiDeltaGrowFactor = 1.05;
    static constexpr num phiDeltaShrinkFactor = 0.95;
    //initialized to the respective constants above in the constructor
    num curminAngleDelta;
    num curmaxAngleDelta;
    num curminScaleDelta;
    num curmaxScaleDelta;


/*
    
    internal counter of performed sweeps. This should be serialized

*/
    uint32_t performedSweeps;


/*
    
    Observables:

*/
    num normPhi;        //averaged norm of field			//not very sensible physically
    Phi meanPhi;		//averaged field
    num meanPhiSquared;	//square of averaged field
    num normMeanPhi;	//norm of averaged field
    num sdwSusc;        //spin-density-wave susceptibility

    checkarray<VecNum, 2> kOcc;     //Fermion occupation number in momentum space for x/y-band; site-index: k-vectors
    VecNum& kOccX;
    VecNum& kOccY;
//    checkarray<VecNum, 2> kOccImag;
//    VecNum& kOccXimag;
//    VecNum& kOccYimag;

    checkarray<VecNum, 2> occ;     //Fermion occupation number in Real space for x/y-band; indexed by site
    VecNum& occX;
    VecNum& occY;
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
    GenMat<MatNum, 2> occCorr;  // index by [XBAND,YBAND] etc.
    // charge-charge correlations
    // <n_i n_j> = <n_{x,i}n_{x,i}> + <n_{x,i}n_{y,i}> + <n_{y,i}n_{x,i}> + <n_{y,i}n_{y,i}>
    MatNum chargeCorr;
    // Fourier transforms, index by wave-vector k = (kx,ky)
    // represented as ksite (like kOcc)
    //
    // We only need to write out the measurements in the form of the
    // fourier-transformed version
    GenMat<VecNum, 2> occCorrFT;
    VecNum chargeCorrFT;

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
        result[0] = phi0(site,timeslice);
        result[1] = phi1(site,timeslice);
        result[2] = phi2(site,timeslice);
        return result;
    }


/*
  
    Functions related to auxiliary / bosonic fields, setting up propK,
    setting up the UdV storage and setting up the Green's function
    
*/
    void setupRandomField();
    // return (coshTerm*, sinhTerm*) for these field values
    std::tuple<num,num> getCoshSinhTermPhi(num phi0, num phi1, num phi2);
    std::tuple<num,num> getCoshSinhTermCDWl(int32_t cdwl);

    //compute new entries of coshTerm* and sinhTerm* from current phi0, phi1, phi2, cdwl
    void updateCoshSinhTermsPhi();
    void updateCoshSinhTermsCDWl();
    void updateCoshSinhTerms();
    //the same, for a single site
    void updateCoshSinhTermsPhi(uint32_t site, uint32_t timeslice);
    void updateCoshSinhTermsCDWl(uint32_t site, uint32_t timeslice);
    void updateCoshSinhTerms(uint32_t site, uint32_t timeslice);
    void setupPropK();          //compute e^(-dtau*K..) matrices by diagonalization
    void setupUdVStorage_and_calculateGreen();


    
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
    MatCpx cbLMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    // with A: NxN, sign = +/- 1, band = XBAND|YBAND: return a matrix equal to A * E^(sign * dtau * K_band)
    template <class Matrix>
    MatCpx cbRMultHoppingExp(const Matrix& A, Band band, int sign, bool invertedCbOrder = false);

    //cbLMultHoppingExp and cbRMultHoppingExp need separate implementations for each CheckerboardMethod,
    //this cannot be realized by a direct partial template specialization, but we need to have a proxy
    //with function overloading
    //compare first solution in winning answer at:
    //http://stackoverflow.com/questions/1501357/template-specialization-of-particular-members
    template <class Matrix>
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_NONE>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD_BERG>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    //functions called by the above:
    template<class Matrix>
    void cb_assaad_applyBondFactorsLeft(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver);
    template<class Matrix>
    void cb_assaad_applyBondFactorsRight(Matrix& result, uint32_t subgroup, num ch_hor, num sh_hor, num ch_ver, num sh_ver);

    //the following take a 4Nx4N matrix A and effectively multiply B(k2,k1)
    //or its inverse to the left or right of it and return the result
    MatCpx checkerboardLeftMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardRightMultiplyBmat(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardLeftMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1);
    MatCpx checkerboardRightMultiplyBmatInv(const MatCpx& A, uint32_t k2, uint32_t k1);

    //helpers for the checkerboardMultiplyFunctions
    MatCpx rightMultiplyBk(const MatCpx& orig, uint32_t k);     //multiply B(k,k-1) from right to orig, return result
    MatCpx rightMultiplyBkInv(const MatCpx& orig, uint32_t k);  //multiply B(k,k-1)^-1 from right to orig, return result
    MatCpx leftMultiplyBk(const MatCpx& orig, uint32_t k);      //multiply B(k,k-1) from left to orig, return result
    MatCpx leftMultiplyBkInv(const MatCpx& orig, uint32_t k);   //multiply B(k,k-1)^-1 from left to orig, return result

/*
  
    Checkerboard or dense computation of B-matrix, some helpers for both cases
    
*/
    //compute B-matrix using dense matrix products or checkerboardLeftMultiplyBmat
    MatCpx computeBmatSDW(uint32_t k2, uint32_t k1);            

    // compute sqrt(dtau) * cdwU * eta_{cdwl[i]}
    template<class Vec>
    VecNum compute_d_for_cdwl(const Vec& cdwl);
    // compute sqrt(dtau) * cdwU * eta_{cdwl} for cdwl
    num compute_d_for_cdwl_site(uint32_t cdwl);


    
/*

    Routines for consistency checks / debugging
    
*/
    // Compute e^( sign * dtau * V(phi-configuration) ) as a dense matrix.
    // for one timeslice. This may be useful for checks.
    MatCpx computePotentialExponential(int sign, VecNum phi0, VecNum phi1, VecNum phi2, VecInt cdwl);

    //several checks, many are normally commented out
    virtual void consistencyCheck();



/*
  
    Routines related to measurements
    
*/    
    //the upper function shiftGreenSymmetric() calls the following
    //helper with functors RightMultiply, LeftMultiply depending on
    //the CheckerboardMethod
    template<class RightMultiply, class LeftMultiply>
    MatCpx shiftGreenSymmetric_impl(RightMultiply, LeftMultiply);

    //measuring observables
    void initMeasurements();				//reset stored observable values (beginning of a sweep)
    void measure(uint32_t timeslice);                   //measure observables for one timeslice
    void finishMeasurements();				//finalize stored observable values (end of a sweep)
    std::set<uint32_t> timeslices_included_in_measurement; 	//for a consistency check -- sweep includes correct #timeslices
    // compute the structure factor from a matrix of real space correlations
    void computeStructureFactor(VecNum& out_k, const MatNum& in_r);
    

    
/*

    Monte Carlo updates, some related data structures
    
*/
    //update the auxiliary field and the green function in the single timeslice
    void updateInSlice(uint32_t timeslice);
    //this template member function calls the right one of those below
    //return acceptance ratio for that sweepm
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
    	MatCpx X;
    	MatCpx Y;
    	MatCpx Rj;
    	MatCpx::fixed<4,4> Sj;
    	MatCpx Cj;
    	MatCpx::fixed<4,4> tempBlock;
    	MatCpx::fixed<4,4> Mj;
    	DelayedUpdatesData(uint32_t N, uint32_t delaySteps)
    		: X(4*delaySteps, 4*N), Y(4*N, 4*delaySteps), Rj(4, 4*N), Cj(4*N, 4) {
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
    changedPhiInt proposeNewPhiBox(uint32_t site, uint32_t timeslice);			//from a box
    changedPhiInt proposeRotatedPhi(uint32_t site, uint32_t timeslice);		//same length, new angles
    changedPhiInt proposeScaledPhi(uint32_t site, uint32_t timeslice);			//new length, same angles
    changedPhiInt proposeRotatedScaledPhi(uint32_t site, uint32_t timeslice);	//new length, new angles
    //generate a new cdwl:
    changedPhiInt proposeNewCDWl(uint32_t site, uint32_t timeslice);			//choose Metropolis-randomly from +-1, +-2
    //more generic helpers
    num deltaSPhi(uint32_t site, uint32_t timeslice, Phi newphi);
    MatCpx::fixed<4,4> get_delta_forsite(Phi newphi, int32_t new_cdwl,
    		uint32_t timeslice, uint32_t site);

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
    	MatNum phi0;
    	MatNum phi1;
    	MatNum phi2;
    	MatNum coshTermPhi;
    	MatNum sinhTermPhi;

    	MatCpx g;

    	std::unique_ptr<checkarray<std::vector<UdVV>, 1>> UdVStorage;

    	//for the cluster update: mark visited sites by 1, else 0
    	MatUint visited;		//as usual: row is spatial index, col is time index
    	typedef std::tuple<uint32_t, uint32_t> SpaceTimeIndex;
    	//for the cluster update: keep book about which sites to try to add next
    	std::stack<SpaceTimeIndex> next_sites;
    	GlobalMoveData(uint32_t N, uint32_t m) :
            phi0(N, m+1), phi1(N, m+1), phi2(N, m+1),
            coshTermPhi(N, m+1), sinhTermPhi(N, m+1),
            g(4*N, 4*N),
            UdVStorage(new checkarray<std::vector<UdVV>, 1>),
            visited(N, m+1), next_sites()
        { }
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
public:
/*

    Serialization

*/
    // only functions that can pass the key to these functions have access
    // -- in this way access is granted only to select DetQMC methods
    template<class Archive>
    void saveContents(SerializeContentsKey const &sck, Archive &ar) {
        Base::saveContents(sck, ar);            //base class
        serializeContentsCommon(sck, ar);
    }

    //after loadContents() a sweep must be performed before any measurements are taken:
    //else the green function would not be in a valid state
    template<class Archive>
    void loadContents(SerializeContentsKey const &sck, Archive &ar) {
        Base::loadContents(sck, ar);            //base class
        serializeContentsCommon(sck, ar);
        //the fields now have a valid state, update UdV-storage to start
        //sweeping again
        setupUdVStorage_and_calculateGreen();
        //now: lastSweepDir == SweepDirection::Up --> the next sweep will be downwards
    }

    template<class Archive>
    void serializeContentsCommon(SerializeContentsKey const& sck, Archive& ar) {
    	(void)sck;
    	ar & acceptedGlobalShifts;
    	ar & attemptedGlobalShifts;
    	ar & acceptedWolffClusterUpdates;
    	ar & attemptedWolffClusterUpdates; 
        ar & acceptedWolffClusterShiftUpdates;  
        ar & attemptedWolffClusterShiftUpdates;         
    	ar & addedWolffClusterSize;
        ar & phi0 & phi1 & phi2;
        ar & cdwl;
        ar & coshTermPhi & sinhTermPhi;
        ar & coshTermCDWl & sinhTermCDWl;
        ar & phiDelta & angleDelta & scaleDelta;
        ar & targetAccRatioLocal_phi & lastAccRatioLocal_phi;
        ar & accRatioLocal_box_RA & accRatioLocal_rotate_RA & accRatioLocal_scale_RA;
        ar & curminAngleDelta & curmaxAngleDelta;
        ar & curminScaleDelta & curmaxScaleDelta;
        ar & performedSweeps;
    }
};





/*

	Static helper member functions

*/
template <CheckerboardMethod CBM>
inline std::string DetSDW<CBM>::bandstr(Band b) {
	return (b == XBAND) ? "x" : (b == YBAND ? "y" : "N");
}
template <CheckerboardMethod CBM>
inline std::string DetSDW<CBM>::spinstr(Spin s) {
	return (s == SPINUP) ? "up" : (s == SPINDOWN ? "dn" : "N");
}
template <CheckerboardMethod CBM>
inline std::string DetSDW<CBM>::bandspinstr(BandSpin bs) {
    switch (bs) {
    case XUP:   return "XUp";
    case XDOWN: return "XDown";
    case YUP:   return "YUp";
    case YDOWN: return "YDown";
    default:    return "N";
    }
}
template <CheckerboardMethod CBM>
inline std::string DetSDW<CBM>::updateMethodstr(UpdateMethod_Type um) {
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
template <CheckerboardMethod CBM>
inline std::string DetSDW<CBM>::spinProposalMethodstr(SpinProposalMethod_Type sp) {
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
template <CheckerboardMethod CBM>
inline num DetSDW<CBM>::cdwl_gamma(int32_t l) {
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
template <CheckerboardMethod CBM>
inline num DetSDW<CBM>::cdwl_eta(int32_t l) {
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
template<CheckerboardMethod CBM>
template<typename Callable>
void DetSDW<CBM>::for_each_band(Callable func) {
    func(XBAND);
    func(YBAND);
}
template<CheckerboardMethod CBM>
template<typename Callable> inline
void DetSDW<CBM>::for_each_site(Callable func) {
    for (uint32_t site = 0; site < pars.N; ++site) {
        func(site);
    }
}
template<CheckerboardMethod CBM>
template<typename Callable> inline
void DetSDW<CBM>::for_each_timeslice(Callable func) {
    for (uint32_t k = 1; k <= m; ++k) {
        func(k);
    }
}
// template<typename CallableSiteTimeslice, typename V> inline
// V sumWholeSystem(CallableSiteTimeslice f, V init);

template<CheckerboardMethod CBM>
template<typename CallableSiteTimeslice, typename V> inline
V DetSDW<CBM>::sumWholeSystem(CallableSiteTimeslice f, V init) {
    for (uint32_t timeslice = 1; timeslice <= pars.m; ++timeslice) {
        for (uint32_t site = 0; site < pars.N; ++site) {
            init += f(site, timeslice);
        }
    }
    return init;
}
template<CheckerboardMethod CBM>
template<typename CallableSiteTimeslice, typename V> inline
V DetSDW<CBM>::averageWholeSystem(CallableSiteTimeslice f, V init) {
    V sum = sumWholeSystem(f, init);
    return sum / num(pars.m * pars.N);
}




/*    

      wrappers to use to instantiate template functions of the base class

*/
template<CheckerboardMethod Checkerboard>
struct DetSDW<Checkerboard>::sdwComputeBmat {
    DetSDW<Checkerboard>* parent;
    sdwComputeBmat(DetSDW<Checkerboard>* parent_) :
        parent(parent_)
        { }
    MatCpx operator()(uint32_t gc, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        return parent->computeBmatSDW(k2, k1);
    }
};

template<CheckerboardMethod Checkerboard>
struct DetSDW<Checkerboard>::sdwLeftMultiplyBmat {
    DetSDW<Checkerboard>* parent;
    sdwLeftMultiplyBmat(DetSDW<Checkerboard>* parent_) :
        parent(parent_)
    { }
    MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (Checkerboard != CB_NONE) {
            return parent->checkerboardLeftMultiplyBmat(mat, k2, k1);
        } else {
            return parent->computeBmatSDW(k2, k1) * mat;
        }
    }
};

template<CheckerboardMethod Checkerboard>
struct DetSDW<Checkerboard>::sdwRightMultiplyBmat {
    DetSDW<Checkerboard>* parent;
    sdwRightMultiplyBmat(DetSDW<Checkerboard>* parent_) :
        parent(parent_)
        { }
    MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (Checkerboard != CB_NONE) {
            return parent->checkerboardRightMultiplyBmat(mat, k2, k1);
        } else {
            return mat * parent->computeBmatSDW(k2, k1);
        }
    }
};

template<CheckerboardMethod Checkerboard>
struct DetSDW<Checkerboard>::sdwLeftMultiplyBmatInv {
    DetSDW<Checkerboard>* parent;
    sdwLeftMultiplyBmatInv(DetSDW<Checkerboard>* parent_) :
        parent(parent_)
        { }
    MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (Checkerboard != CB_NONE) {
            return parent->checkerboardLeftMultiplyBmatInv(mat, k2, k1);
        } else {
            return arma::inv(parent->computeBmatSDW(k2, k1)) * mat;
        }
    }
};

template<CheckerboardMethod Checkerboard>
struct DetSDW<Checkerboard>::sdwRightMultiplyBmatInv {
    DetSDW<Checkerboard>* parent;
    sdwRightMultiplyBmatInv(DetSDW<Checkerboard>* parent_) :
        parent(parent_)
        { }
    MatCpx operator()(uint32_t gc, const MatCpx& mat, uint32_t k2, uint32_t k1) {
        (void)gc;
        assert(gc == 0);
        if (Checkerboard != CB_NONE) {
            return parent->checkerboardRightMultiplyBmatInv(mat, k2, k1);
        } else {
            return mat * arma::inv(parent->computeBmatSDW(k2, k1));
        }
    }
};





#endif /* DETSDW_H_ */
