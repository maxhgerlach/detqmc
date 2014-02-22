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
#include "neighbortable.h"
#include "RunningAverage.h"
#include "checkarray.h"
#include "normaldistribution.h"

typedef std::complex<num> cpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<cpx> SpMatCpx;
typedef arma::Col<cpx> VecCpx;
typedef arma::Cube<cpx> CubeCpx;

class SerializeContentsKey;

std::unique_ptr<DetModel> createDetSDW(RngWrapper& rng, ModelParams pars);

enum CheckerboardMethod {
	CB_NONE,				//regular, dense matrix products
	CB_SANTOS,				//checkerboard, four break-ups in e^{K_*} as described in: R. R. dos Santos, Braz. J. Phys 33, 36 (2003).
	CB_ASSAAD,				//checkerboard, two break-ups in e^{K_*} as described in: F. F. Assaad, in Quantum Simulations Complex Many-Body Syst. From Theory to Algorithms, edited by J. Grotendorst, D. Marx, and A. Muramatsu (FZ-J�lich, J�lich, Germany, 2002).
	CB_ASSAAD_BERG,			//checkerboard, two break-ups, making sure all multiplications are symmetric, as described by Erez Berg
};
std::string cbmToString(CheckerboardMethod cbm);


// template parameters: evaluate time-displaced Green functions? do a checker-board decomposition of
// which kind?
template <bool TimeDisplaced, CheckerboardMethod Checkerboard>
class DetSDW: public DetModelGC<1, cpx, TimeDisplaced> {
    DetSDW(RngWrapper& rng, const ModelParams& pars);
public:
    friend std::unique_ptr<DetModel> createDetSDW(RngWrapper& rng, ModelParams pars);
    virtual ~DetSDW();
    virtual uint32_t getSystemN() const;
    virtual MetadataMap prepareModelMetadataMap() const;

    virtual void thermalizationOver();

    virtual void sweep(bool takeMeasurements);
    virtual void sweepThermalization();
    virtual void sweepSimple(bool takeMeasurements);
    virtual void sweepSimpleThermalization();

	//Perform the correction of the Green's function to ensure an
	//effectively symmetric Trotter decomposition.  This returns the shifted matrix for the current timeslice.
	//This should be done before measurements.
    virtual MatCpx shiftGreenSymmetric();
protected:
    typedef DetModelGC<1, cpx, TimeDisplaced> Base;
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

    typedef VecNum::fixed<3> Phi;       //value of the three-component field at a single site and timeslice
    
    //complex 4x4 identity matrix
    const MatCpx::fixed<4,4> eye4cpx;

    RngWrapper& rng;
    NormalDistribution normal_distribution;

    static const uint32_t d = 2;
    static const uint32_t z = 2*d;
    const bool checkerboard;
    const std::string checkerboardMethod;
    const uint32_t L;
    const uint32_t N;
    const num r;
    const num txhor;
    const num txver;
    const num tyhor;
    const num tyver;
    const num mu;
    const num c;
    const num u;
    const num lambda;

    enum Band {XBAND = 0, YBAND = 1};
    static inline std::string bandstr(Band b) {
    	return (b == XBAND) ? "x" : (b == YBAND ? "y" : "N");
    }
    enum Spin {SPINUP = 0, SPINDOWN = 1};
    static inline std::string spinstr(Spin s) {
    	return (s == SPINUP) ? "up" : (s == SPINDOWN ? "dn" : "N");
    }
    enum BC_Type { PBC, APBC_X, APBC_Y, APBC_XY };
    BC_Type bc;
    enum UpdateMethod_Type { ITERATIVE, WOODBURY, DELAYED };
    static inline std::string updateMethodstr(UpdateMethod_Type um) {
    	switch (um) {
    	case ITERATIVE:
    		return "iterative";
    	case WOODBURY:
    		return "woodbury";
    	case DELAYED:
    		return "delayed";
    	default:
    		return "invalid";
    	}
    }
    UpdateMethod_Type updateMethod;
    enum SpinProposalMethod_Type { BOX, ROTATE_THEN_SCALE, ROTATE_AND_SCALE };
    static inline std::string spinProposalMethodstr(SpinProposalMethod_Type sp) {
    	switch (sp) {
    	case BOX:
    		return "box";
    	case ROTATE_THEN_SCALE:
    		return "rotate_then_scale";
    	case ROTATE_AND_SCALE:
    		return "rotate_and_scale";
    	default:
    		return "invalid";
    	}
    }
    SpinProposalMethod_Type spinProposalMethod;
    uint32_t delaySteps;					//for delayed updates

    const bool rescale;
    const uint32_t rescaleInterval;
    const num rescaleGrowthFactor;
    const num rescaleShrinkFactor;
    uint32_t acceptedRescales;
    uint32_t attemptedRescales;

    uint32_t repeatUpdateInSlice;

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
    PeriodicSquareLatticeNearestNeighbors spaceNeigh;
    PeriodicChainNearestNeighbors<1> timeNeigh;

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

    MatCpx& g;
    // MatCpx& gFwd;
    // MatCpx& gBwd;

    //three component sdw-order parameter,
    //column indexes timeslice, row indexes site
    MatNum phi0;
    MatNum phi1;
    MatNum phi2;
    //evaluation of element-wise functions of phi:
    MatNum phiCosh;         // cosh(dtau * |phi|)
    MatNum phiSinh;         // sinh(dtau * |phi|) / |phi|


    num phiDelta;       //MC step size for field components (box update)
    num angleDelta;		//  for rotation update (this is the minimal cos\theta)
    num scaleDelta;		//  for scaling update (the standard deviation of the gaussian radius update)
    //used to adjust phiDelta/angleDelta/scaleDelta; acceptance ratios for local field updates
    num targetAccRatioLocal;
    num lastAccRatioLocal;
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
    bool adaptScaleDelta;

    uint32_t performedSweeps;		//internal counter of performed sweeps. This should be serialized

    //Observables:
    num normPhi;        //averaged norm of field
    Phi meanPhi;		//averaged field
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

    //these for_each functions don't really work well with the class template
    template<typename Callable>
    void for_each_band(Callable func) {
        func(XBAND);
        func(YBAND);
    }
    template<typename Callable>
    void for_each_site(Callable func) {
        for (uint32_t site = 0; site < N; ++site) {
            func(site);
        }
    }
    template<typename Callable>
    void for_each_timeslice(Callable func) {
        for (uint32_t k = 1; k <= m; ++k) {
            func(k);
        }
    }
    template<typename CallableSiteTimeslice, typename V>
    V sumWholeSystem(CallableSiteTimeslice f, V init) {
        for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
            for (uint32_t site = 0; site < N; ++site) {
                init += f(site, timeslice);
            }
        }
        return init;
    }
    template<typename CallableSiteTimeslice, typename V>
    V averageWholeSystem(CallableSiteTimeslice f, V init) {
        V sum = sumWholeSystem(f, init);
        return sum / num(m * N);
    }

    uint32_t coordsToSite(uint32_t x, uint32_t y) {
        return y*L + x;
    }

    void setupRandomPhi();
    void setupPropK();          //compute e^(-dtau*K..) matrices by diagonalization


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
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_SANTOS>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbLMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
    							  const Matrix& A, Band band, int sign, bool invertedCbOrder = false);
    template <class Matrix>
    MatCpx cbRMultHoppingExp_impl(std::integral_constant<CheckerboardMethod, CB_ASSAAD>,
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
    template<class Matrix>
    void cb_santos_applyBondFactorsLeft(Matrix& result, NeighDir neigh, uint32_t subgroup, num ch, num sh);
    template<class Matrix>
    void cb_santos_applyBondFactorsRight(Matrix& result, NeighDir neigh, uint32_t subgroup, num ch, num sh);

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


    MatCpx computeBmatSDW(uint32_t k2, uint32_t k1);            //compute B-matrix using dense matrix products or checkerboardLeftMultiplyBmat

    // Compute e^( sign * dtau * V(phi-configuration) ) as a dense matrix.
    // for one timeslice. This may be useful for checks.
    MatCpx computePotentialExponential(int sign, VecNum phi0, VecNum phi1, VecNum phi2);


    //the upper function calls the following helper with functors RightMultiply, LeftMultiply depending
    //on the CheckerboardMethod
    template<class RightMultiply, class LeftMultiply>
    MatCpx shiftGreenSymmetric_impl(RightMultiply, LeftMultiply);



    //update the auxiliary field and the green function in the single timeslice
    void updateInSlice(uint32_t timeslice);
    //this template member function calls the right one of those below
    template<class CallableProposeSpin>
    void callUpdateInSlice_for_updateMethod(uint32_t timeslice, CallableProposeSpin proposeSpin) {
    	switch(updateMethod) {
    	case ITERATIVE:
    		updateInSlice_iterative(timeslice, proposeSpin);
    		break;
    	case WOODBURY:
    		updateInSlice_woodbury(timeslice, proposeSpin);
    		break;
    	case DELAYED:
    		updateInSlice_delayed(timeslice, proposeSpin);
    		break;
    	}
    };
    //specific implementations of updateInSlice, the callable parameter is for one of the field proposal routines below
    //CallableProposeSpin should have the signature tuple<bool,Phi>(uint32_t site, uint32_t timeslice)
    //the returned bool is false if the returned field needs to be rejected (e.g. negative magnitude)
    template<class CallableProposeSpin>
    void updateInSlice_iterative(uint32_t timeslice, CallableProposeSpin proposeSpin);
    template<class CallableProposeSpin>
    void updateInSlice_woodbury(uint32_t timeslice, CallableProposeSpin proposeSpin);
   	//this one uses a nested struct because for some status
    template<class CallableProposeSpin>
    void updateInSlice_delayed(uint32_t timeslice, CallableProposeSpin proposeSpin);
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
    typedef std::tuple<bool,Phi> boolPhi;
    boolPhi proposeNewField(uint32_t site, uint32_t timeslice);				//from a box
    boolPhi proposeRotatedField(uint32_t site, uint32_t timeslice);			//same length, new angles
    boolPhi proposeScaledField(uint32_t site, uint32_t timeslice);			//new length, same angles
    boolPhi proposeRotatedScaledField(uint32_t site, uint32_t timeslice);	//new length, new angles
    num deltaSPhi(uint32_t site, uint32_t timeslice, Phi newphi);
    MatCpx::fixed<4,4> get_deltanonzero(Phi newphi, uint32_t timeslice, uint32_t site);

    //Try a timeslice-global move, where all the phi-fields of a timeslice are multiplied
    //by a common factor.
    void attemptTimesliceRescaleMove(uint32_t timeslice, num factor);
    num deltaSPhiTimesliceRescale(uint32_t timeslice, num factor);

    //compute the total value of the action associated with the field phi
    num phiAction();


    //measuring observables
    void initMeasurements();				//reset stored observable values (beginning of a sweep)
    void measure(uint32_t timeslice);		//measure observables for one timeslice
    void finishMeasurements();				//finalize stored observable values (end of a sweep)


    // for now: check that phiCosh and phiSinh are set correctly
    virtual void consistencyCheck();

    //wrappers to use to instantiate template functions of the base class
    struct sdwComputeBmat {
        DetSDW<TimeDisplaced,Checkerboard>* parent;
        sdwComputeBmat(DetSDW<TimeDisplaced,Checkerboard>* parent) :
            parent(parent)
        { }
        MatCpx operator()(uint32_t gc, uint32_t k2, uint32_t k1) {
            (void)gc;
            assert(gc == 0);
            return parent->computeBmatSDW(k2, k1);
        }
    };

    struct sdwLeftMultiplyBmat {
        DetSDW<TimeDisplaced,Checkerboard>* parent;
        sdwLeftMultiplyBmat(DetSDW<TimeDisplaced,Checkerboard>* parent) :
            parent(parent)
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

    struct sdwRightMultiplyBmat {
        DetSDW<TimeDisplaced,Checkerboard>* parent;
        sdwRightMultiplyBmat(DetSDW<TimeDisplaced,Checkerboard>* parent) :
            parent(parent)
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

    struct sdwLeftMultiplyBmatInv {
        DetSDW<TimeDisplaced,Checkerboard>* parent;
        sdwLeftMultiplyBmatInv(DetSDW<TimeDisplaced,Checkerboard>* parent) :
            parent(parent)
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

    struct sdwRightMultiplyBmatInv {
        DetSDW<TimeDisplaced,Checkerboard>* parent;
        sdwRightMultiplyBmatInv(DetSDW<TimeDisplaced,Checkerboard>* parent) :
            parent(parent)
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


public:
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
        setupUdVStorage_skeleton(sdwComputeBmat(this));
        //now: lastSweepDir == SweepDirection::Up --> the next sweep will be downwards
    }

    template<class Archive>
    void serializeContentsCommon(SerializeContentsKey const& sck, Archive& ar) {
    	(void)sck;
    	ar & acceptedRescales & attemptedRescales;
        ar & phi0 & phi1 & phi2;
        ar & phiCosh & phiSinh;
        ar & phiDelta & angleDelta & scaleDelta;
        ar & targetAccRatioLocal & lastAccRatioLocal;
        ar & accRatioLocal_box_RA & accRatioLocal_rotate_RA & accRatioLocal_scale_RA;
        ar & curminAngleDelta & curmaxAngleDelta;
        ar & curminScaleDelta & curmaxScaleDelta;
        ar & performedSweeps;
    }
};



#endif /* DETSDW_H_ */
