#ifndef DETSDWPARAMS_H
#define DETSDWPARAMS_H

#include "detmodelparams.h"

enum CheckerboardMethod {
    CB_NONE,                //regular, dense matrix products
    CB_ASSAAD_BERG, //checkerboard, two break-ups, making sure all multiplications are symmetric, as described by Erez Berg
};



//use this in place of
//  ModelParams<DetSDW<CB_NONE>> and ModelParams<DetSDW<CB_ASSAAD_BERG>> 


struct ModelParamsDetSDW {
    std::string model;          // should be "sdw"

    bool turnoffFermions;      // normally false, if true: simulate a pure O(opdim) model, without fermion determinants

    bool turnoffFermionMeasurements; // normally false. If turnoffFermions is true, but this is false, we simulate a model
                                     // with fermions, but don't compute observables from the Green's function.

    bool dumpGreensFunction;    // dump the various blocks of the Green's function in real space when measuring.  Defaults to false, use very sparingly!
    
    bool checkerboard;               //use a checkerboard decomposition for computing the propagator
    std::string updateMethod_string; //"iterative", "woodbury", or "delayed"
    enum UpdateMethod_Type { ITERATIVE, WOODBURY, DELAYED };
    UpdateMethod_Type updateMethod;
    std::string spinProposalMethod_string;  //"box", "rotate_then_scale", or "rotate_and_scale"
    enum SpinProposalMethod_Type { BOX, ROTATE_THEN_SCALE, ROTATE_AND_SCALE };
    SpinProposalMethod_Type spinProposalMethod;
    bool adaptScaleVariance;         //valid unless spinProposalMethod=="box" -- this controls if the variance of the spin updates should be adapted during thermalization
    uint32_t delaySteps;             //parameter in case updateMethod is "delayed"

    bool overRelaxation;        // for no-fermion simulations:  additional non-ergodic, microcanonical over relaxation sweeps
    uint32_t repeatOverRelaxation;
    std::string repeatOverRelaxation_string; // how many overRelaxationSweeps to carry out in a row during a single sweep; only if the fermions are turned off.  Pass \"systemSize\" to use a number growing with system size: N * beta/dtau <-- this is way too much, though.  Default: 1

    uint32_t opdim;             // order parameter dimension: 1, 2 or 3 (default: 3)
    bool phi2bosons;            // if this is true: run calculations with a simple theory (r/2)\sum_i \phi_i^2 -- ignores parameter u and spatial terms in the bosonic action
    bool phiFixed;              // if this is true: set a constant, fixed phi field -- all phi_i(tau) equal to (1 [, 0[, 0]])
    num r;
    num c;                      // currently fixed to 1.0
    num u;
    num lambda;	//fermion-boson coupling strength
    num txhor;	//hopping constants depending on direction and band
    num txver;	
    num tyhor;	
    num tyver;	
    num cdwU;	  //hoping to get a CDW transition
    num mu;       // chemical potential for both orbitals
    num mux, muy; // *alternatively*: specify separately for the two
                  // orbitals, if not specified we will set:
                  //    mux := muy := mu .
                  // Otherwise mux and muy supersede mu.
    bool weakZflux;             // apply the weakest possible magnetic field (+- 1/L^2) in z-direction.
    uint32_t L;
    uint32_t N;                 // L*L, set in check()
    uint32_t d;                 // should be 2
    num beta;
    uint32_t m;     //either specify number of timeslices 'm'
    num dtau;       //or timeslice separation 'dtau'
    uint32_t s;     //separation of timeslices where the Green function is calculated
                    //from scratch
    num accRatio;   //target acceptance ratio for tuning spin update box size

    std::string bc_string; //boundary conditions: "pbc", "apbc-x", "apbc-y" or "apbc-xy"
    enum BC_Type { PBC, APBC_X, APBC_Y, APBC_XY };
    BC_Type bc;

    uint32_t globalUpdateInterval; //attempt global move every # sweeps
    bool globalShift;              //perform a global constant shift move?
    bool wolffClusterUpdate;       //perform a Wolff single cluster update?
    bool wolffClusterShiftUpdate;  // perform a combined global constant shift and Wolff single cluster update

    uint32_t repeatWolffPerSweep;
    std::string repeatWolffPerSweep_string; // how many Wolff cluster flips we attempt in a row during a single sweep; this is mostly useful if the fermions are turned off, as the fermion determinant is only taken into consideration after the whole series of updates.    Pass \"systemSize\" to use a number growing with system size: N * beta/dtau <-- this is way too much, though.  Default: 1

    uint32_t repeatUpdateInSlice;  //how often to repeat updateInSlice for eacht timeslice per sweep, default: 1

    std::set<std::string> specified;

    ModelParamsDetSDW() :
        model("sdw"), turnoffFermions(false), turnoffFermionMeasurements(false),
        dumpGreensFunction(false),
        checkerboard(),
        updateMethod_string("woodbury"), updateMethod(WOODBURY),
        spinProposalMethod_string("box"), spinProposalMethod(BOX),
        adaptScaleVariance(), delaySteps(), overRelaxation(false),
        repeatOverRelaxation(1), repeatOverRelaxation_string(""),
        opdim(3), phi2bosons(false), phiFixed(false), r(), c(1.0), u(1.0), lambda(),
        txhor(), txver(), tyhor(), tyver(), cdwU(), mu(), mux(0.), muy(0.), weakZflux(false), L(), N(), d(2),
        beta(), m(), dtau(), s(), accRatio(), bc_string("pbc"), bc(PBC), globalUpdateInterval(),
        globalShift(), wolffClusterUpdate(), wolffClusterShiftUpdate(), repeatWolffPerSweep(1), repeatWolffPerSweep_string(""),
        repeatUpdateInSlice(),
        specified()
    { }

    void check();
    MetadataMap prepareMetadataMap() const;

    void set_exchange_parameter_value(num val) {
        r = val;
        specified.insert("r");
    }
    num  get_exchange_parameter_value() { return r; }    
private:
    friend class boost::serialization::access;

    template<class Archive>
        void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar  & model & turnoffFermions & turnoffFermionMeasurements
            & dumpGreensFunction
            & checkerboard
            & updateMethod_string & updateMethod
            & spinProposalMethod_string & spinProposalMethod
            & adaptScaleVariance & delaySteps
            & overRelaxation
            & repeatOverRelaxation & repeatOverRelaxation_string
            & opdim & phi2bosons & phiFixed & r & c & u & lambda
            & txhor & txver & tyhor & tyver
            & cdwU & mu & mux & muy & weakZflux
            & L & N & d & beta & m & dtau & s & accRatio
            & bc_string & bc
            & globalUpdateInterval & globalShift & wolffClusterUpdate
            & wolffClusterShiftUpdate
            & repeatWolffPerSweep & repeatWolffPerSweep_string
            & repeatUpdateInSlice
            & specified;
    }
    inline std::string updateMethodstr(UpdateMethod_Type um) const {
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
    inline std::string spinProposalMethodstr(SpinProposalMethod_Type sp) const {
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
    
};


#endif /* DETSDWPARAMS_H */
