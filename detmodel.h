/*
 * detmodel.h
 *
 *  Created on: Feb 18, 2013
 *      Author: gerlach
 */

#ifndef DETMODEL_H_
#define DETMODEL_H_


#include <functional>
#include <utility>
#include <memory>               // unique_ptr
#include <algorithm>            // minmax_element
#include <vector>
#include <tuple>
#include <armadillo>
#include <cassert>
#include <type_traits>          // std::is_same

#include "tools.h"
#include "toolsdebug.h"
#include "rngwrapper.h"
#include "checkarray.h"
#include "dataserieswritersucc.h"
#include "detmodelparams.h"
#include "detmodelloggingparams.h"
#include "observable.h"
#include "udv.h"
#include "metadata.h"
#include "timing.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/base_object.hpp"
#pragma GCC diagnostic pop
#include "boost_serialize_array.h"
#include "boost_serialize_armadillo.h"

typedef double num;
typedef std::complex<double> cpx;

typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;
typedef arma::Col<int32_t> VecInt;
typedef arma::Mat<int32_t> MatInt;
typedef arma::Col<uint32_t> VecUint;
typedef arma::Mat<uint32_t> MatUint;
typedef arma::Col<cpx> VecCpx;
typedef arma::Mat<cpx> MatCpx;
typedef arma::SpMat<num> SpMatNum;
typedef std::tuple<MatNum,MatNum,MatNum,MatNum> MatNum4;
typedef UdV<num> UdVnum;



// The previous function template, which was to be specialized
// explicitly, led to errors that were hard to understand and fix.
// To avoid this we do not use any generic template at all and
// just *require* that model implementations are accompanied by
// a function
//     void createReplica(std::unique_ptr<Model>& replica_out,
//                        RngWrapper& rng, ModelParams pars)



// // This template should not be used.  Provide explicit specializations
// // for each implementation of a model
// template<class Model, class Params = ModelParams<Model> >
// std::unique_ptr<Model> createReplica(RngWrapper& rng, Params pars) {
//     BOOST_MPL_ASSERT(( not_defined<Model> ));
//     pars.check();
//     return std::unique_ptr<Model>(new Model(rng, pars));
// }


// Compute min{1, exp(-Delta), the acceptance probability for a
// replica exchange in parallel tempering
// --
// currently only specialized for DetSDW -- it may be that this will
// have to be generalized and redefined for other models where the
// relevant terms in the action are not simply
//   parameter * action_contributio
template<class Model>
num get_replica_exchange_probability(
    num parameter_1, num action_contribution_1,
    num parameter_2, num action_contribution_2)
{
    (void)parameter_1;
    (void)parameter_2;
    (void)action_contribution_1;
    (void)action_contribution_2;    
    BOOST_MPL_ASSERT(( not_defined<Model> ));
    return 0.0;
}





//base class for a model to be simulated by determinantal quantum Monte Carlo
//


//purely abstract base class
class DetModel {
public:
    virtual ~DetModel() { };
    virtual uint32_t getSystemN() const = 0;

    //Create a MetadataMap describing the parameters of the
    //simulated model
    virtual MetadataMap prepareModelMetadataMap() const = 0;

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
    //These values are valid after a sweep, in which measurements where taken
    virtual std::vector<ScalarObservable> getScalarObservables() = 0;
    virtual std::vector<VectorObservable> getVectorObservables() = 0;
    virtual std::vector<KeyValueObservable> getKeyValueObservables() = 0;

    //perform a sweep updating the auxiliary field with costly re-computations
    //of Green functions from scratch
    //if takeMeasurements==true: perform measurements of all observables
    virtual void sweepSimple(bool takeMeasurements) = 0;
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepSimpleThermalization() = 0;


    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    //if takeMeasurements==true: perform measurements of all observables
    virtual void sweep(bool takeMeasurements) = 0;
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters
    virtual void sweepThermalization() = 0;

    //notify DetModel that thermalization has finished
    //do nothing by default
    virtual void thermalizationOver() {
    }
public:
    // For serialization. To be called by DetQMC methods
    template<class Archive>
    void saveContents(Archive &) {
    }
    template<class Archive>
    void loadContents(Archive &) {
    }
};



//GreenComponents is the number of independent sectors of the Green's function,
//e.g. in the S=1/2-Hubbard model it is 2 for spin up and spin down
//
//ValueType can be a complex number if the Green function is not purely real
//
//if TimeDisplaced==true: generate code that evaluates time-displaced green
//functions in the sweep
//
//This provides template functions like sweep_skeleton<>() that expect callable template
//arguments for routines that compute B-matrices etc.  A derived class that provides
//them should instantiate them and call them in their implementations of virtual
//functions like sweep()

template<uint32_t GreenComponents, typename ValueType = num, bool TimeDisplaced = false>
class DetModelGC: public DetModel {
protected:
    template<class ModelParams>
    DetModelGC(const ModelParams& pars, uint32_t greenComponentSize,
               const DetModelLoggingParams& loggingParams = DetModelLoggingParams());
public:
    virtual ~DetModelGC()
    { }

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
    virtual std::vector<ScalarObservable> getScalarObservables();
    virtual std::vector<VectorObservable> getVectorObservables();
    virtual std::vector<KeyValueObservable> getKeyValueObservables();

    //perform a sweep updating the auxiliary field with costly re-computations
    //of Green functions from scratch
    //
    //  if takeMeasurements == true : perform observable measurements
    //
    //  Callable_GC_k2_k1: take arguments green component, timeslices k2 > k1,
    //  and give the corresponding B-matrix
    //  Callable_UpdateInSlice_k: argument timeslice k, update fields for this timeslice
    //  Callable_init: no arguments, init observable measurements for this sweep
    //  Callable_measure_k: argument timeselice k = 1,...,m, take measurement data for timeslice k
    //  Callable_finish: no arguments, finalize observable measurements for this sweep
    template<class Callable_GC_k2_k1, class Callable_UpdateInSlice_k,
    		 class Callable_init, class Callable_measure_k,
    		 class Callable_finish >
    void sweepSimple_skeleton(bool takeMeasurements,
                              Callable_GC_k2_k1 computeBmat,
                              Callable_UpdateInSlice_k updateInSlice,
                              Callable_init initMeasurement,
                              Callable_measure_k measure, Callable_finish finishMeasurement);
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters, but does not take any measurements ever
    template<class Callable_GC_k2_k1, class Callable_UpdateInSlice_k>
    void sweepSimpleThermalization_skeleton(Callable_GC_k2_k1 computeBmat,
    		Callable_UpdateInSlice_k updateInSlice);


    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //  /* at some point Will give equal-time and time-displaced Green functions if TimeDisplace == ture */.
    //if takeMeasurements == true : perform observable measurements
     //
    //*_Callable_GC_mat_k2_k1: take arguments green-component, some matrix,
    //                         time slices k2 > k1
    //      -> return left/right product of matrix with Bmat or Bmat-inverse
    //      useful if a checkerboard-breakup is performed
    //Callable_UpdateInSlice_k: argument timeslice k, update fields for this timeslice
    //Callable_init: no arguments, init observable measurements for this sweep
    //Callable_measure_k: argument timeselice k, take measurement data for timeslice k
    //Callable_finish: no arguments, finalize observable measurements for this sweep
    //
    //optional:
    //  Callable_GlobalUpdate: no arguments, this is called before each sweep-down
    //          and may be used to provide a global update encompassing all
    //          imaginary time slices
    //          By default: do nothing
    //optional:
    //  Callable_GreenConsistency:
    //    Arguments: const Mat& g1, const Mat& g2, SweepDirection cursweepdir.
    //    Perform some sort of consistency check, comparing the two matrices,
    //    also be informed about whether we are sweeping up or down.
    //    By default: do nothing
    template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
             class c_Callable_GC_mat_k2_k1, class d_Callable_GC_mat_k2_k1,
             class Callable_UpdateInSlice_k,
             class Callable_init, class Callable_measure_k, class Callable_finish,
             class Callable_GlobalUpdate = VoidNoOp,
             class Callable_GreenConsistency = VoidNoOp>
    void sweep_skeleton(bool takeMeasurements,
    					a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
                        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
                        c_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
                        d_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
                        Callable_UpdateInSlice_k updateInSlice,
                        Callable_init initMeasurement, Callable_measure_k measure,
                        Callable_finish finishMeasurement,
                        Callable_GlobalUpdate globalUpdate = VoidNoOp(),
                        Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());
    //the same to be called during thermalization, may do the same or iteratively
    //adjust parameters, but does not take any measurements ever
    template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
             class c_Callable_GC_mat_k2_k1, class d_Callable_GC_mat_k2_k1,
             class Callable_UpdateInSlice_k,
             class Callable_GlobalUpdate = VoidNoOp,
             class Callable_GreenConsistency = VoidNoOp>
    void sweepThermalization_skeleton(
                        a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
                        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
                        c_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
                        d_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
                        Callable_UpdateInSlice_k updateInSlice,
                        Callable_GlobalUpdate globalUpdate = VoidNoOp(),
                        Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());
protected:
    typedef arma::Mat<ValueType> MatV;
    typedef arma::Col<ValueType> VecV;
    typedef arma::Cube<ValueType> CubeV;
    typedef UdV<ValueType> UdVV;
    typedef std::tuple<MatV,MatV,MatV,MatV> MatV4;

//    //update the auxiliary field and the green function in the single timeslice
//    virtual void updateInSlice(uint32_t timeslice) = 0;
//    //separate function to be called during thermalization, by default just do the
//    //same; a derived class may override this to introduce an adaptive behavior
//    virtual void updateInSliceThermalization(uint32_t timeslice) {
//        updateInSlice(timeslice);
//    }

    //Given B(beta, tau) = V_l d_l U_l and B(tau, 0) = U_r d_r V_r
    //calculate a tuple of four NxN matrices (a,b,c,d) with
    // a = G(0), b = -(1-G(0))*B^(-1)(tau,0), c = B(tau,0)*G(0), d = G(tau)
    //b is the backward time-displaced Green function; c the forward time-
    //displaced Green function; d is the equal-time Green function
    //todo: get rid of this MatNum4 business
    MatV4 greenFromUdV_timedisplaced(const UdVV& UdV_l, const UdVV& UdV_r) const;
    //use a faster method that does not yield information about the time-displaced
    //Green functions.
    // Uses B(beta, tau) = U_l d_l V_l   and     B(tau, 0) = U_r d_r V_r,
    // computes G(tau) = [Id + B(tau,0).B(beta,tau)]^{-1}
    //                 = [Id + U_r d_r V_r U_l d_l V_l]^{-1}
    //                 = (V_t_L V_t_x) D_x^{-1} (U_R U_x)^{dagger}
        // and stores the singular values of G^{-1} [their product yields
    // the absolute value of the inverse determinant of G]
    void greenFromUdV(MatV& green_out, VecNum& green_inv_sv, const UdVV& UdV_l, const UdVV& UdV_r) const;
    //The following is useful to compute G(\beta) = [1 + B(\beta, 0)]^{-1}
    void greenFromEye_and_UdV(MatV& green_out, VecNum& green_inv_sv, const UdVV& UdV_r) const;

    //compute Green function from UdV-decomposed matrices L/R
    //for a single timeslice and update the member variables green --
    //and if desired -- greenFwd and greenBwd.
    //Also updates green_inv_sv.
    void updateGreenFunctionUdV(uint32_t gc, const UdVV& UdV_L, const UdVV& UdV_R);
    void updateGreenFunction_Eye_UdV(uint32_t gc, const UdVV& UdV_R);    

    //for each greenComponent call a function with the greenComponent as a parameter
    template<typename Callable>
    void for_each_gc(Callable func) {
        for (uint32_t gc = 0; gc < GreenComponents; ++gc) {
            func(gc);
        }
    }

    // //call in a derived class:
    // //  Callable_GC_k2_k1: take arguments green component, timeslices k2 > k1,
    // //  and give the corresponding B-matrix
    // //
    // //This will setup the UdV storage used to compute Green's functions from scratch
    // //in the following sweep-down and also compute the Green's function G(\beta)
    // template<class Callable_GC_k2_k1>
    // void setupUdVStorage_and_calculateGreen_skeleton(Callable_GC_k2_k1 computeBmat);


    //call in a derived class:
    //Callable_GC_mat_k2_k1: take arguments green-component, some matrix,
    //                       time slices k2 > k1
    //      -> return left product of matrix with Bmat: Bmat(k2,k1) * matrix
    //
    //This will setup the UdV storage used to compute Green's functions from scratch
    //in the following sweep-down and also compute the Green's function G(\beta)
    template<class Callable_GC_mat_k2_k1>
    void setupUdVStorage_and_calculateGreen_skeleton(Callable_GC_mat_k2_k1 leftMultiplyBmat);
    // this is the same, but computes the Green's function G(k \dtau)
    // at an arbitrary timeslice k.  Note that this leaves the UdV
    // storage etc in a stage that is unsuitable for a continued
    // sweep.  Use this only for consistency checks etc and restore
    // the proper state before continuing the regular program flow.
    template<class Callable_GC_mat_k2_k1>
    void setupUdVStorage_and_calculateGreen_forTimeslice_skeleton(uint32_t timeslice,
                                                                  Callable_GC_mat_k2_k1 leftMultiplyBmat);
    


    //helpers for sweep_skeleton(), sweepThermalization_skeleton():
    //
    //Callable_GC_mat_k2_k1: take arguments green-component, some matrix,
    //                       time slices k2 > k1
    //      -> return left/right product of matrix with Bmat or Bmat-inverse
    //optional:
    //  Callable_GreenConsistency:
    //    Arguments: const Mat& g1, const Mat& g2, SweepDirection cursweepdir.
    //    Perform some sort of consistency check, comparing the two matrices,
    //    also be informed about whether we are sweeping up or down.
    //    By default: do nothing

    template<class Callable_GC_mat_k2_k1,
             class Callable_GreenConsistency = VoidNoOp>
    void advanceDownGreen(Callable_GC_mat_k2_k1 rightMultiplyBmat,
                          uint32_t l, uint32_t gc,
                          Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());

    template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1>
    void wrapDownGreen(a_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
                       b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
                       uint32_t k, uint32_t gc);

    template<class Callable_GC_mat_k2_k1,
             class Callable_GreenConsistency = VoidNoOp>
    void advanceUpGreen(Callable_GC_mat_k2_k1 leftMultiplyBmat,
                        uint32_t l, uint32_t gc,
                        Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());

//    template<class Callable_GC_mat_k2_k1>
//    void advanceUpUpdateStorage(Callable_GC_mat_k2_k1 leftMultiplyBmat,
//                                uint32_t l, uint32_t gc);

    template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1>
    void wrapUpGreen(a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
                     b_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
                     uint32_t k, uint32_t gc);

    //these receive as a template parameter the function to call for updates in a slice,
    //as well as B-Mat multiplicators like above
    template <class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
              class CallableUpdateInSlice,
              class Callable_init, class Callable_measure_k, class Callable_finish,
              class Callable_GreenConsistency = VoidNoOp>
    void sweepUp(bool takeMeasurements,
                 a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
                 b_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
                 CallableUpdateInSlice funcUpdateInSlice,
                 Callable_init initMeasurement, Callable_measure_k measure,
                 Callable_finish finishMeasurement,
                 Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());
    template <class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
              class CallableUpdateInSlice,
              class Callable_init, class Callable_measure_k, class Callable_finish,
              class Callable_GreenConsistency = VoidNoOp>
    void sweepDown(bool takeMeasurements,
                   a_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
                   b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
                   CallableUpdateInSlice funcUpdateInSlice,
                   Callable_init initMeasurement, Callable_measure_k measure,
                   Callable_finish finishMeasurement,
                   Callable_GreenConsistency greenConsistencyCheck = VoidNoOp());

    // This method is called after each sweep.
    // A derived class, which implements the model, may overload it to check
    // its internal state for consistency.  An exception should be thrown if
    // it fails.
    virtual void consistencyCheck() {
    	// default: do nothing
    }


    //Green component size, e.g. sz == N for the Hubbard model
    const uint32_t sz;

    //some simulation parameters are already relevant for member functions implemented
    //in this base class, the rest will only be used in derived classes
    const bool timedisplaced;
    const num beta;     //inverse temperature
    const uint32_t m;   //number of imaginary time discretization steps (time slices) beta*m=dtau
    const uint32_t s;   //maximum interval between time slices where the Green-function is calculated from scratch
    const uint32_t n;   //number of time slices where the Green-function is calculated from scratch == ceil(m/s)
    const num dtau;     // beta / m

    // this struct contains parameters related to logging that should
    // be done in this class
    DetModelLoggingParams loggingParams;
    std::unique_ptr<DoubleVectorWriterSuccessive> svLogging, svMaxLogging, svMinLogging;


//    //equal-imaginary-time and time-displaced Green's functions
//    //slices indexed k=0..m correspond to time slices at dtau*k,
//    //which are then indexed by sites in row and column.
//    //Most code, however, only uses timeslices k >= 1 ! Don't rely on g*.slice(0)
//    //being valid.
//    //The Green functions for k=0 are conceptually equal to those for k=m.
//    checkarray<CubeV, GreenComponents> green;
//    checkarray<CubeV, GreenComponents> greenFwd;
//    checkarray<CubeV, GreenComponents> greenBwd;

    // During the sweep: hold the matrix elements of the equal-time Green's function for the
    // current timeslice
    checkarray<MatV, GreenComponents> green;
    uint32_t currentTimeslice;					//currently green is valid for this timeslice
    // This stores the singular values of G^{-1} (assuming that det(G) > 0).  This is only valid after
    // updateGreenFunction[_Eye_]UdV
    checkarray<VecNum, GreenComponents> green_inv_sv;
    

    //The UdV-instances in UdVStorage will not move around after setup, so storing
    //the (rather big) objects in the vector is fine.
    //However, for instance for deciding on doing a global update we need the possibility
    //to swap the whole vector of UdV's. For this reason: handle the whole container over
    //a unique_ptr.
    //Remember that to recover the decomposed matrices: m == U * diag(d) * trans(V)
    //The conjugate-transpose still needs to be taken.
    const UdVV eye_UdV;
    const MatV eye_gc;
    std::unique_ptr<checkarray<std::vector<UdVV>, GreenComponents>> UdVStorage;

    enum class SweepDirection: int {Up = 1, Down = -1};
    SweepDirection lastSweepDir;

    //observable handling -- these contain information about observables (such as their names)
    //as well as reference to their current value, which will be shared with simulation management
    //in a different class. The values referenced there are to be updated here in the replica class.
    std::vector<ScalarObservable> obsScalar;
    std::vector<VectorObservable> obsVector;
    std::vector<KeyValueObservable> obsKeyValue;

public:
    // serialization by DetQMC::serializeContents
    template<class Archive>
    void saveContents(Archive &ar) {
        DetModel::saveContents(ar);        //base class

        //the following commented lines are for contents we no longer
        //serialize as they can be reconstructed from the field configuration
        //easily with setupUdVstorage and a sweep
//      ar & green & greenFwd & greenBwd;
//      ar & UdVStorage;
//      ar & lastSweepDir;
    }

    template<class Archive>
    void loadContents(Archive &ar) {
        DetModel::loadContents(ar);        //base class
        //UdV-storage, green, green_inv_sv, greenFwd, greenBwd still need to be recast into a valid state
        //by a derived class!
        //TODO: this is a mess!
    }
};


template<uint32_t GC, typename V, bool TimeDisplaced>
template<class ModelParams>
DetModelGC<GC,V,TimeDisplaced>::DetModelGC(const ModelParams& pars, uint32_t greenComponentSize,
                                           const DetModelLoggingParams& loggingParams_ /* default argument */) :
    sz(greenComponentSize),
    timedisplaced(TimeDisplaced),
    beta(pars.beta), m(pars.m), s(pars.s),
    n(uint32_t(std::ceil(double(m) / s))),
    dtau(pars.dtau),
    loggingParams(loggingParams_),
    svLogging(), svMaxLogging(), svMinLogging(),
    green(), //greenFwd(), greenBwd(),
    currentTimeslice(),
    green_inv_sv(),
    eye_UdV(sz), eye_gc(arma::eye<MatV>(sz, sz)),
    UdVStorage(new checkarray<std::vector<UdVV>, GC>),
    lastSweepDir(SweepDirection::Up),
    obsScalar(), obsVector(), obsKeyValue()
{
    //init Green's functions with zeros
    for(uint32_t gc = 0; gc < GC; ++gc) {
        green[gc].zeros(greenComponentSize, greenComponentSize);
        green_inv_sv[gc].zeros(greenComponentSize);            
//        if (TimeDisplaced) {
//            greenFwd[gc].zeros(greenComponentSize, greenComponentSize, m+1);
//            greenBwd[gc].zeros(greenComponentSize, greenComponentSize, m+1);
//        }
    }

    if (loggingParams.logSV) {
        svLogging = std::unique_ptr<DoubleVectorWriterSuccessive>(
            new DoubleVectorWriterSuccessive(
                loggingParams.logSV_filename,
                false // append to file = false: always start a new file for this
                )
            );
        svLogging->addHeaderText("Attention: this file is recreated and the log restarted for each run of the program. It is not continued if the simulation is resumed from a saved state.");
        svLogging->addHeaderText("Here we log the logarithmic range of the singular values: log(max sv) - log(min sv), each time the (inverse) Green's function is computed from the singular values.");
        svLogging->writeHeader();
        svMaxLogging = std::unique_ptr<DoubleVectorWriterSuccessive>(
            new DoubleVectorWriterSuccessive(
                loggingParams.logSV_max_filename,
                false // append to file = false: always start a new file for this
                )
            );
        svMaxLogging->addHeaderText("Attention: this file is recreated and the log restarted for each run of the program. It is not continued if the simulation is resumed from a saved state.");
        svMaxLogging->addHeaderText("Here we log the max logarithmic singular values: log(max sv), each time the (inverse) Green's function is computed from the singular values.");
        svMaxLogging->writeHeader();
        svMinLogging = std::unique_ptr<DoubleVectorWriterSuccessive>(
            new DoubleVectorWriterSuccessive(
                loggingParams.logSV_min_filename,
                false // append to file = false: always start a new file for this
                )
            );
        svMinLogging->addHeaderText("Attention: this file is recreated and the log restarted for each run of the program. It is not continued if the simulation is resumed from a saved state.");
        svMinLogging->addHeaderText("Here we log the min logarithmic singular values: log(min sv), each time the (inverse) Green's function is computed from the singular values.");
        svMinLogging->writeHeader();
    }

//  // Default functors for multiplication with B-matrices
//  for_each_gc( [this](uint32_t gc) {
//      leftMultiplyBmat[gc] = [this, gc](const MatV A, uint32_t k2, uint32_t k1) -> MatV {
//          return computeBmat[gc](k2, k1) * A;
//      };
//      rightMultiplyBmat[gc] = [this, gc](const MatV A, uint32_t k2, uint32_t k1) -> MatV {
//          return A * computeBmat[gc](k2, k1);
//      };
//      leftMultiplyBmatInv[gc] = [this, gc](const MatV A, uint32_t k2, uint32_t k1) -> MatV {
//          return arma::inv(computeBmat[gc](k2, k1)) * A;
//      };
//      rightMultiplyBmatInv[gc] = [this, gc](const MatV A, uint32_t k2, uint32_t k1) -> MatV {
//          return A * arma::inv(computeBmat[gc](k2, k1));
//      };
//  } );
}

template<uint32_t GC, typename V, bool TimeDisplaced>
std::vector<ScalarObservable> DetModelGC<GC,V,TimeDisplaced>::getScalarObservables() {
    return obsScalar;
}

template<uint32_t GC, typename V, bool TimeDisplaced>
std::vector<VectorObservable> DetModelGC<GC,V,TimeDisplaced>::getVectorObservables() {
    return obsVector;
}

template<uint32_t GC, typename V, bool TimeDisplaced>
std::vector<KeyValueObservable> DetModelGC<GC,V,TimeDisplaced>::getKeyValueObservables() {
    return obsKeyValue;
}

template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_mat_k2_k1>
void DetModelGC<GC,V,TimeDisplaced>::setupUdVStorage_and_calculateGreen_forTimeslice_skeleton(
    uint32_t timeslice, Callable_GC_mat_k2_k1 leftMultiplyBmat) {

    // handle the timeslice == beta case separately
    if (timeslice == m) {
        setupUdVStorage_and_calculateGreen_skeleton(leftMultiplyBmat);
        return;
    }

    
    timing.start("setupUdVStorage");

    auto setup = [this, timeslice, &leftMultiplyBmat](uint32_t gc) -> uint32_t {
        std::vector<UdVV>& storage = (*UdVStorage)[gc];
        storage = std::vector<UdVV>(n + 1);

        uint32_t lk = uint32_t(std::floor(num(timeslice) / s));
        uint32_t k_lkp1 = ((lk < n - 1) ? (s*(lk+1)) : (m));

        // std::cout << "timeslice: " << timeslice << " lk: " << lk << " k_lkp1: " << k_lkp1 << "\n";
        
        // std::cout << "(" << k_lkp1 << ", " << timeslice << ")\n";
        udvDecompose(storage[0], leftMultiplyBmat(gc, eye_gc, k_lkp1, timeslice));

        uint32_t storageCounter = 0;
        for (uint32_t l = lk + 1; l <= n - 1; ++l) {
            // std::cout << "l = " << l << ", storageCounter = " << storageCounter << "\n";
            const MatV&   U_l   = storage[storageCounter].U;
            const VecNum& d_l   = storage[storageCounter].d;
            const MatV&   V_t_l = storage[storageCounter].V_t;
            const uint32_t k_l   = s*l;
            const uint32_t k_lp1 = ((l < n - 1) ? (s*(l+1)) : (m));
            // std::cout << "(" << k_lp1 << ", " << k_l << ")\n";
            MatV B_lp1_times_U_l = leftMultiplyBmat(gc, U_l, k_lp1, k_l);
            udvDecompose<V>(storage[storageCounter+1], B_lp1_times_U_l * arma::diagmat(d_l));
            storage[storageCounter+1].V_t =  V_t_l * storage[storageCounter+1].V_t;
            ++storageCounter;
        }

        uint32_t target = ( (lk * s == timeslice) ? lk - 1 : lk );
        // std::cout << "target = (lk * s ==  timeslice) ? lk - 1 : lk ... "
        // << target << " = ( (" << lk << " == " << s << " * " << timeslice << ") ? "
        // << lk - 1 << " : " << lk << " );\n";

        for (uint32_t l = 0; l <= target; ++l) {
            // std::cout << "l = " << l << ", storageCounter = " << storageCounter << "\n";
            const MatV&   U_l   = storage[storageCounter].U;
            const VecNum& d_l   = storage[storageCounter].d;
            const MatV&   V_t_l = storage[storageCounter].V_t;
            const uint32_t k_l   = s*l;
            const uint32_t k_lp1 = ((l < lk) ? (s*(l+1)) : (timeslice));
            // std::cout << "(" << k_lp1 << ", " << k_l << ")\n";
            MatV B_lp1_times_U_l = leftMultiplyBmat(gc, U_l, k_lp1, k_l);
            udvDecompose<V>(storage[storageCounter+1], B_lp1_times_U_l * arma::diagmat(d_l));
            storage[storageCounter+1].V_t =  V_t_l * storage[storageCounter+1].V_t;            
            ++storageCounter;
        }

        return storageCounter;  // return the highest index in the storage
    };

    for (uint32_t gc = 0; gc < GC; ++gc) {
        uint32_t index = setup(gc);
        updateGreenFunction_Eye_UdV(gc, (*UdVStorage)[gc][index]);
    }
    
    timing.stop("setupUdVStorage");
}



template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_mat_k2_k1>
void DetModelGC<GC,V,TimeDisplaced>::setupUdVStorage_and_calculateGreen_skeleton(
        Callable_GC_mat_k2_k1 leftMultiplyBmat) {
    timing.start("setupUdVStorage");
    auto setup = [this, &leftMultiplyBmat](uint32_t gc) {
        std::vector<UdVV>& storage = (*UdVStorage)[gc];
        storage = std::vector<UdVV>(n + 1);

        storage[0] = eye_UdV; 
        // storage[1] = udvDecompose(computeBmat(gc, s, 0));
        udvDecompose(storage[1], leftMultiplyBmat(gc, eye_gc, s, 0));

        for (uint32_t l = 1; l <= n - 1; ++l) {
            const MatV&   U_l   = storage[l].U;
            const VecNum& d_l   = storage[l].d;
            const MatV&   V_t_l = storage[l].V_t;
            const uint32_t k_l   = s*l;
            const uint32_t k_lp1 = ((l < n - 1) ? (s*(l+1)) : (m));
//            MatV B_lp1_times_U_l = computeBmat(gc, k_lp1, k_l) * U_l;
            MatV B_lp1_times_U_l = leftMultiplyBmat(gc, U_l, k_lp1, k_l);
            udvDecompose<V>(storage[l+1], B_lp1_times_U_l * arma::diagmat(d_l));
            storage[l+1].V_t =  V_t_l * storage[l+1].V_t;
        }
    };

    for_each_gc(setup);

    for (uint32_t gc = 0; gc < GC; ++gc) {
        updateGreenFunction_Eye_UdV(gc, (*UdVStorage)[gc][n]);
    }
    currentTimeslice = m;

    lastSweepDir = SweepDirection::Up;
    timing.stop("setupUdVStorage");
}



//warning: the thermalization version below is almost a copy of this -- without measurements
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_k2_k1, class Callable_UpdateInSlice_k,
		 class Callable_init, class Callable_measure_k, class Callable_finish>
void DetModelGC<GC,V,TimeDisplaced>::sweepSimple_skeleton(
		bool takeMeasurements,
		Callable_GC_k2_k1 computeBmat, Callable_UpdateInSlice_k updateInSlice,
		Callable_init initMeasurement,
		Callable_measure_k measure, Callable_finish finishMeasurement) {
    if (takeMeasurements) {
        initMeasurement();
    }
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for_each_gc( [this, timeslice, &computeBmat](uint32_t gc) {
                green[gc] = arma::inv(arma::eye(sz,sz) + computeBmat(gc, timeslice, 0) *
                                      computeBmat(gc, m, timeslice));
            });
        updateInSlice(timeslice);
        if (takeMeasurements) {
            measure(timeslice);
        }
    }
    if (takeMeasurements) {
    	finishMeasurement();
    }
    consistencyCheck();
}

//warning: this is almost a copy of sweepSimple() defined above
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_k2_k1, class Callable_UpdateInSlice_k>
void DetModelGC<GC,V,TimeDisplaced>::sweepSimpleThermalization_skeleton(
        Callable_GC_k2_k1 computeBmat, Callable_UpdateInSlice_k updateInSliceThermalization) {
    for (uint32_t timeslice = 1; timeslice <= m; ++timeslice) {
        for_each_gc( [this, timeslice, &computeBmat](uint32_t gc) {
            green[gc] =
                    arma::inv(arma::eye(sz,sz) + computeBmat(gc, timeslice, 0) *
                                                  computeBmat(gc, m, timeslice));
        });
        updateInSliceThermalization(timeslice);
    }
}



//use a faster method that does not yield information about the time-displaced
//Green functions.
// Uses B(beta, tau) = U_l d_l V_l   and     B(tau, 0) = U_r d_r V_r,
// computes G(tau) = [Id + B(tau,0).B(beta,tau)]^{-1}
//                 = [Id + U_r d_r V_r U_l d_l V_l]^{-1}
//                 = (V_t_L V_t_x) D_x^{-1} (U_R U_x)^{dagger}
template<uint32_t GC, typename V, bool TimeDisplaced>
void DetModelGC<GC,V,TimeDisplaced>::greenFromUdV(
		MatV& green_out,
                VecNum& green_inv_sv,
		const UdVV& UdV_l,
		const UdVV& UdV_r) const {
    timing.start("greenFromUdV");
    const MatV&   U_l   = UdV_l.U;
    const VecNum& d_l   = UdV_l.d;
    const MatV&   V_t_l = UdV_l.V_t;
    const MatV&   U_r   = UdV_r.U;
    const VecNum& d_r   = UdV_r.d;
    const MatV&   V_t_r = UdV_r.V_t;

    using arma::diagmat; using arma::trans;

    MatV VU_rl_product = trans(V_t_r) * U_l;
    MatV UtVt_rl_product = trans(U_r) * V_t_l;


    // here we get just the singular values of G^{-1}
    MatV U_temp, V_t_temp;
    udvDecompose<V>(U_temp, green_inv_sv, V_t_temp,
        UtVt_rl_product +
        diagmat(d_r) * VU_rl_product * diagmat(d_l)
        );

    
    if (loggingParams.logSV) {
        auto min_max_pair = std::minmax_element(green_inv_sv.begin(), green_inv_sv.end());
        num log_min_sv = std::log(*min_max_pair.first);
        num log_max_sv = std::log(*min_max_pair.second);
        // for(auto& i : green_inv_sv)
        //     std::cout << i << ' ';
        // std::cout << '\n';

        svLogging->writeData( log_max_sv - log_min_sv );
        svMinLogging->writeData ( log_min_sv );
        svMaxLogging->writeData ( log_max_sv );        
    }
    
    
    MatV Vt_product = V_t_l * V_t_temp;
    MatV U_product  = U_r * U_temp;
    
    green_out = Vt_product *
        diagmat(1.0 / green_inv_sv) *
        trans(U_product);

    timing.stop("greenFromUdV");
}



template<uint32_t GC, typename V, bool TimeDisplaced>
void DetModelGC<GC,V,TimeDisplaced>::greenFromEye_and_UdV(
		MatV& green_out,
                VecNum& green_inv_sv,
		const UdVV& UdV_r) const {
    timing.start("greenFromUdV");
    //Here we consider the special case U_l*d_l*V_t_l.t() = 1
    const MatV&   U_r   = UdV_r.U;
    const VecNum& d_r   = UdV_r.d;
    const MatV&   V_t_r = UdV_r.V_t;

    using arma::diagmat; using arma::trans;

    // here we get just the singular values of G^{-1}
    MatV U_temp, V_t_temp;
    udvDecompose<V>(U_temp, green_inv_sv, V_t_temp,
        trans(U_r) * V_t_r + diagmat(d_r)
        );


    if (loggingParams.logSV) {
        auto min_max_pair = std::minmax_element(green_inv_sv.begin(), green_inv_sv.end());
        num log_min_sv = std::log(*min_max_pair.first);
        num log_max_sv = std::log(*min_max_pair.second); 
        svLogging->writeData( log_max_sv - log_min_sv );        
        svMinLogging->writeData ( log_min_sv );
        svMaxLogging->writeData ( log_max_sv );        
    }
    

    MatV V_t_product = V_t_r * V_t_temp;
    MatV U_product = U_r * U_temp;
    
    green_out = V_t_product *
    		diagmat(1.0 / green_inv_sv) *
                trans(U_product);

    timing.stop("greenFromUdV");
}



template<uint32_t GC, typename V, bool TimeDisplaced>
typename DetModelGC<GC,V,TimeDisplaced>::MatV4 DetModelGC<GC,V,TimeDisplaced>::greenFromUdV_timedisplaced(
        const UdVV& UdV_l, const UdVV& UdV_r) const {
    timing.start("greenFromUdV_timedisplaced");

    //Ul vs Vl to be compatible with labeling in the notes
    const MatV&   Ul = UdV_l.V;   //!
    const VecNum& dl = UdV_l.d;
    const MatV&   Vl = UdV_l.U;   //!
    const MatV&   Ur = UdV_r.U;
    const VecNum& dr = UdV_r.d;
    const MatV&   Vr = UdV_r.V;

    uint32_t sz_ = Ul.n_rows;

    //submatrix view helpers for 2*N x 2*N matrices
    auto upleft = [sz_](MatV& m) {
        return m.submat(0,0, sz_-1,sz_-1);
    };
    auto upright = [sz_](MatV& m) {
        return m.submat(0,sz_, sz_-1,2*sz_-1);
    };
    auto downleft = [sz_](MatV& m) {
        return m.submat(sz_,0, 2*sz_-1,sz_-1);
    };
    auto downright = [sz_](MatV& m) {
        return m.submat(sz_,sz_, 2*sz_-1,2*sz_-1);
    };

    MatV temp(2*sz_,2*sz_);
    upleft(temp)    = arma::inv(Vr * Vl);
    upright(temp)   = arma::diagmat(dl);
    downleft(temp)  = arma::diagmat(-dr);
    downright(temp) = arma::inv(Ul * Ur);
    UdVV tempUdV = udvDecompose<V>(temp);

    MatV left(2*sz_,2*sz_);
    upleft(left) = arma::inv(Vr);
    upright(left).zeros();
    downleft(left).zeros();
    downright(left) = arma::inv(Ul);

    MatV right(2*sz_,2*sz_);
    upleft(right) = arma::inv(Vl);
    upright(right).zeros();
    downleft(right).zeros();
    downright(right) = arma::inv(Ur);

    MatV result = (left * arma::inv(tempUdV.V)) * arma::diagmat(1.0 / tempUdV.d)
                    * (arma::inv(tempUdV.U) * right);

    timing.stop("greenFromUdV_timedisplaced");

    return MatV4(upleft(result), upright(result),
                   downleft(result), downright(result));
}


template<uint32_t GC, typename V, bool TimeDisplaced>
void DetModelGC<GC,V,TimeDisplaced>::updateGreenFunctionUdV(
        uint32_t gc, const UdVV& UdV_L, const UdVV& UdV_R)
{
    if (TimeDisplaced) {
//        std::tie(std::ignore, greenBwd[gc].slice(targetSlice),
//                greenFwd[gc].slice(targetSlice), green[gc].slice(targetSlice))
//            = greenFromUdV_timedisplaced(UdV_L, UdV_R);
    } else {
        greenFromUdV(green[gc], green_inv_sv[gc], UdV_L, UdV_R);
    }
}

template<uint32_t GC, typename V, bool TimeDisplaced>
void DetModelGC<GC,V,TimeDisplaced>::updateGreenFunction_Eye_UdV(
    uint32_t gc, const UdVV& UdV_R) {
    if (TimeDisplaced) {
        // no-op
    } else {
        greenFromEye_and_UdV(green[gc], green_inv_sv[gc], UdV_R);
    }
}

//compute the green function in timeslice s*(l-1) from scratch with the help
//of the B-matrices computed before in the last up-sweep
//
//preconditions: storage[l] contains B(beta, l*s*dtau)
//               storage[l - 1] contains B((l-1)*s*dtau, 0)
//
//postconditions: storage[l - 1] contains B(beta, (l-1)*s*dtau)
//                green is G((l-1)*s*dtau)
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_mat_k2_k1,
         class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::advanceDownGreen(
        Callable_GC_mat_k2_k1 rightMultiplyBmat,
        uint32_t l, uint32_t gc,
        Callable_GreenConsistency greenConsistencyCheck)
{
    timing.start("advanceDownGreen");

    //This is the point where the function should be called in the
    //sweep, even though we do not actually use green explicitly here.
    //The sweep is now set-up in such a way, that we do one
    //superfluous wrap-up step, this means the advance-functions serve
    //as a refresh of the current time slice.
    assert(currentTimeslice == s*(l-1));

    std::vector<UdVV>& storage = (*UdVStorage)[gc];

    const uint32_t k_l   = ((l < n) ? (s*l) : (m));
    const uint32_t k_lm1 = s*(l-1);

    //UdV_L will correspond to B(beta,k_lm1*dtau)
    UdVV UdV_L;
    if (l < n) {
        //U_l, d_l, V_l correspond to B(beta,k_l*dtau) [set in the last step]
        const MatV&   U_l   = storage[l].U;
        const VecNum& d_l   = storage[l].d;
        const MatV&   V_t_l = storage[l].V_t;

        udvDecompose<V>(UdV_L,
            arma::diagmat(d_l) *
            rightMultiplyBmat(gc, trans(V_t_l), k_l, k_lm1)
            );
        UdV_L.U = U_l * UdV_L.U;
    } else {
        // special case l==n, can compute UdV_L from scratch
        udvDecompose<V>(UdV_L, rightMultiplyBmat(gc, eye_gc, k_l, k_lm1));
    }

    // //Accuracy check:
    MatV g_wrapped;
    if ( not std::is_same<Callable_GreenConsistency, VoidNoOp>::value ) {
        g_wrapped = green[gc];
    }

    if (l - 1 > 0) {
        //UdV_R corresponds to B(k_lm1*dtau,0) [set in last sweep]
        const UdVV& UdV_R = storage[l - 1];

        updateGreenFunctionUdV(gc, UdV_L, UdV_R);
    } else {
        updateGreenFunction_Eye_UdV(gc, UdV_L);
    }
    
    //Accuracy check:
    // print_matrix_rel_diff(g_wrapped, green[gc], "Adv-Down" + numToString(currentTimeslice));
    greenConsistencyCheck(g_wrapped, green[gc], SweepDirection::Down);

    storage[l - 1] = UdV_L;

    currentTimeslice = s*(l-1);

    timing.stop("advanceDownGreen");
}

////compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
////also compute time-displaced Green functions
//template<uint32_t GC, typename V, bool TimeDisplaced>
//void DetModelGC<GC,V,TimeDisplaced>::wrapDownGreen_timedisplaced(uint32_t k, uint32_t gc) {
//  timing.start("wrapDownGreen_timedisplaced");
////    MatV B_k = computeBmat[greenComponent](k, k - 1);
////    MatV B_k_inv = computeBmatInverse[greenComponent](k, k - 1);
////    green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
//
//  green[gc].slice(k - 1) =
//          leftMultiplyBmatInv[gc](rightMultiplyBmat[gc](green[gc].slice(k), k, k-1), k, k-1);
//  //DEBUG: avoid chaining
////    MatV intermed = rightMultiplyBmat[gc](green[gc].slice(k), k, k-1);
////    green[gc].slice(k - 1) = leftMultiplyBmatInv[gc](intermed, k, k-1);
//
////    greenFwd[greenComponent].slice(k - 1) = arma::inv(B_k) * greenFwd[greenComponent].slice(k);
//  greenFwd[gc].slice(k - 1) = leftMultiplyBmatInv[gc](greenFwd[gc].slice(k), k, k-1);
//  greenBwd[gc].slice(k - 1) = rightMultiplyBmat[gc](greenBwd[gc].slice(k), k, k-1);
//  timing.stop("wrapDownGreen_timedisplaced");
//}
//
////compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
////only equal-time Green functions
//template<uint32_t GC, typename V, bool TimeDisplaced>
//void DetModelGC<GC,V,TimeDisplaced>::wrapDownGreen(uint32_t k, uint32_t gc) {
//  timing.start("wrapDownGreen");
////    MatV B_k = computeBmat[greenComponent](k, k - 1);
////    green[greenComponent].slice(k - 1) = arma::inv(B_k) * green[greenComponent].slice(k) * B_k;
//
//  green[gc].slice(k - 1) =
//          leftMultiplyBmatInv[gc](rightMultiplyBmat[gc](green[gc].slice(k), k, k-1), k, k-1);
//  //DEBUG: avoid chaining
////    MatV slice = green[gc].slice(k);
////    MatV intermed = rightMultiplyBmat[gc](slice, k, k-1);
//  //DEBUG: explicit without checkerboard
////    MatV Bmat = computeBmat[gc](k, k-1);
////    MatV intermed = slice * Bmat;
////    green[gc].slice(k - 1) = leftMultiplyBmatInv[gc](intermed, k, k-1);
//
//  timing.stop("wrapDownGreen");
//}

// compute the green function at k-1 by wrapping the one at k (accumulates rounding errors),
// store the result in green.
// [if required also compute time-displaced green functions]
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1>
void DetModelGC<GC,V,TimeDisplaced>::wrapDownGreen(
        a_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
        uint32_t k, uint32_t gc)
{
    timing.start("wrapDownGreen");

    assert(currentTimeslice == k);
    // //DEBUG
    // CHECK_NAN(green[0]);
    // MatV mat1 = rightMultiplyBmat(gc, green[gc], k, k-1);
    // CHECK_NAN(mat1);        
    // green[gc] = leftMultiplyBmatInv(gc, mat1,
    //                                 k, k-1);
    // CHECK_NAN(green[gc]);    
    // //DEBUG

    //ORIG
    green[gc] = leftMultiplyBmatInv(gc, rightMultiplyBmat(gc, green[gc], k, k-1),
    								k, k-1);
                                                                
    if (TimeDisplaced) {
//        greenFwd[gc].slice(k - 1) = leftMultiplyBmatInv(gc, greenFwd[gc].slice(k), k, k-1);
//        greenBwd[gc].slice(k - 1) = rightMultiplyBmat(gc, greenBwd[gc].slice(k), k, k-1);
    }

    currentTimeslice = k-1;

    timing.stop("wrapDownGreen");
}


//update the green function in timeslice s*(l+1) [or m] from scratch with the help
//of B-matrices computed before
//
// preconditions: storage[l + 1] contains B(beta, (l+1)*s*dtau)		[or B(\beta,\beta)]
//                storage[l] contains B(l*s*dtau, 0)
//
// postconditions: storage[l + 1] contains B((l+1)*s*dtau, 0)
//                 green is G((l+1)*s*dtau)
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class Callable_GC_mat_k2_k1,
         class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::advanceUpGreen(
        Callable_GC_mat_k2_k1 leftMultiplyBmat,
        uint32_t l, uint32_t gc,
        Callable_GreenConsistency greenConsistencyCheck)
{
    timing.start("advanceUpGreen");

    std::vector<UdVV>& storage = (*UdVStorage)[gc];

    const uint32_t k_l = s*l;
    const uint32_t k_lp1 = ((l < n - 1) ? (s*(l+1)) : (m));

    //This is the point where the function should be called in the
    //sweep, even though we do not actually use green explicitly here.
    //The sweep is now set-up in such a way, that we do one
    //superfluous wrap-up step, this means the advance-functions serve
    //as a refresh of the current time slice.
    assert(currentTimeslice == k_lp1);

    // //Accuracy check:
    MatV g_wrapped;
    if ( not std::is_same<Callable_GreenConsistency, VoidNoOp>::value ) {
        g_wrapped = green[gc];
    }
    
    //from the last step the following are B(k_l*dtau, 0):
    const MatV&   U_l   = storage[l].U;
    const VecNum& d_l   = storage[l].d;
    const MatV&   V_t_l = storage[l].V_t;

    //UdV_temp will be the new B(k_lp1*dtau, 0):
    UdVV UdV_temp;
    udvDecompose<V>(UdV_temp,
                    leftMultiplyBmat(gc, U_l, k_lp1, k_l) * arma::diagmat(d_l));
    UdV_temp.V_t = V_t_l * UdV_temp.V_t;
    
    if (k_lp1 != m) {
        //The following is B(beta, k_lp1*dtau), valid from the last sweep
        const UdVV& UdV_lp1 = storage[l + 1];

        updateGreenFunctionUdV(gc, UdV_lp1, UdV_temp);
    } else {
        updateGreenFunction_Eye_UdV(gc, UdV_temp);
    }

    storage[l + 1] = UdV_temp;
    
    //Accuracy check:
    // print_matrix_rel_diff(g_wrapped, green[gc], "Adv-Up" + numToString(currentTimeslice));
    greenConsistencyCheck(g_wrapped, green[gc], SweepDirection::Up);

    currentTimeslice = k_lp1;

    timing.stop("advanceUpGreen");
}

////Given B(l*s*dtau, 0) from the last step in the storage, compute
////B([(l+1)*s or m]*dtau, 0) and put it into storage
//template<uint32_t GC, typename V, bool TimeDisplaced>
//template<class Callable_GC_mat_k2_k1>
//void DetModelGC<GC,V,TimeDisplaced>::advanceUpUpdateStorage(
//        Callable_GC_mat_k2_k1 leftMultiplyBmat,
//        uint32_t l, uint32_t gc)
//{
//    timing.start("advanceUpUpdateStorage");
//
//    std::vector<UdVV>& storage = (*UdVStorage)[gc];
//
//    const uint32_t k_l = s*l;
//    const uint32_t k_lp1 = ((l < n - 1) ? (s*(l+1)) : (m));
//
//    //from the last step the following are B(k_l*dtau, 0):
//    const MatV& U_l = storage[l].U;
//    const VecV& d_l = storage[l].d;
//    const MatV& V_l = storage[l].V;
//    //the new B(k_lp1*dtau, 0):
//    storage[l+1] = udvDecompose<V>(leftMultiplyBmat(gc, U_l, k_lp1, k_l)
//                                   * arma::diagmat(d_l));
//
//    storage[l+1].V *= V_l;
//
//    timing.stop("advanceUpUpdateStorage");
//};

////compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
////also handle the time-displaced Green functions
//template<uint32_t GC, typename V, bool TimeDisplaced>
//void DetModelGC<GC,V,TimeDisplaced>::wrapUpGreen_timedisplaced(uint32_t k, uint32_t gc) {
//  timing.start("wrapUpGreen_timedisplaced");
////    MatV B_kp1 = computeBmat[greenComponent](k + 1, k);
////    green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
//
//  green[gc].slice(k + 1) =
//           leftMultiplyBmat[gc](rightMultiplyBmatInv[gc](green[gc].slice(k), k+1, k), k+1, k);
//  //DEBUG: avoid chaining
////    MatV intermed = rightMultiplyBmatInv[gc](green[gc].slice(k), k+1, k);
////    green[gc].slice(k + 1) = leftMultiplyBmat[gc](intermed, k+1, k);
//
////    greenFwd[greenComponent].slice(k + 1) = B_kp1 * greenFwd[greenComponent].slice(k);
//  greenFwd[gc].slice(k + 1) = leftMultiplyBmat[gc](greenFwd[gc].slice(k), k+1, k);
////    greenBwd[greenComponent].slice(k + 1) = greenBwd[greenComponent].slice(k) * arma::inv(B_kp1);
//  greenBwd[gc].slice(k + 1) = rightMultiplyBmatInv[gc](greenBwd[gc].slice(k), k+1, k);
//  timing.stop("wrapUpGreen_timedisplaced");
//}
//
////compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
////only compute equal-time Green functions
//template<uint32_t GC, typename V, bool TimeDisplaced>
//void DetModelGC<GC,V,TimeDisplaced>::wrapUpGreen(uint32_t k, uint32_t gc) {
//  timing.start("wrapUpGreen");
////    MatV B_kp1 = computeBmat[greenComponent](k + 1, k);
////    green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);
//
//  green[gc].slice(k + 1) =
//          leftMultiplyBmat[gc](rightMultiplyBmatInv[gc](green[gc].slice(k), k+1, k), k+1, k);
//  //DEBUG: avoid chaining
////    MatV intermed = rightMultiplyBmatInv[gc](green[gc].slice(k), k+1, k);
////    green[gc].slice(k + 1) = leftMultiplyBmat[gc](intermed, k+1, k);
//
//  timing.stop("wrapUpGreen");
//}

//compute the green function at k+1 by wrapping the one at k (accumulates rounding errors),
//store the result in green.
//if necessary, also handle the time-displaced Green functions
template<uint32_t GC, typename V, bool TimeDisplaced>
template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1>
void DetModelGC<GC,V,TimeDisplaced>::wrapUpGreen(
        a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
        uint32_t k, uint32_t gc)
{
    timing.start("wrapUpGreen");
//  MatV B_kp1 = computeBmat[greenComponent](k + 1, k);
//  green[greenComponent].slice(k + 1) = B_kp1 * green[greenComponent].slice(k) * arma::inv(B_kp1);

    assert(currentTimeslice == k);

    green[gc] = leftMultiplyBmat(gc, rightMultiplyBmatInv(gc, green[gc], k+1, k), k+1, k);

    if (TimeDisplaced) {
//        //  greenFwd[greenComponent].slice(k + 1) = B_kp1 * greenFwd[greenComponent].slice(k);
//        greenFwd[gc].slice(k + 1) = leftMultiplyBmat(gc, greenFwd[gc].slice(k), k+1, k);
//        //  greenBwd[greenComponent].slice(k + 1) = greenBwd[greenComponent].slice(k) * arma::inv(B_kp1);
//        greenBwd[gc].slice(k + 1) = rightMultiplyBmatInv(gc, greenBwd[gc].slice(k), k+1, k);
    }

    currentTimeslice = k + 1;

    timing.stop("wrapUpGreen");
}

template<uint32_t GC, typename V, bool TimeDisplaced>
template <class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
	  class CallableUpdateInSlice,
	  class Callable_init, class Callable_measure_k, class Callable_finish,
          class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::sweepUp(
        bool takeMeasurements,
        a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
        CallableUpdateInSlice funcUpdateInSlice,
        Callable_init initMeasurement, Callable_measure_k measure,
        Callable_finish finishMeasurement,
        Callable_GreenConsistency greenConsistencyCheck)
{
    if (takeMeasurements) {
        initMeasurement();
    }

    auto updateInSliceAndMaybeMeasure = [&,this](uint32_t timeslice) -> void {
        assert(currentTimeslice == timeslice);
        funcUpdateInSlice(timeslice);
        if (takeMeasurements) {
            measure(timeslice);
        }
    };

    //Precondition for the following:
    //We need to have computed the Green function for time slice k=0 so that the first
    //wrap-up step is correct [done by sweepDown],
    //and our storage contains (UdV)_l = B(beta, l*s*dtau) for l = 0, ..., n

    //set storage at k=0 to unity for the up-coming sweep:
    for (uint32_t gc = 0; gc < GC; ++gc) {
        (*UdVStorage)[gc][0] = eye_UdV;
    }
    //sweep up:
    for (uint32_t l = 0; l <= n - 2; ++l) {
        for (uint32_t k = l*s + 1; k <= (l+1)*s; ++k) { // s wrap-up steps
            for (uint32_t gc = 0; gc < GC; ++gc) {
                wrapUpGreen(leftMultiplyBmat, rightMultiplyBmatInv, k - 1, gc);
            }
            updateInSliceAndMaybeMeasure(k);
        }
        for (uint32_t gc = 0; gc < GC; ++gc) {
            advanceUpGreen(leftMultiplyBmat, l, gc, greenConsistencyCheck);        //new version
        }
    }
    //special handling for the highest time-slices
    for (uint32_t k = (n - 1)*s + 1; k <= m; ++k) {
        for (uint32_t gc = 0; gc < GC; ++gc) {
            wrapUpGreen(leftMultiplyBmat, rightMultiplyBmatInv, k - 1, gc);
        }
        updateInSliceAndMaybeMeasure(k);
    }
    //refresh Green's function at highest time-slice
    for (uint32_t gc = 0; gc < GC; ++gc) {
        advanceUpGreen(leftMultiplyBmat, n - 1, gc, greenConsistencyCheck);        //new version
    }

    if (takeMeasurements) {
        finishMeasurement();
    }

    consistencyCheck();
}


template<uint32_t GC, typename V, bool TimeDisplaced>
template <class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
          class CallableUpdateInSlice,
          class Callable_init, class Callable_measure_k, class Callable_finish,
          class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::sweepDown(
        bool takeMeasurements,
        a_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
        CallableUpdateInSlice funcUpdateInSlice,
        Callable_init initMeasurement, Callable_measure_k measure,
        Callable_finish finishMeasurement,
        Callable_GreenConsistency greenConsistencyCheck)
{
    if (takeMeasurements) {
        initMeasurement();
    }

    auto updateInSliceAndMaybeMeasure = [&,this](uint32_t timeslice) -> void {
        assert(currentTimeslice == timeslice);
        // //DEBUG
        // CHECK_NAN(green[0]);
        // //DEBUG
        funcUpdateInSlice(timeslice);
        if (takeMeasurements) {
            measure(timeslice);
        }
    };

    //Precondition for the following:
    //We need to have computed the Green function for time slice k=m (at tau=beta)
    //[this is done by setupUdVStorage_and_calculateGreen_skeleton as well as by
    //sweepUp], and our storage contains (UdV)_l = B(l*s*dtau, 0) for l = 1, ..., n

    // Handle timeslices between l=n (-> k=m) and k=s*(n-1) + 1,
    // in contrast to the lower values of l, these may be less than s
    for (uint32_t k = m; k >= (n-1)*s + 1; --k) {
        updateInSliceAndMaybeMeasure(k);
        for (uint32_t gc = 0; gc < GC; ++gc) {
            wrapDownGreen(leftMultiplyBmatInv, rightMultiplyBmat, k, gc);
        }
    }

    // Handle the remaining timeslices, including advanceDown steps to refresh
    // the Green's function

    // for (uint32_t gc = 0; gc < GC; ++gc) {
    //     (*UdVStorage)[gc][n] = eye_UdV;
    // }

    for (uint32_t l = n - 1; l >= 1; --l) {
        for (uint32_t gc = 0; gc < GC; ++gc) {
            advanceDownGreen(rightMultiplyBmat, l + 1, gc, greenConsistencyCheck);
        }
        for (uint32_t k = l*s; k >= (l-1)*s + 1; --k) {
            updateInSliceAndMaybeMeasure(k);
            for (uint32_t gc = 0; gc < GC; ++gc) {
                wrapDownGreen(leftMultiplyBmatInv, rightMultiplyBmat, k, gc);
            }
        }
    }
    // refresh the Green's function at k=0 so we are ready to sweep-up
    for (uint32_t gc = 0; gc < GC; ++gc) {
        advanceDownGreen(rightMultiplyBmat, 1, gc, greenConsistencyCheck);
    }

    if (takeMeasurements) {
        finishMeasurement();
    }

    consistencyCheck();
}

template<uint32_t GC, typename V, bool TimeDisplaced>
template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
         class c_Callable_GC_mat_k2_k1, class d_Callable_GC_mat_k2_k1,
         class Callable_UpdateInSlice_k,
         class Callable_init, class Callable_measure_k, class Callable_finish,
         class Callable_GlobalUpdate,
         class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::sweep_skeleton(
        bool takeMeasurements,
        a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
        c_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
        d_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
        Callable_UpdateInSlice_k updateInSlice,
        Callable_init initMeasurement, Callable_measure_k measure,
        Callable_finish finishMeasurement,
        Callable_GlobalUpdate globalUpdate,
        Callable_GreenConsistency greenConsistencyCheck)
{
    timing.start("sweep");

    if (lastSweepDir == SweepDirection::Up) {
        globalUpdate();
        sweepDown(takeMeasurements,
                  leftMultiplyBmatInv, rightMultiplyBmat,
                  updateInSlice,
                  initMeasurement, measure, finishMeasurement,
                  greenConsistencyCheck);
        lastSweepDir = SweepDirection::Down;
    } else if (lastSweepDir == SweepDirection::Down) {
        sweepUp(takeMeasurements,
                leftMultiplyBmat, rightMultiplyBmatInv,
                updateInSlice,
                initMeasurement, measure, finishMeasurement,
                greenConsistencyCheck);
        lastSweepDir = SweepDirection::Up;
    }

    timing.stop("sweep");
}

template<uint32_t GC, typename V, bool TimeDisplaced>
template<class a_Callable_GC_mat_k2_k1, class b_Callable_GC_mat_k2_k1,
         class c_Callable_GC_mat_k2_k1, class d_Callable_GC_mat_k2_k1,
         class Callable_UpdateInSlice_k, class Callable_GlobalUpdate,
         class Callable_GreenConsistency>
void DetModelGC<GC,V,TimeDisplaced>::sweepThermalization_skeleton(
        a_Callable_GC_mat_k2_k1 leftMultiplyBmat,
        b_Callable_GC_mat_k2_k1 rightMultiplyBmat,
        c_Callable_GC_mat_k2_k1 leftMultiplyBmatInv,
        d_Callable_GC_mat_k2_k1 rightMultiplyBmatInv,
        Callable_UpdateInSlice_k updateInSliceThermalization,
        Callable_GlobalUpdate globalUpdate,
        Callable_GreenConsistency greenConsistencyCheck)
{
    timing.start("sweep");

    if (lastSweepDir == SweepDirection::Up) {
        globalUpdate();
        sweepDown(false,                                    //no measurements
                  leftMultiplyBmatInv, rightMultiplyBmat,
                  updateInSliceThermalization,
                  VoidNoOp(), VoidNoOp(), VoidNoOp(),       //no measurements
                  greenConsistencyCheck
                  );
        lastSweepDir = SweepDirection::Down;
    } else if (lastSweepDir == SweepDirection::Down) {
        sweepUp(false,
                leftMultiplyBmat, rightMultiplyBmatInv,     //no measurements
                updateInSliceThermalization,
                VoidNoOp(), VoidNoOp(), VoidNoOp(),         //no measurements
                greenConsistencyCheck
                );
        lastSweepDir = SweepDirection::Up;
    }

    timing.stop("sweep");
}




//compute e^{-scalar * real matrix}, matrix must be symmetric
MatNum computePropagator(num scalar, const MatNum& matrix);
//compute e^{-scalar * compex matrix}, matrix must be hermitian
MatCpx computePropagator(num scalar, const MatCpx& matrix);




#endif /* DETMODEL_H_ */
