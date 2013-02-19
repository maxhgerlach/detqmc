/*
 * detmodel.h
 *
 *  Created on: Feb 18, 2013
 *      Author: gerlach
 */

#ifndef DETMODEL_H_
#define DETMODEL_H_


#include <array>
#include <functional>
#include <utility>
#include <memory>
#include <vector>
#include <tuple>
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wconversion"
#include <armadillo>
#pragma GCC diagnostic warning "-Weffc++"
#pragma GCC diagnostic warning "-Wconversion"



//base class for a model to be simulated by determinantal quantum Monte Carlo
//
//GreenComponents is the number of independent sectors of the Green's function,
//e.g. in the Hubbard model it is 2 for spin up and spin down


typedef arma::Col<num> VecNum;
typedef arma::Mat<num> MatNum;
typedef arma::Cube<num> CubeNum;
typedef arma::Mat<int> MatInt;
typedef arma::SpMat<num> SpMatNum;


template<unsigned GreenComponents>
class DetModel {
public:
	virtual ~DetModel()
	{ }

	virtual unsigned getSystemN() const;

	//Create a MetadataMap describing the parameters of the
	//simulated model
	virtual MetadataMap prepareModelMetadataMap();

	//perform measurements of all observables
    virtual void measure();

    //get values of observables normalized by system size, the structures returned
    //contain references to the current values measured by DetHubbard.
	virtual std::vector<ScalarObservable> getScalarObservables();
	virtual std::vector<VectorObservable> getVectorObservables();
	virtual std::vector<KeyValueObservable> getKeyValueObservables();

    //perform a sweep updating the auxiliary field with costly recomputations
    //of Green functions from scratch
    virtual void sweepSimple();

    //perform a sweep as suggested in the text by Assaad with stable computation
    //of Green functions, alternate between sweeping up and down in imaginary time.
    //Will give equal-time and time-displaced Green functions.
    virtual void sweep();

protected:
    //equal-imaginary-time and time-displaced Green's functions
	//slices indexed k=0..m correspond to time slices at dtau*k,
	//which are then indexed by sites in row and column.
	//Most code, however, only uses timeslices k >= 1 ! Don't rely on g*.slice(0)
	//being valid.
	//The Green functions for k=0 are conceptually equal to those for k=m.
    std::array<CubeNum, GreenComponents> green;
    std::array<CubeNum, GreenComponents> greenFwd;
    std::array<CubeNum, GreenComponents> greenBwd;

    //The UdV-instances in UdVStorage will not move around much after setup, so storing
	//the (rather big) objects in the vector is fine
    typedef UdV<MatNum, VecNum> UdVnum;
	std::array<std::vector<UdVnum>, GreenComponents> UdVStorage;

    //functions to compute B-matrices for the different Green function sectors
    typedef std::function<MatNum(unsigned k2, unsigned k1)> FuncComputeBmat;
    std::array<FuncComputeBmat, GreenComponents> computeBmat;

    enum class SweepDirection: int {Up = 1, Down = -1};
	SweepDirection lastSweepDir;

	//observable handling -- these contain information about observables (such as their names)
	//as well as reference to their current value, which will be shared with simulation management
	//in a different class. The values reference there are to be updated here in the replica class.
	std::vector<ScalarObservable> obsScalar;
	std::vector<VectorObservable> obsVector;
	std::vector<KeyValueObservable> obsKeyValue;
};

#endif /* DETMODEL_H_ */
