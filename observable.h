/*
 * observable.h
 *
 *  Created on: Feb 4, 2013
 *      Author: gerlach
 */

#ifndef OBSERVABLE_H_
#define OBSERVABLE_H_

#include <string>
#include <functional>
#include <armadillo>
#include "parameters.h"

//a simple handle holding an observable value reference to be shared by a replica class
//implementing a model and a simulation class handling the calculation of expectation values
//and so on

template <typename ObsType>
struct Observable {
	std::reference_wrapper<const ObsType> valRef;		//this will be a reference to an object held by the replica!
	std::string name;
	std::string shortName;

	Observable(const ObsType& valRef_, const std::string& name_, const std::string& short_)
		: valRef(std::cref(valRef_)), name(name_), shortName(short_)
	{ }
	virtual ~Observable() { }
};


typedef Observable<num> ScalarObservable;


struct VectorObservable : public Observable<arma::Col<num>> {
	unsigned vectorSize;

	VectorObservable(const arma::Col<num>& valRef_, unsigned vectorSize_,
			const std::string& name_, const std::string& short_)
		: Observable<arma::Col<num>>(valRef_, name_, short_),
		  vectorSize(vectorSize_)
	{ }
};


struct KeyValueObservable : public VectorObservable {
	arma::Col<num> keys;
	std::string keyName;

	KeyValueObservable(const arma::Col<num>& valRef_, const arma::Col<num>& keys_,
			const std::string& keyName_,
			const std::string& name_, const std::string& short_)
		: VectorObservable(valRef_, keys.n_rows, name_, short_), keys(keys_),
		  keyName(keyName_)
	{ }
};

#endif /* OBSERVABLE_H_ */
