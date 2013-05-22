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

//share

template <typename ObsType>
struct Observable {
	typedef std::reference_wrapper<const ObsType> RefType;
	RefType valRef;
	std::string name;
	std::string shortName;

	Observable(RefType v, const std::string& name_, const std::string& short_)
		: valRef(v), name(name_), shortName(short_)
	{ }
	virtual ~Observable() { }
};


typedef Observable<num> ScalarObservable;


struct VectorObservable : public Observable<arma::Col<num>> {
	uint32_t vectorSize;

	VectorObservable(RefType v, uint32_t vectorSize_,
			const std::string& name_, const std::string& short_)
		: Observable<arma::Col<num>>(v, name_, short_),
		  vectorSize(vectorSize_)
	{ }
};


struct KeyValueObservable : public VectorObservable {
	arma::Col<num> keys;
	std::string keyName;

	KeyValueObservable(RefType v, const arma::Col<num>& keys_,
			const std::string& keyName_,
			const std::string& name_, const std::string& short_)
		: VectorObservable(v, keys_.n_elem, name_, short_), keys(keys_),
		  keyName(keyName_)
	{ }
};

#endif /* OBSERVABLE_H_ */
