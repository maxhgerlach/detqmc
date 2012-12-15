/*
 * exceptions.h
 *
 *  Created on: Dec 6, 2012
 *      Author: gerlach
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <exception>
#include <string>
#include "tools.h"

class WrongObsIndex : public std::exception {
	int oi;
public:
	WrongObsIndex(int obsIndex) : oi(obsIndex) {}
    virtual const char* what() const throw () {
        return ("Observable index " + numToString(oi) +
        		" not supported").c_str();
    }
};


class ParameterMissing : public std::exception {
	std::string par;
public:
	ParameterMissing(const std::string& par) : par(par) {}
	virtual ~ParameterMissing() throw () { }
	virtual const char* what() const throw () {
		return ("Parameter " + par + " not given.").c_str();
	}
};


template <typename ValType>
class ParameterWrong : public std::exception {
	std::string par;
	ValType val;
public:
	ParameterWrong(const std::string& par,
			const ValType& val) : par(par), val(val) {}
	virtual const char* what() const throw () {
		return ("Parameter " + par + " has incorrect value "
				+ numToString(val));
	}
};


#endif /* EXCEPTIONS_H_ */
