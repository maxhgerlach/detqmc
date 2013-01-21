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
	bool vec;
public:
	WrongObsIndex(int obsIndex, bool vector=false) : oi(obsIndex), vec(vector) {}
    virtual const char* what() const throw () {
        return ((vec ? "Vector " : "") + std::string("Observable index ") + numToString(oi) +
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


class ConfigurationError : public std::exception {
	std::string message;
public:
	ConfigurationError(const std::string& msg) : message(msg) {}
	virtual ~ConfigurationError() throw () { }
	virtual const char* what() const throw () {
		return message.c_str();
	}
};


#endif /* EXCEPTIONS_H_ */
