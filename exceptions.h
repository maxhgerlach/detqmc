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


class ParameterWrong : public std::exception {
	std::string message;
public:
	template <typename ValType>
	ParameterWrong(const std::string& par, const ValType& val) :
			message("Parameter " + par + " has incorrect value "
					+ numToString(val)) {
	}
	ParameterWrong(const std::string& msg) : message(msg) {
	}
	virtual ~ParameterWrong() throw () { }
	virtual const char* what() const throw () {
		return message.c_str();
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


class SerializationError : public std::exception {
	std::string message;
public:
	SerializationError(const std::string& msg) : message(msg) {}
	virtual ~SerializationError() throw () { }
	virtual const char* what() const throw () {
		return message.c_str();
	}
};



class ReadError : public std::exception {
    std::string filename;
public:
    ReadError(std::string filename) throw ():
        filename(filename) {
    }
    ~ReadError() throw () { }
    virtual const char* what() const throw () {
        return ("Can't read from file " + filename).c_str();
    }
};

class KeyUndefined : public std::exception {
    std::string key;
public:
    KeyUndefined() : key("")
    { }
    KeyUndefined(std::string k) : key(k)
    { }

    virtual const char* what() const throw () {
        return ("key undefined: " + key).c_str();
    }

    ~KeyUndefined() throw() { }
};


#endif /* EXCEPTIONS_H_ */
