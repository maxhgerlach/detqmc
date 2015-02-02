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

class GeneralError : public std::exception {
    std::string message;
public:
    GeneralError(const std::string& msg, const char* file = 0, int line = -1)
    {
        if (file and line >= 0) {
            message = file + (":" + numToString(line)) + " -> " + msg;
        } else {
            message = msg;
        }
    }
    virtual ~GeneralError() throw () { }
    virtual const char* what() const throw () {
        return message.c_str();
    }
};

class WrongObsIndex : public GeneralError {
public:
    WrongObsIndex(int obsIndex, bool vector=false, const char* file = 0, int line = -1)
        : GeneralError(((vector ? "Vector " : "") + std::string("Observable index ") + numToString(obsIndex) +
                        " not supported"), file, line)
        { }
};


class ParameterMissing : public GeneralError {
public:
    ParameterMissing(const std::string& par, const char* file = 0, int line = -1)
        : GeneralError("Parameter " + par + " not given.", file, line)
        { }
};


class ParameterWrong : public GeneralError {
public:
    template <typename ValType>
    ParameterWrong(const std::string& par, const ValType& val, const char* file = 0, int line = -1)
        : GeneralError("Parameter " + par + " has incorrect value "
                       + numToString(val), file, line)
        { }
    ParameterWrong(const std::string& message, const char* file = 0, int line = -1)
        : GeneralError(message, file, line)
        { }
};


class ConfigurationError : public GeneralError {
public:
    ConfigurationError(const std::string& msg, const char* file = 0, int line = -1)
        : GeneralError(msg, file, line)
        { }
};


class SerializationError : public GeneralError {
public:
    SerializationError(const std::string& msg, const char* file = 0, int line = -1)
        : GeneralError(msg, file, line)
        {}
};



class ReadError : public GeneralError {
public:
    ReadError(std::string filename, const char* file = 0, int line = -1) throw ()
        : GeneralError("Can't read from file " + filename, file, line)
        { }
};

class KeyUndefined : public GeneralError {
public:
    KeyUndefined(const std::string& key = "", const char* file = 0, int line = -1)
        : GeneralError("key undefined: " + key, file, line)
        { }
};


#endif /* EXCEPTIONS_H_ */
