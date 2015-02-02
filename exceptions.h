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
#define throw_GeneralError(msg) throw GeneralError(msg, __FILE__, __LINE__)

class WrongObsIndex : public GeneralError {
public:
    WrongObsIndex(int obsIndex, bool vector=false, const char* file = 0, int line = -1)
        : GeneralError(((vector ? "Vector " : "") + std::string("Observable index ") + numToString(obsIndex) +
                        " not supported"), file, line)
        { }
};
#define throw_WrongObsIndex(obsIndex, vector) throw WrongObsIndex(obsIndex, vector, __FILE__, __LINE__)


class ParameterMissing : public GeneralError {
public:
    ParameterMissing(const std::string& par, const char* file = 0, int line = -1)
        : GeneralError("Parameter " + par + " not given.", file, line)
        { }
};
#define throw_ParameterMissing(par) throw ParameterMissing(par, __FILE__, __LINE__)


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
#define throw_ParameterWrong_message(message) throw ParameterWrong(message, __FILE__, __LINE__)
#define throw_ParameterWrong(par, val) throw ParameterWrong(par, val, __FILE__, __LINE__)


class ConfigurationError : public GeneralError {
public:
    ConfigurationError(const std::string& msg, const char* file = 0, int line = -1)
        : GeneralError(msg, file, line)
        { }
};
#define throw_ConfigurationError(message) throw ConfigurationError(message, __FILE__, __LINE__)


class SerializationError : public GeneralError {
public:
    SerializationError(const std::string& msg, const char* file = 0, int line = -1)
        : GeneralError(msg, file, line)
        {}
};
#define throw_SerializationError(message) throw SerializationError(message, __FILE__, __LINE__)



class ReadError : public GeneralError {
public:
    ReadError(std::string filename, const char* file = 0, int line = -1) throw ()
        : GeneralError("Can't read from file " + filename, file, line)
        { }
};
#define throw_ReadError(message) throw ReadError(message, __FILE__, __LINE__)


class KeyUndefined : public GeneralError {
public:
    KeyUndefined(const std::string& key = "", const char* file = 0, int line = -1)
        : GeneralError("key undefined: " + key, file, line)
        { }
};
#define throw_KeyUndefined(key) throw KeyUndefined(key, __FILE__, __LINE__)


#endif /* EXCEPTIONS_H_ */
