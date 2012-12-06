/*
 * exceptions.h
 *
 *  Created on: Dec 6, 2012
 *      Author: gerlach
 */

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

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


#endif /* EXCEPTIONS_H_ */
