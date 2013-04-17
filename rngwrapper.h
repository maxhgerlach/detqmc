//Code by Max Henner Gerlach, 2010--2012,
//used for the diploma thesis "Directional Ordering in the Classical Compass Model in Two and Three Dimensions"
//contact: maxgerlach@gmail.com

/*
 * rngwrapper.h
 *
 *  Created on: Oct 12, 2011
 *      Author: gerlach
 *
 */

// C++ Wrapper around DSFMT random number generator

#ifndef RNGWRAPPER_H_
#define RNGWRAPPER_H_

#include <string>
extern "C" {
#include "dsfmt/dSFMT.h"
}

#include <boost/serialization/string.hpp>

class RngWrapper {
    unsigned long seed;
    unsigned int processIndex;
    dsfmt_t dsfmt;
public:
    RngWrapper(unsigned long seed_ = 0, unsigned processIndex_ = 0);
    virtual ~RngWrapper() {}

    std::string getName();

    //return a floating point random number from (0, 1)
    double rand01() {
        //dSFTM
        return dsfmt_genrand_open_open(&dsfmt);
    }

    //return a floating point random number from (low, high)
    double randRange(double low, double high) {
    	return low + (high - low) * rand01();
    }

    //return an integer random number from the discrete uniform distribution over the integers
    //{low, low + 1, . . . , high}. One call to rand01()
    int randInt(int low, int high) {
        return low + static_cast<int> ((high - low + 1.0) * rand01());
    };

    //saving/loading state without Boost serialization
    void saveState();
    void loadState();

private:
    //for serialization with Boost
	friend class boost::serialization::access;

    std::string stateToString();
    void stringToState(const std::string& stateString);

    template<class Archive>
	void save(Archive& ar, const unsigned int version) {
    	ar << seed << processIndex;
    	ar << stateToString();
    }

    template<class Archive>
	void load(Archive& ar, const unsigned int version) {
    	ar >> seed >> processIndex;
    	std::string stateString;
    	ar >> stateString;
    	stringToState(stateString);
    }

};

#endif /* RNGWRAPPER_H_ */
