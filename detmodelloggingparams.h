#ifndef DETMODELLOGGINGPARAMS_H
#define DETMODELLOGGINGPARAMS_H

#include <vector>
#include <string>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"
#include "boost/assign/std/vector.hpp"    // 'operator+=()' for vectors
#pragma GCC diagnostic pop            


#include "metadata.h"
#include "exceptions.h"



// struct holding parameters concering the logging of specific data in DetModelGC derived classes

struct DetModelLoggingParams {
    bool logSV;                 // log Green's function singular value range
    std::string logSV_filename; // the filename associated to that log  -- we set this in createReplica()

    std::set<std::string> specified;


    DetModelLoggingParams() :
        logSV(false), logSV_filename(""),
        specified()
    { }
    void check();
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar  & logSV & logSV_filename
            & specified;
    }
};



#endif //DETMODELLOGGINGPARAMS_H
