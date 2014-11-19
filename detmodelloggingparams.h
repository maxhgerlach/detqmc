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
    bool logSV;                 // log Green's function singular value range, max, min
    std::string logSV_filename; // the filename associated to that range log  -- we set this in createReplica()
    std::string logSV_max_filename; // for the max sv
    std::string logSV_min_filename; // for the min sv
    
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
