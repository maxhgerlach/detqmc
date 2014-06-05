#ifndef DETQMCPTPARAMS_H
#define DETQMCPTPARAMS_H

#include <vector>
#include <string>
#include <set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/serialization/string.hpp"
#include "boost/serialization/set.hpp"
#include "boost/serialization/vector.hpp"
#pragma GCC diagnostic pop

#include "metadata.h"

// Collect various structs defining various parameters. [more structs
// in the other ...params.h headers]

// The set specified included in each struct contains string representations
// of all parameters actually specified.  This allows throwing an exception
// at the appropriate point in the program if a parameter is missing.


// Representation of replica exchange / parallel tempering Monte Carlo
// simulation parameters for DetQMC
struct DetQMCPTParams {
    uint32_t exchangeInterval;
    std::vector<num> controlParameterValues;
    std::string controlParameterName;
    
    std::set<std::string> specified; // used to record names of specified parameters

    DetQMCPTParams() :
        exchangeInterval(0), controlParameterValues(), controlParameterName(""), specified()
    { }

    // check consistency
    void check();

    // to export human readable form of parameters
    MetadataMap prepareMetadataMap() const;
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar & exchangeInterval & controlParameterValues & controlParameterName
           & specified;
    }
};

#endif /* DETQMCPTPARAMS_H */
