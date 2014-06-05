#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"
#include "boost/algorithm/string/join.hpp"
#pragma GCC diagnostic pop

#include "tools.h"
#include "detqmcptparams.h"
#include "exceptions.h"

void DetQMCPTParams::check() {
    //check parameters
    using namespace boost::assign;
    std::vector<std::string> neededPars;
    neededPars += "exchangeInterval", "controlParameterValues", "controlParameterName";
    for (auto p = neededPars.cbegin(); p != neededPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

    if (controlParameterValues.size() == 0) {
        throw ParameterWrong("No PT control parameters specified");
    }
}

MetadataMap DetQMCPTParams::prepareMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(VAR)
    META_INSERT(exchangeInterval);
    META_INSERT(controlParameterName);
#undef META_INSERT
    std::string valuesString = boost::algorithm::join(controlParameterValues, " ");    
    meta["controlParameterValues"] = valuesString;
    return meta;
}
