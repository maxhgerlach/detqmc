#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"
#pragma GCC diagnostic pop

#include "tools.h"
#include "detqmcptparams.h"
#include "exceptions.h"

void DetQMCParams::check() {
    //check parameters
    using namespace boost::assign;
    std::vector<std::string> neededPars;
    neededPars += "exchangeInterval", "controlParameterValues";
    for (auto p = neededPars.cbegin(); p != neededPars.cend(); ++p) {
        if (specified.count(*p) == 0) {
            throw ParameterMissing(*p);
        }
    }

    if (controlParameterValues.size() == 0) {
        throw ParameterWrong("No PT control parameters specified");
    }
}

MetadataMap DetQMCParams::prepareMetadataMap() const {
    MetadataMap meta;
#define META_INSERT(VAR) meta[#VAR] = numToString(VAR)
    META_INSERT(exchangeInterval);
#undef META_INSERT
    //controlParameterValues
    return meta;
}
