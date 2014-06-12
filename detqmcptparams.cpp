#include <functional>           // std::function
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wshadow"
#include "boost/assign/std/vector.hpp"
#include "boost/algorithm/string/join.hpp"
#include "boost/range/adaptor/transformed.hpp"
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
    using boost::algorithm::join;
    using boost::adaptors::transformed;
    // std::string valuesString = join(
    //     controlParameterValues | transformed(numToString<num>),
    //     " ");
    
    /// The argument to transformed needs to fulfill some qualities,
    /// e.g. it must be copy-constructible, the simple statement above
    /// does not compile on the Intel Compiler, C++11 lambdas do not
    /// work either (!
    /// http://stackoverflow.com/questions/11872558/using-boost-adaptors-with-c11-lambdas),
    /// apparently they need to be a functor with member result_type.
    /// Work around: std::function
    /// -- This may or may not be better in a more recent version of
    /// Boost than 1.51 --
    std::function<std::string(num)> func = [](num v) {
        return numToString<num>(v);
    };
    std::string valuesString = join(
        controlParameterValues | transformed(func),
        " ");
    meta["controlParameterValues"] = valuesString;
    return meta;
}
