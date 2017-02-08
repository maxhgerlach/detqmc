/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0.  See the enclosed file LICENSE for a copy or if
 * that was not distributed with this file, You can obtain one at
 * http://mozilla.org/MPL/2.0/.
 *
 * Copyright 2017 Max H. Gerlach
 * 
 * */

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

    bool checkAndLogDetRatio;   // verify correctness of local spin update transition probabilities
    std::string logDetRatio_filename;

    bool checkAndLogGreen;      // verify correctness of updated Green's functions in local spin update
    std::string logGreen_filename;

    bool logGreenConsistency;   // check green consistency between wrapping / advancing, log differences

    bool checkCheckerboardConsistency; // compare results obtained by
                                       // multiplying dense / sparse
                                       // hopping matrix exponentials.
                                       // this one just prints to stdout

    std::set<std::string> specified;

    DetModelLoggingParams() :
        logSV(false), logSV_filename(""), logSV_max_filename(""), logSV_min_filename(""),
        checkAndLogDetRatio(false), logDetRatio_filename(""),
        checkAndLogGreen(false), logGreen_filename(""),
        logGreenConsistency(false),
        checkCheckerboardConsistency(false),
        specified()
    { }
    void check();
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const uint32_t version) {
        (void)version;
        ar  & logSV & logSV_filename & logSV_max_filename & logSV_min_filename
            & checkAndLogDetRatio & logDetRatio_filename
            & checkAndLogGreen & logGreen_filename
            & logGreenConsistency
            & checkCheckerboardConsistency
            & specified;
    }
};



#endif //DETMODELLOGGINGPARAMS_H
