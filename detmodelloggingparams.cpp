#include "detmodelloggingparams.h"


void DetModelLoggingParams::check() {
    if (specified.count("logSV") and not specified.count("logSV_filename")) {
        throw ParameterMissing("logSV_filename");
    }
}
