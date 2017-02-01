#!/usr/bin/env python

import re
import sys
import numpy as np
from scripthelpers import getHeader, parseHeader

# meta files for the replica exchange process contain one column with
# the control parameter index ("cpi") and one with the relevant data.
# this script exchanges the cpi with the actual control parameter values



if __name__ == "__main__":
    assert len(sys.argv) == 3
    filename = sys.argv[1]
    outputfilename = sys.argv[2]

    header = getHeader(filename)
    
    data = np.loadtxt(filename)
    cpi = data[:,0]

    cpv = np.array([float(s) for s in parseHeader(header)["controlParameterValues"].split(" ")])
    assert len(cpv) == len(cpi)

    data[:,0] = cpv

    np.savetxt(outputfilename, data, header="".join(header).replace("control parameter index", "control parameter value"))
