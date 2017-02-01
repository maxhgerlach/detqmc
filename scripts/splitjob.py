#!/usr/bin/env python
from collections import OrderedDict
import itertools
import os
import os.path
import errno
import re
import sys
import subprocess
from subprocess import Popen, PIPE, check_output
import argparse



# take a string with configuration settings confOption = confValue1, confValue2, ... \n
# return a ordered dictionary mapping confOption -> [List of confValues]
# If the list ends in a "," not followed by another confValue, this confOption 
# will be used to construct the jobname regardless of it having only one possible
# value or not.
def parseConf(header):
    d = OrderedDict()
    # Regexp to match "key = anything":
    p = re.compile(r"\s* ([^\s=]+) \s* = (.+)", re.VERBOSE)
    for line in header:
        if line[0] == "#":
            continue
        m = p.match(line)
        try:
            d[m.group(1)] = [s.strip() for s in m.group(2).split(',')]
        except:
            pass
    return d
 

def writeConf(jobconfdict, jobfilename):
    with open(jobfilename, 'w') as jobfile:
        for key in jobconfdict:
            jobline = key + " = " + ", ".join(jobconfdict[key]) + '\n'
            jobfile.write(jobline)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # only positional arguments
    parser.add_argument("jobtemplate", type=str, help=".job filename of template job that should be split up")
    parser.add_argument("splitvar", type=str, help="name of job option, defined with multiple values in the job template file, from which multiple split jobs should be created")
    args = parser.parse_args()
    splitvar = args.splitvar
    jobtemplate = args.jobtemplate

    if not jobtemplate.endswith(".job"):
        raise Exception("jobtemplate filename needs to end with `.job`")

    with open(jobtemplate, 'r') as templatefile:
        templateconf = parseConf(templatefile)

    splitvalues = templateconf[splitvar]

    for sv in splitvalues:
        if sv != '':
            jobconf = templateconf.copy()
            # the following will lead to a trailing comma, to be included in the job name
            jobconf[splitvar] = [sv, ""]
            
            basename, dotjob, _ = jobtemplate.rpartition(".job")
            jobfilename = basename + "_" + splitvar + sv + dotjob

            writeConf(jobconf, jobfilename)


        
