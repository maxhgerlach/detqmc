#!/usr/bin/env python
import os
import sys
import os.path
import itertools
import subprocess
from glob import glob
import argparse
from scripthelpers import getHeader, parseHeader, is_file_older_than_a_day


evalOnlyNew = False
evalOnlyRecent = False



# parse command line arguments and start program
################################################
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", type=str, default="",
                    help="consider directories with this prefix for tauint evaluation")
parser.add_argument("variable", type=str,
                    help="quantity to use as primary independent variable")
parser.add_argument("-i",  nargs='*', type=str,
                    default=[],
                    help="list of additional variables not to be discerned if"+
                    "they occur with multiple values")
parser.add_argument("-s",  nargs='*', type=str,
                    help="list of variables for which multiple values should be discerned "+
                    "even though they normally are not; these are always added to output filenames too.")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
parser.add_argument("-n", nargs='*', type=str,
                    help="list of observables not to collect or evaluate (will also be added to @options)")
parser.add_argument("--options", type=str, default="--noexp",
                    help="options to be passed to deteval, default: --noexp")
parser.add_argument("--simindexjoined", dest='simindexjoined', action='store_true', help="only consider subdirectories with data joined from multiple simulations (by joinall.py, jointimeseries).  If this is not given, data with simindex=joined is ignored.  If this is given, simindex is added to the list otherwise given with option -s.")
parser.set_defaults(simindexjoined=False)
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)
parser.add_argument("--evalOnlyRecent", dest='evalOnlyRecent', action='store_true', help="skip already evaluated subdirectories, if data is older than 24 hours")
parser.add_argument("--no-evalOnlyRecent", dest='evalOnlyRecent', action='store_false', help="do not skip already evaluated subdirectories, if data is older than 24 hours")
parser.set_defaults(evalOnlyRecent=False)

args = parser.parse_args()
variable = args.variable
options = args.options
prefix = args.prefix
nonMultiVal = [ # compilation related keys:
    "ARMA_VERSION", "BOOST_LIB_VERSION", "buildDate", "buildHost", "buildTime", "cppflags", "cxxflags", 
    "gitBranch", "gitRevisionHash",
    # etc.:
    "averageAcceptedWolffClusterSize",
    "globalShiftAccRatio",
    "wolffClusterUpdateAccRatio",
    "wolffClusterShiftUpdateAccRatio",
    "controlParameterValues",
    "controlParameterCount",
    "eval-jackknife-blocks",
    "sweeps", "thermalization",
    "turnoffFermionMeasurements",
    "join-total-samples",
    "join-subsample",
    "join-discard",
    "join-read"] \
    + ["m", "dtau", "rngSeed", "N", "alpha", "eval-samples"]
extraNonMultiVal = args.i if args.i else []
nonMultiVal.extend(extraNonMultiVal)
extraMultiVal = args.s if args.s else []
excludedSubdirs = args.x if args.x else []
noncollectObs = args.n if args.n else []
if noncollectObs:
    options += " -n " + " ".join(noncollectObs)    

only_simindexjoined = args.simindexjoined
if only_simindexjoined:
    extraMultiVal.append("simindex")
evalOnlyNew = args.evalOnlyNew
evalOnlyRecent = args.evalOnlyRecent


if set(extraMultiVal) & set(extraNonMultiVal):
    # intersection
    raise Exception("Arguments to options -i and -s are not disjoint")


# also include directories further down in the hierarchy,
# exclude everything below a tree FAILED or belonging to the list of subdirectories
# specified with '-x'
# http://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
subdir_candidates = []
for root, dirs, files in os.walk('.', topdown=True, followlinks=True):
    dirs[:] = [d for d in dirs if d not in excludedSubdirs and d not in ['FAILED', 'wrong-r-range', 'histograms']
               and not d.startswith('eval-discard-')]
    if root != '.':
        subdir_candidates.append(root)
print "collected subdir_candidates"
sys.stdout.flush()
   
subdirs = [f for f in subdir_candidates if (f.startswith(prefix) or
                                            f.startswith("./" + prefix))
           and (variable in f or
                (variable in ["m", "beta", "temp"]
                 and ("m" in f or "beta" in f or "temp" in f)))]
print "pruned subdir_candidates for prefix"
sys.stdout.flush()


# potentially append "_" for the output prefix
if prefix != "":
    if prefix[-1:] == "-":
        prefix = prefix[:-1] + "_" # replace "-" in tail by "_"
    elif prefix[-1:] != "_":
        # append a "_" if it is not there already
        prefix = prefix + "_"



#collect info.dat contents (only common entries), potentially prune non simindexjoined
######################################################################################
infodata = {}
for sd in subdirs:
    header = parseHeader(getHeader(sd + "/info.dat"))
    if not header is None:
        if only_simindexjoined:
            if not "simindex" in header or header["simindex"] != "joined":
                # ignore individual simindex data
                continue
            infodata[sd] = header                
        else:
            if "simindex" in header and header["simindex"] == "joined":
                # only take into account individual simindex data
                continue
            infodata[sd] = header
subdirs = infodata.keys()
# print subdirs
        
print "collected subdir info.dat contents"
sys.stdout.flush()
infodata = {k:v for k,v in infodata.iteritems() if v is not None}
commoninfodata = dict(set.intersection(*(set(d.iteritems()) for d in infodata.itervalues())))
commoninfodata_filename = prefix + "info-common.dat"
outf = open(commoninfodata_filename, "w")
for k,v in commoninfodata.items():
    outf.write('# %s = %s\n' % (k, v))
outf.close()
print "wrote common info.dat contents to file", commoninfodata_filename
sys.stdout.flush()


#evaluate subdirectories
########################
maindirectory = os.getcwd()
def eval_subdir(sd):
    "if possible call deteval for sd, otherwise return it as skipped"
    print sd,
    if not glob(sd + "/*.series"):
        print "=> no time series found: skip"
        sys.stdout.flush()
        return sd
    elif evalOnlyRecent and os.path.exists(sd + "/eval-tauint.values") and is_file_older_than_a_day(sd + "/info.dat"):
        print "=> already evaluated and data is older than 24 hours: skip evaluation, but take into account"
        sys.stdout.flush()
        return ''
    elif evalOnlyNew and os.path.exists(sd + "/eval-tauint.values"):
        print "=> already evaluated: skip evaluation, but take into account"
        sys.stdout.flush()
        return ''
    else:
        os.chdir(sd)
        print "=> evaluate"
        sys.stdout.flush()
        stdout_and_stderr = subprocess.check_output("deteval " + options + " ; exit 0", shell=True,
                                                    stderr=subprocess.STDOUT)
        print stdout_and_stderr
        if not os.path.exists("eval-tauint.values"):
            print "=> did not write eval-tauint.values, skip after all"
            skip_if_sd = sd
        else:
            skip_if_sd = ''
        os.chdir(maindirectory)
        return skip_if_sd

import multiprocessing as mp
pool_size = os.getenv("SLURM_CPUS_PER_TASK")
if pool_size is None:
    pool_size = max(4, mp.cpu_count())
else:
    pool_size = int(pool_size)
pool = mp.Pool(processes=pool_size)
skipped_subdirs = pool.map(eval_subdir, subdirs)
for sd in skipped_subdirs:
    if sd != '':
        subdirs.remove(sd)


#collect eval-results
#####################

# helper:
def addControlParameterCount(meta_dict):
    if meta_dict is None:
        return None
    if "controlParameterValues" in meta_dict:
        meta_dict["controlParameterCount"] = str(len(meta_dict["controlParameterValues"].split()))
    return meta_dict

# map: subdirectory -> metadata dictionary [for replica exchange simulations: count controlParameterValues]
metadata = {sd: addControlParameterCount(parseHeader(getHeader(sd + "/eval-tauint.values"))) for sd in subdirs}
# prune subdirectories with empty metadata
metadata = {sd:meta for sd,meta in metadata.iteritems() if meta is not None}
# go over all the metadata dictionaries, each time take the keys of those dictionaries, then find all 
# the common ones (set intersection)
commonkeys = set.intersection(*(set(d.iterkeys()) for d in metadata.itervalues())) # * : reverse of zip
# map commonkeys -> metadata ; only if metadata also equal
commonmetadata = dict(set.intersection(*(set(d.iteritems()) for d in metadata.itervalues())))

# mykeys: all the metadata-keys that have set values atleast in one subdirectory.
# this means that for some subdirectories that metadata may not be set. In this case use
# a default value "NV" ("NOVALUE")
mykeys = set.union(*(set(d.iterkeys()) for d in metadata.itervalues())) # * : reverse of zip
for sd, meta in metadata.items():
    for k in mykeys:
        if not k in meta.keys():
            meta[k] = "NV"
    
#print commonmetadata
mapKeySubdirValue = { k : {sd : metadata[sd][k] for sd in subdirs if sd in metadata.keys()} for k in mykeys }
multivalKeys = [k for k in mykeys if len(set(mapKeySubdirValue[k].values())) > 1 or k in extraMultiVal]

#don't consider some keys:
multivalKeys = [k for k in multivalKeys if k in extraMultiVal or k not in nonMultiVal]
#also don't leave our variable key in that list
multivalKeys = [k for k in multivalKeys if k != variable]
#special handling for temp/beta/m
temperature_variables = ["temp", "beta", "m"]
if variable in temperature_variables:
    for to_be_removed in [v for v in temperature_variables if v != variable and v in multivalKeys]:
        multivalKeys.remove(to_be_removed)

print multivalKeys
    
multivalValueLists = [sorted(set(mapKeySubdirValue[k].values())) for k in multivalKeys]  # sorted(set(..)) to uniquify
#print multivalValueLists
multivalValueTuples = [t for t in itertools.product(*multivalValueLists)]
#print multivalValueTuples

#print len(multivalValueTuples)

#mapKeyValueSubdirs = { k : {v: [s for s in subdirs if mapKeySubdirValue[k][s] == v]} for k in commonkeys }

for tup in multivalValueTuples:
    if len(multivalKeys) == 0:
        relsubdirs = [s for s in subdirs if metadata.has_key(s)]              # take all with metadata
    else:
        def allKeysMatched(subdir):
            matched = True
            for k,v in zip(multivalKeys, tup):
                if metadata[subdir][k] != v:
                    matched = False
                    break
            return matched
        relsubdirs = filter(allKeysMatched, subdirs)

    print  "_".join(["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)])
    #now relsubdirs only differ in the values of the chosen variable
    #
    #map: observableName -> variableValue -> tauint-estimate
    mapToTauint = {}
    for sd in relsubdirs:
        print " ", sd,
        for k in multivalKeys:
            print k + str(": ") + metadata[sd][k],
        print
        
        #if the job gave any error messages, print them:
        errorfile = "%s/errors.log" % sd
        if os.path.exists(errorfile) and os.path.getsize(errorfile) > 0:
            with open(errorfile, 'r') as f:
                print sd, " errors:"
                for line in f:
                    print line
                print
        #parse evaluated tauint file:
        if variable == 'temp':
            varvalue = 1.0 / float(metadata[sd]['beta'])
        else:
            varvalue = float(metadata[sd][variable])
        tauintFile = open(sd + "/eval-tauint.values", 'r')
        rows = (row.strip().split() for row in tauintFile if row[0] != '#')
        for row in rows:
            if row:
                obs = row[0]
                tau = row[1]
                mapToTauint.setdefault(obs, {})[varvalue] = tau
    for obsname in mapToTauint.keys():
        if obsname in noncollectObs:
            continue

        outfname = prefix + "eval-tauint-" + variable + "-" + obsname + "_" + "_".join(
            ["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)]) + ".values"
        outf = open(outfname, 'w')
        outmeta = commonmetadata.copy()
        outmeta["key"] = variable
        outmeta["observable"] = obsname
        for k,v in zip(multivalKeys, tup):
            outmeta[k] = v
        for k,v in outmeta.items():
            outf.write('# %s = %s\n' % (k, v))
        outf.write('## %s,\t %s\n' % (variable, obsname))
        for varvalue in sorted(mapToTauint[obsname].keys()):
            tau = mapToTauint[obsname][varvalue]
            if tau != "nan":
                outf.write('%s\t%s' % (varvalue, tau))
            else:
                outf.write('# %s -- no value' % varvalue)
            outf.write('\n')
        outf.close()
