#!/usr/bin/env python
import os
import os.path
import operator                 # itemgetter
import itertools
import subprocess
from glob import glob
import argparse
from scripthelpers import mkdir_p, getHeader, parseHeader

evalOnlyNew = False
#evalOnlyNew = True


# parse command line arguments and start program
################################################
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", type=str, default="",
                    help="consider directories with this prefix for evaluation")
parser.add_argument("variable", type=str,
                    help="quantity to use as primary independent variable")
parser.add_argument("-i",  nargs='*', type=str,
                    default=[],
                    help="list of additional variables not to be discerned if"+
                    "they occur with multiple values")
parser.add_argument("-s",  nargs='*', type=str,
                    help="list of variables for which multiple values should be discerned "+
                    "even though they normally are not")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
parser.add_argument("--options", type=str, default="-j 20",
                    help="options to be passed to detevalbc, default: -j 20")
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)


args = parser.parse_args()
evalOnlyNew = args.evalOnlyNew
variable = args.variable
options = args.options
prefix = args.prefix
nonMultiVal = ["totalWallTimeSecs",
               "averageAcceptedWolffClusterSize",
               "globalShiftAccRatio",
               "wolffClusterUpdateAccRatio",
               "wolffClusterShiftUpdateAccRatio",
               "controlParameterValues",
               "controlParameterCount",
               "sweeps", "thermalization",
               "sweepsDone", "sweepsDoneThermalization"] \
    + ["m", "dtau", "rngSeed", "N", "alpha", "eval-samples",
       "eval-samples_pbc", "eval-samples_apbc-x", "eval-samples_apbc-y", "eval-samples_apbc-xy"]
extraNonMultiVal = args.i if args.i else []
nonMultiVal.extend(extraNonMultiVal)
extraMultiVal = args.s if args.s else []
excludedSubdirs = args.x if args.x else []

if set(extraMultiVal) & set(extraNonMultiVal):
    # intersection
    raise Exception("Arguments to options -i and -s are not disjoint")


# also include directories further down in the hierarchy,
# exclude everything below a tree FAILED or belonging to the list of subdirectories
# specified with '-x'
# http://stackoverflow.com/questions/19859840/excluding-directories-in-os-walk
subdir_candidates = []
for root, dirs, files in os.walk('.', topdown=True, followlinks=True):
    dirs[:] = [d for d in dirs if d not in excludedSubdirs and d not in ['FAILED']]
    if root != '.':
        subdir_candidates.append(root)
   
subdirs = [f for f in subdir_candidates if (f.startswith(prefix) or
                                            f.startswith("./" + prefix))
           and glob(f + "/*.series")  # exclude subdirs without time series -- e.g. pt top directories
           and (variable in f or
                (variable in ["m", "beta", "temp"]
                 and ("m" in f or "beta" in f or "temp" in f)))]

# potentially append "_" for the output prefix
if prefix != "":
    if prefix[-1:] == "-":
        prefix = prefix[:-1] + "_" # replace "-" in tail by "_"
    elif prefix[-1:] != "_":
        # append a "_" if it is not there already
        prefix = prefix + "_"



#collect info.dat contents
##########################
infodata = {}
for sd in subdirs:
    header = parseHeader(getHeader(sd + "/info.dat"))
    if not header is None:
        infodata[sd] = header
infodata = {k:v for k,v in infodata.iteritems() if v is not None}


# We need to group the subdirectories in groups that differ only by the boundary conditions
###########################################################################################
groupable_list = []
# keys in the following links are removed from the metadata before sorting by the metadata
exclude_from_metadata = ["bc", # in one group we want to have all 4 possible b.c.'s
                         # compilation related keys:
                         "ARMA_VERSION", "BOOST_LIB_VERSION", "buildDate", "buildHost", "buildTime", "cppflags", "cxxflags", "gitBranch", "gitRevisionHash",
                         # rng seed will be different:
                         "rngSeed",
                         # target sweeps could have been modified for some
                         "sweeps", "thermalization",
                         # current state of simulation will be different
                         "sweepsDone", "sweepsDoneThermalization", "totalWallTimeSecs"]
for sd, metadict in infodata.iteritems():
    bc = metadict["bc"]
    interesting_metadict = { k: metadict[k] for k in metadict if k not in exclude_from_metadata }
    groupable_list.append( (sd, bc, sorted(interesting_metadict.iteritems())) )

groupable_list.sort(key = operator.itemgetter(2))

# for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(2)):
#     output_directory = None
#     for entry in entries:
#         # print bc, sd
#         print entry[1], entry[0]
#         # generate target directory name -- boundary condition renamed to "averaged"
#         if output_directory is None:
#             output_directory = entry[0].replace(entry[1], "averaged")
#     print "->", output_directory
#     print


# import sys; sys.exit(0)


#evaluate subdirectory-groups
#############################

# turn into list of lists, one sub-list for each group

bc_subdir_groups = [ [ (entry[1], entry[0]) for entry in entries ]
                     for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(2)) ]

# maindirectory = os.getcwd()
# def eval_subdir(sd):
#     "if possible call deteval for sd, otherwise return it as skipped"
#     if evalOnlyNew and os.path.exists(sd + "/eval-results.values"):
#         print "already evaluated: skip", sd
#         return sd
#     elif not glob(sd + "/*.series"):
#         print "no time series found: skip", sd
#         return sd
#     else:
#         os.chdir(sd)
#         print sd
#         stdout_and_stderr = subprocess.check_output("deteval " + options + " ; exit 0", shell=True,
#                                                     stderr=subprocess.STDOUT)
#         print stdout_and_stderr
#         os.chdir(maindirectory)
#         return ''

def show_group(tuple_list):
    output_directory = None
    if len(tuple_list) != 4:
        print "faulty group: ", tuple_list
        return ""
    for bc, sd in tuple_list:
#        print bc, sd
        # generate target directory name -- boundary condition renamed to "averaged"
        if output_directory is None:
            output_directory = sd.replace(bc, "averaged")
    return output_directory

def eval_bc_subdirs(tuple_list):
    " tuple_list : [(bc, subdir), ...]. return output directory if detevalbc has been run "
    if len(tuple_list) != 4:
        print "faulty group: ", tuple_list
        return ""
    output_directory = tuple_list[0][1].replace("bc"+tuple_list[0][0], "bc"+"averaged")
    if evalOnlyNew and os.path.exists(output_directory + "/eval-results.values"):
        print "already evaluated: skip", output_directory
        return output_directory
    mkdir_p(output_directory)
    commandline = "detevalbc " + options + " --outputDirectory " + output_directory
    for bc, sd in tuple_list:
        commandline = commandline + " " + sd
    commandline = commandline + " ; exit 0"
    print commandline
    print output_directory
    stdout_and_stderr = subprocess.check_output(commandline, shell=True,
                                                stderr=subprocess.STDOUT)
    print stdout_and_stderr
    return output_directory

import multiprocessing as mp
pool_size = max(4, os.getenv("SLURM_CPUS_PER_TASK"))
if pool_size is None:
    pool_size = mp.cpu_count()
else:
    pool_size = int(pool_size)
pool = mp.Pool(processes=pool_size)

output_directories = pool.map(eval_bc_subdirs, bc_subdir_groups, 1)
output_directories = [od for od in output_directories if od != '']
# print output_directories


#collect eval-results from the output directories
#################################################

# helper:
def addControlParameterCount(meta_dict):
    if "controlParameterValues" in meta_dict:
        meta_dict["controlParameterCount"] = str(len(meta_dict["controlParameterValues"].split()))
    return meta_dict

# map: subdirectory -> metadata dictionary [for replica exchange simulations: count controlParameterValues]
metadata = {sd: addControlParameterCount(parseHeader(getHeader(sd + "/eval-results.values"))) for sd in output_directories}
# prune subdirectories with empty metadata
metadata = {sd:meta for sd,meta in metadata.iteritems() if meta is not None}
# go over all the metadata dictionaries, each time take the keys of those dictionaries, then find all 
# the common ones (set intersection)
commonkeys = set.intersection(*(set(d.iterkeys()) for d in metadata.itervalues())) # * : reverse of zip
# map commonkeys -> metadata ; only if metadata also equal
commonmetadata = dict(set.intersection(*(set(d.iteritems()) for d in metadata.itervalues())))
try:
    del commonmetadata['jkBlocks'] # remove metadata that is no longer valid for the eval-results
except KeyError:
    pass

# mykeys: all the metadata-keys that have set values at least in one subdirectory.
# this means that for some subdirectories that metadata may not be set. In this case use
# a default value "NOVALUE"
mykeys = set.union(*(set(d.iterkeys()) for d in metadata.itervalues())) # * : reverse of zip
for sd, meta in metadata.items():
    for k in mykeys:
        if not k in meta.keys():
            meta[k] = "NV"
    
#print commonmetadata
mapKeySubdirValue = { k : {sd : metadata[sd][k] for sd in output_directories if sd in metadata.keys()} for k in mykeys }
multivalKeys = [k for k in mykeys if len(set(mapKeySubdirValue[k].values())) > 1]

# add "bc" explicitly as a multivalKey --> files will have "bcaveraged" written on them
if "bc" not in multivalKeys:
    multivalKeys.append("bc")

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

for tup in multivalValueTuples:
    if len(multivalKeys) == 0:
        relsubdirs = [s for s in output_directories if metadata.has_key(s)]              # take all with metadata
    else:
        def allKeysMatched(subdir):
            matched = True
            for k,v in zip(multivalKeys, tup):
                if metadata[subdir][k] != v:
                    matched = False
                    break
            return matched
        relsubdirs = filter(allKeysMatched, output_directories)

    print  "_".join(["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)])
    #now relsubdirs only differ in the values of the chosen variable
    #
    #map: observableName -> variableValue -> (average, error)
    mapToResults = {}
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
        #parse results:
        if variable == 'temp':
            varvalue = 1.0 / float(metadata[sd]['beta'])
        else:
            varvalue = float(metadata[sd][variable])
        resultFile = open(sd + "/eval-results.values", 'r')
        rows = (row.strip().split() for row in resultFile if row[0] != '#')
        for row in rows:
            if row:
                obs = row[0]
                val = row[1]
                try:
                    err = row[2]
                except:
                    err = None
                mapToResults.setdefault(obs, {})[varvalue] = (val, err)
    for obsname in mapToResults.keys():
        outfname = prefix + "eval-results-" + variable + "-" + obsname + "_" + "_".join(
            ["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)]) + ".values"
        outf = open(outfname, 'w')
        outmeta = commonmetadata.copy()
        outmeta["key"] = variable
        outmeta["observable"] = obsname
        for k,v in zip(multivalKeys, tup):
            outmeta[k] = v
        for k,v in outmeta.items():
            outf.write('# %s = %s\n' % (k, v))
        outf.write('## %s,\t %s,\t error\n' % (variable, obsname))
        for varvalue in sorted(mapToResults[obsname].keys()):
            v, e = mapToResults[obsname][varvalue]
            if v != "nan":
                outf.write('%s\t%s' % (varvalue, v))
                if e is not None:
                    outf.write('\t%s' % e)
            else:
                outf.write('# %s -- no value' % varvalue)
            outf.write('\n')
        outf.close()
