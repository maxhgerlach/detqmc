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


def is_one_of_these_files(file_in_question, list_of_files):
    result = False
    for other_file in list_of_files:
        if os.path.samefile(file_in_question, other_file):
            result = True
            break
    return result


    
# parse command line arguments and start program
################################################
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", type=str, default="",
                    help="consider directories with this prefix for evaluation")
# always use the control parameter ("r") as variable
parser.add_argument("-i",  nargs='*', type=str,
                    default=[],
                    help="list of additional variables not to be discerned if"+
                    "they occur with multiple values")
parser.add_argument("-n", nargs='*', type=str,
                    default=[],
                    help="list of observables not to collect")
parser.add_argument("-s",  nargs='*', type=str,
                    help="list of variables for which multiple values should be discerned "+
                    "even though they normally are not; these are always added to output filenames too.")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
default_options = "--direct --cp-auto-range 0.001 -b 1000 -j 20 --non-iterative -i 100000 --no-tau"
parser.add_argument("--options", type=str, default=default_options,
                    help="options to be passed to mrpt, default: " + default_options)
parser.add_argument("--simindexjoined", dest='simindexjoined', action='store_true', help="only consider subdirectories with data joined from multiple simulations (by joinall.py, jointimeseries)")
parser.set_defaults(simindexjoined=False)
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)
parser.add_argument("--evalOnlyRecent", dest='evalOnlyRecent', action='store_true', help="skip already evaluated subdirectories, if data is older than 24 hours")
parser.add_argument("--no-evalOnlyRecent", dest='evalOnlyRecent', action='store_false', help="do not skip already evaluated subdirectories, if data is older than 24 hours")
parser.set_defaults(evalOnlyRecent=False)

args = parser.parse_args()
options = args.options
prefix = args.prefix
nonMultiVal = ["ARMA_VERSION",
               "buildDate",
               "buildTime",
               "gitBranch",
               "gitRevisionHash",
               "totalWallTimeSecs",
               "averageAcceptedWolffClusterSize",
               "globalShiftAccRatio",
               "wolffClusterUpdateAccRatio",
               "wolffClusterShiftUpdateAccRatio",
               "controlParameterValues",
               "controlParameterCount",
               "sweeps", "thermalization",
               "sweepsDone", "sweepsDoneThermalization",
               "turnoffFermionMeasurements",
               "join-total-samples",
               "join-subsample",
               "join-discard",
               "join-read"] \
    + ["m", "dtau", "rngSeed", "N", "alpha"]

extraNonMultiVal = args.i if args.i else []
nonMultiVal.extend(extraNonMultiVal)
nonCollectObs = args.n if args.n else []
extraMultiVal = args.s if args.s else []
excludedSubdirs = args.x if args.x else []
only_simindexjoined = args.simindexjoined
evalOnlyNew = args.evalOnlyNew
evalOnlyRecent = args.evalOnlyRecent

if set(extraMultiVal) & set(extraNonMultiVal):
    # intersection
    raise Exception("Arguments to options -i and -s are not disjoint")


# subdir candidates: only take directories in the top level
files_depth_1 = glob( "./*"  )
subdir_candidates = filter(lambda f: os.path.isdir(f), files_depth_1)
print "collected subdir_candidates"
sys.stdout.flush()

def is_simindexjoined_subdir(sd):
    header = parseHeader(getHeader(sd + "/info.dat"))
    if not header is None:
        if "simindex" in header and header["simindex"] == "joined":
            return True
    return False

subdirs = [ f for f in subdir_candidates if ((f.startswith(prefix) or
                                              f.startswith("./" + prefix)) and
                                             (not only_simindexjoined
                                              or is_simindexjoined_subdir(f)) and
                                             not is_one_of_these_files(f, excludedSubdirs)) ]
print "pruned subdir_candidates for prefix and only_simindexjoined and excludedSubdirs"
# print subdirs
# exit()

sys.stdout.flush()


# potentially append "_" for the output prefix
if prefix != "":
    if prefix[-1:] == "-":
        prefix = prefix[:-1] + "_" # replace "-" in tail by "_"
    elif prefix[-1:] != "_":
        # append a "_" if it is not there already
        prefix = prefix + "_"



#evaluate subdirectories
########################
maindirectory = os.getcwd()
def eval_subdir(sd):
    "if possible call mrpt for sd, otherwise return it as skipped"
    print sd,
    if evalOnlyRecent and os.path.exists(sd + "/mrpt-dos.dat") and is_file_older_than_a_day(sd + "/info.dat"):
        print "=> already evaluated and data is older than 24 hours: skip evaluation, but take into account"
        sys.stdout.flush()
        return ''
    elif evalOnlyNew and os.path.exists(sd + "/mrpt-dos.dat"):
        print "=> already evaluated: skip evaluation, but take into account"
        sys.stdout.flush()
        return ''
    elif not glob(sd + "/p*/*.series"):
        print "=> no time series found: skip"
        sys.stdout.flush()
        return sd
    else:
        os.chdir(sd)
        print "=> evaluate"
        sys.stdout.flush()
        stdout_and_stderr = subprocess.check_output("mrpt " + options + " --info ./info.dat --savez z.dat --loadz z.dat  ./p*/associatedEnergy.series ./p*/normMeanPhi.series; exit 0",
                                                    shell=True,
                                                    stderr=subprocess.STDOUT)
        print stdout_and_stderr
        os.chdir(maindirectory)
        return ''

# mrpt already runs parallelized, only run one job at a time
skipped_subdirs = map(eval_subdir, subdirs)
for sd in skipped_subdirs:
    if sd != '':
        subdirs.remove(sd)



# collect mrpt results and metadata
###################################

# helper:
def addControlParameterCount(meta_dict):
    if "controlParameterValues" in meta_dict:
        meta_dict["controlParameterCount"] = str(len(meta_dict["controlParameterValues"].split()))
    return meta_dict

# map: subdirectory -> metadata dictionary [for replica exchange simulations: count controlParameterValues]
metadata = {sd: addControlParameterCount(parseHeader(getHeader(sd + "/info.dat"))) for sd in subdirs}
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

variable = commonmetadata["controlParameterName"]

# mykeys: all the metadata-keys that have set values atleast in one subdirectory.
# this means that for some subdirectories that metadata may not be set. In this case use
# a default value "NOVALUE"
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




commonHeader = "".join(['# %s = %s\n' % (k, v) for k,v in commonmetadata.items()])

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

    if len(relsubdirs) == 0:
        # no data matching this key tuple
        continue

    # MRPT: we should only have a single subdir matching each
    # multivalValueTuple, the mrpt results already collected data for
    # different control parameter values
    #
    # NOTE: We might want to
    # generalize this at some point, to be able to combine data from
    # separate PT runs for different parameter ranges -- this is not
    # essential, though
    if len(relsubdirs) != 1:
        print relsubdirs, "not good for", tup
    assert len(relsubdirs) == 1

    multivalString =  "_".join(["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)])
    print multivalString

    sd = relsubdirs[0]

    mrpt_direct_files = glob(sd + "/mrpt-direct-*.values")
    mrpt_files = [f for f in glob(sd + "/mrpt-*.values") if f not in mrpt_direct_files]

    def collect_mrpt_file(filename, mrpt_prefix):
        # get observable name
        observable_name = parseHeader(getHeader(filename))["observable"]

        if observable_name in nonCollectObs:
            # skip this observable
            return
        
        output_filename = prefix + mrpt_prefix + variable + "-" + observable_name + "_" + \
                          multivalString + ".values"

        with open(output_filename, 'w') as output_file:
            # prepend commonmetadata, add key = variable
            output_file.write(commonHeader)
            for k,v in zip(multivalKeys, tup):
                output_file.write("# %s = %s" % (k, v) + "\n")
            output_file.write("# key = " + variable + "\n")

            # copy rest of file contents
            with open(filename, "r") as input_file:
                for line in input_file:
                    output_file.write(line) 

    for f in mrpt_files:
        collect_mrpt_file(f, "mrpt-")
    for f in mrpt_direct_files:
        collect_mrpt_file(f, "mrpt-direct-")


