#!/usr/bin/env python
import os
import sys
import os.path
import operator                 # itemgetter
import itertools
import subprocess
from glob import glob
import argparse
from scripthelpers import mkdir_p, getHeader, parseHeader



evalOnlyNew = False
#evalOnlyNew = True


def getCommonDictionary(dictionaries):
    return dict(set.intersection(
        *(set(d.iteritems()) for d in dictionaries)))

def writeMetadictToFile(filename, meta):
    with open(filename, 'w') as f:
        for key in meta:
            f.write("# %s = %s\n" % (key, meta[key]))


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
parser.add_argument("-s",  nargs='*', type=str,
                    help="list of variables for which multiple values should be discerned "+
                    "even though they normally are not; these are always added to output filenames too.")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
default_options = "--cp-auto-range 0.001 -b 1000 -j 20 --non-iterative -i 100000 --no-tau"
parser.add_argument("--options", type=str, default=default_options,
                    help="options to be passed to mrptbc, default: " + default_options)
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)

args = parser.parse_args()
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
    + ["s", "m", "dtau", "rngSeed", "N", "alpha"]

extraNonMultiVal = args.i if args.i else []
nonMultiVal.extend(extraNonMultiVal)
extraMultiVal = args.s if args.s else []
excludedSubdirs = args.x if args.x else []
evalOnlyNew = args.evalOnlyNew

if set(extraMultiVal) & set(extraNonMultiVal):
    # intersection
    raise Exception("Arguments to options -i and -s are not disjoint")


# subdir candidates: only take directories in the top level
files_depth_1 = glob( "./*"  )
subdir_candidates = filter(lambda f: os.path.isdir(f), files_depth_1)
print "collected subdir_candidates"
sys.stdout.flush()

subdirs = [ f for f in subdir_candidates if (f.startswith(prefix) or
                                            f.startswith("./" + prefix)) ]

print "pruned subdir_candidates for prefix"
sys.stdout.flush()

# remove subdirs not containing time series files
subdirs = [ sd for sd in subdirs if glob(sd + "/p*/*.series") ]
print "pruned subdirs not containing time series files"
sys.stdout.flush()




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

# turn into list of lists, one sub-list for each group

# [ [("pbc", "subdir1"), ("apbc-x", subdir2), ...], ...]

bc_subdir_groups = [ [ (entry[1], entry[0]) for entry in entries ]
                     for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(2)) ]




#evaluate subdirectory-groups
#############################


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
    """if possible call mrpt for subdirs in tuple_list, return output_directory if not skipped from evaluation
    
       tuple_list : [(bc, subdir), ...]. return output directory if mrptbc has been run for those """

    print tuple_list,
    
    if len(tuple_list) != 4:
        print "=> faulty group, skip "
        return ""
        
    four_subdirs = [subdir for (bc, subdir) in tuple_list]

    # already done earlier: removed directories not containing any time series files

    output_directory = tuple_list[0][1].replace("bc"+tuple_list[0][0], "bc"+"averaged")
    if evalOnlyNew and glob(output_directory + "/mrpt-*.values"):
        print "already evaluated: skip, but still take into account:", output_directory
        return output_directory
    mkdir_p(output_directory)

    print "=> evaluate"

    # generate info.dat file with common metadata
    bc_infodat_files = ["%s/info.dat" % subdir for subdir in four_subdirs]
    bc_infodat_meta  = [parseHeader(getHeader(f)) for f in bc_infodat_files]
    combined_metadata = getCommonDictionary(bc_infodat_meta)
    combined_metadata["bc"] = "averaged"
    writeMetadictToFile("%s/info.dat" % output_directory, combined_metadata)

    # run mrptbc
    
    commandline = "mrptbc " + options + " --outputDirectory " + output_directory

    for bc, sd in tuple_list:
        commandline += " --info-" + bc + " " + sd + "/info.dat"

    for bc, sd in tuple_list:
        commandline += " %s/p*/associatedEnergy.series %s/p*/normMeanPhi.series" % (sd, sd)
        
    commandline += " ; exit 0"
    print commandline
    print output_directory
    stdout_and_stderr = subprocess.check_output(commandline, shell=True,
                                                stderr=subprocess.STDOUT)
    print stdout_and_stderr
    return output_directory


# mrpt already runs parallelized, only run one job at a time
output_directories = map(eval_bc_subdirs, bc_subdir_groups)
output_directories = [od for od in output_directories if od != '']





# collect mrptbc results and metadata from the output directories
#################################################################

# helper:
def addControlParameterCount(meta_dict):
    if "controlParameterValues" in meta_dict:
        meta_dict["controlParameterCount"] = str(len(meta_dict["controlParameterValues"].split()))
    return meta_dict

# map: subdirectory -> metadata dictionary [for replica exchange simulations: count controlParameterValues]
metadata = {sd: addControlParameterCount(parseHeader(getHeader(sd + "/info.dat"))) for sd in output_directories}
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
mapKeySubdirValue = { k : {sd : metadata[sd][k] for sd in output_directories if sd in metadata.keys()} for k in mykeys }
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
        print relsubdirs
    assert len(relsubdirs) == 1

    multivalString =  "_".join(["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)])
    print multivalString

    sd = relsubdirs[0]

    mrpt_direct_files = glob(sd + "/mrpt-direct-*.values")
    mrpt_files = [f for f in glob(sd + "/mrpt-*.values") if f not in mrpt_direct_files]

    def collect_mrpt_file(filename, mrpt_prefix):
        # get observable name
        observable_name = parseHeader(getHeader(filename))["observable"]
        
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
        collect_mrpt_file(f, "mrpt-bcaveraged-")
    for f in mrpt_direct_files:
        collect_mrpt_file(f, "mrpt-direct-bcaveraged-")


