#!/usr/bin/env python
import numpy as np
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


def dictContainedInDict(smallDict, largeDict):
    # convert to item pairs and check for containment
    return all(item in largeDict.items() for item in smallDict.items())


# parse command line arguments and start program
################################################
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", type=str, default="",
                    help="consider directories with this prefix for evaluation")
parser.add_argument("variable", type=str,
                    help="quantity to use as primary independent variable")
parser.add_argument("-i",  nargs='+', type=str,
                    default=[],
                    help="list of additional variables not to be discerned if" +
                    "they occur with multiple values")
parser.add_argument("-f", nargs='*', type=str,
                    help="list of 4 column text files with entries: [L1 L2 cpMin cpMax] -- search intersection for these lattice sizes in that cp range; each file needs to come with a metadata block '# beta = ...' etc. uniquely specifying a group of dictionaries")
parser.add_argument("-s",  nargs='*', type=str,
                    help="list of variables for which multiple values should be discerned "+
                    "even though they normally are not; these are always added to output filenames too.")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
default_options = "-b 1000 -j 20 --non-iterative -i 100000 --no-tau"
parser.add_argument("--options", type=str, default=default_options,
                    help="options to be passed to mrpt-binderratio-intersect, default: " + default_options)
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)

args = parser.parse_args()
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
    + ["s", "m", "dtau", "rngSeed", "N", "alpha"]

extraNonMultiVal = args.i if args.i else []
nonMultiVal.extend(extraNonMultiVal)
extraMultiVal = args.s if args.s else []
excludedSubdirs = args.x if args.x else []
evalOnlyNew = args.evalOnlyNew
controlFiles = args.f

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


# collect controlFiles metadata
controlFilesMetaData = {}
for cf in controlFiles:
    header = parseHeader(getHeader(cf))
    if not header is None:
        controlFilesMetaData[cf] = header
    else:
        print "control file", cf, "does not contain metadata"
        



# group subdirectories such that the metadata for one group only differs in lattice size L
##########################################################################################
groupable_list = []
# keys in the following links are removed from the metadata before sorting by the metadata
exclude_from_metadata = ["L", "N", "bc", # we want to group by L and bc
                         # compilation related keys:
                         "ARMA_VERSION", "BOOST_LIB_VERSION", "buildDate", "buildHost", "buildTime", "cppflags", "cxxflags", "gitBranch", "gitRevisionHash",
                         # rng seed will be different:
                         "rngSeed",
                         # control parameter values do not matter
                         "controlParameterCount", "controlParameterValues",
                         # target sweeps could have been modified for some
                         "sweeps", "thermalization",
                         # current state of simulation will be different
                         "sweepsDone", "sweepsDoneThermalization", "totalWallTimeSecs"]
for sd, metadict in infodata.iteritems():
    L = metadict["L"]
    bc = metadict["bc"]
    interesting_metadict = { k: metadict[k] for k in metadict if k not in exclude_from_metadata }
    groupable_list.append( (sd, L, bc, sorted(interesting_metadict.iteritems())) )

# first sort by non-L, non-bc metadata: items of a L-bc-group will be consecutive
groupable_list.sort(key = operator.itemgetter(3))

# turn into list of lists, one sub-list for each group
# each sub-list is sorted by L
# [ [(L1, bca, "subdir1a"), (L1, bcb, "subdir1b"), ..., (L2, bca, "subdir2a"), (L2, bcb, "subdir2b")], ...]

L_bc_subdir_groups = [ sorted([ (int(entry[1]), entry[2], entry[0]) for entry in entries ], key = operator.itemgetter(0) )
                       for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(3)) ]



#evaluate subdir-groups
#######################

all_bc = ["pbc", "apbc-x", "apbc-y", "apbc-xy"]

def find_intersection_for_subdirs(tuple_list):
    """if possible call mrpt-binderratio-intersect for subdirs in tuple_list,
       return output_directory if not skipped from evaluation
    
       tuple_list : [(L, bc, subdir), ...]. return output directory if successful """

    if len(tuple_list) % 4 != 0:
        print "=> number of subdirectories should be multiple of 4, skipping"
        return ""
    
    # map: L -> bc -> subdir
    map_L_bc_subdir = { }
    for L, bc, subdir in tuple_list:
        map_L_bc_subdir.setdefault(int(L), { })[bc] = subdir
    
    print map_L_bc_subdir

    for L in map_L_bc_subdir:
        if len(map_L_bc_subdir[L].keys()) != 4:
            print "Wrong number of bc subdirs for L =", L, "=> skipping"
            return ""

    my_subdirs = [subdir for (L, bc, subdir) in tuple_list]
    # already done earlier: removed directories not containing any time series files

    # Find control_file for our group of subdirectories, if we have one.
    # Otherwise, skip this directory.
    my_cf = None
    my_cf_meta = None
    for cf, meta in controlFilesMetaData.items():
        if dictContainedInDict(meta, infodata[my_subdirs[0]]):
            my_cf = cf
            my_cf_meta = meta
            break
    if my_cf is None:
        return ""
    print "control file:", my_cf
    for sd in my_subdirs[1:]:
        if not dictContainedInDict(my_cf_meta, infodata[sd]):
            print "Error: control file does not match subdirectory:", sd
            return ""

    # get information (4 columns) about where to look for L-pair
    # Binder cumulant intersections
    L1_L2_cpMin_cpMax = np.loadtxt(my_cf, ndmin=2)

    output_directory = prefix + "mrptbc-binderratio-intersect"
    for key, value in my_cf_meta.items():
        output_directory += "_" + key + value
    output_directory += "_bcaveraged"
    mkdir_p(output_directory)
    print output_directory

    print "=> evaluate"

    # generate info.dat file with common metadata
    L_infodat_files = ["%s/info.dat" % subdir for subdir in my_subdirs]
    L_infodat_meta  = [parseHeader(getHeader(f)) for f in L_infodat_files]
    combined_metadata = getCommonDictionary(L_infodat_meta)
    combined_metadata["bc"] = "averaged"
    writeMetadictToFile("%s/info.dat" % output_directory, combined_metadata)
    
    # make mrptbc-binderratio-intersect calls
    for L1, L2, cpMin, cpMax in L1_L2_cpMin_cpMax:
        L1 = int(L1)
        L2 = int(L2)
        if evalOnlyNew and glob(output_directory + "/mrptbc-binder-intersect-l%dl%d.dat" % (L1,L2)):
            print "already evaluated: skip, but still take into account: (%d, %d)" % (L1,L2)
            continue
        
        sd1_bc = {}
        sd2_bc = {}
        info1_bc = {}
        info2_bc = {}

        for bc in all_bc:
            sd1_bc[bc] = map_L_bc_subdir[L1][bc]
            sd2_bc[bc] = map_L_bc_subdir[L2][bc]          
            info1_bc[bc] = sd1_bc[bc] + "/info.dat"
            info2_bc[bc] = sd2_bc[bc] + "/info.dat"

        print "Finding Binder-ratio intersection for L1=", L1, ", L2=", L2

        commandline = "mrptbc-binderratio-intersect " + options \
                      + " --outputDirectory " + output_directory \
                      + " --cp-range %f %f" % (cpMin, cpMax)
        for bc in all_bc:
            commandline += " --info1-%s %s --info2-%s %s" % (bc, info1_bc[bc], bc, info2_bc[bc])

        for bc in all_bc:
            for sd in [sd1_bc[bc], sd2_bc[bc]]:
                commandline += " %s/p*/associatedEnergy.series %s/p*/normMeanPhi.series" % (sd, sd)
        
        commandline += " ; exit 0"

        print "commandline:", commandline

        stdout_and_stderr = subprocess.check_output(commandline, shell=True,
                                                    stderr=subprocess.STDOUT)
        print stdout_and_stderr

    return output_directory


# mrpt* already runs parallelized, only run one job at a time
output_directories = map(find_intersection_for_subdirs, L_bc_subdir_groups)
output_directories = [od for od in output_directories if od != '']



# collect Binder-ratio intersection results
###########################################

# map: output_directory -> metadata dictionary
metadata = {sd: parseHeader(getHeader(sd + "/info.dat")) for sd in output_directories}
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


# mykeys: all the metadata-keys that have set values atleast in one output directory.
# this means that for some subdirectories that metadata may not be set. In this case use
# a default value "NOVALUE"
mykeys = set.union(*(set(d.iterkeys()) for d in metadata.itervalues())) # * : reverse of zip
for sd, meta in metadata.items():
    for k in mykeys:
        if not k in meta.keys():
            meta[k] = "NV"

mapKeySubdirValue = { k : {sd : metadata[sd][k] for sd in metadata} for k in mykeys }
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
multivalValueTuples = [t for t in itertools.product(*multivalValueLists)]


controlParameterName = commonmetadata["controlParameterName"]
controlParameterNameError = controlParameterName + "Error"
commonHeader = "".join(['# %s = %s\n' % (k, v) for k,v in commonmetadata.items()])


for tup in multivalValueTuples:
    ## find all subdirectories for this multival tuple    
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
    #now relsubdirs only differ in the values of the chosen variable

    multivalString =  "_".join(["%s%s" % (k,v) for (k,v) in zip(multivalKeys, tup)])
    print multivalString

    ## find all assembled Binder crossing data
    #map: (L1, L2) -> variableValue -> (intersection control param location, error)
    mapToIntersectionResults = {}
    for sd in relsubdirs:
        if variable == 'temp':
            varvalue = 1.0 / float(metadata[sd]['beta'])
        else:
            varvalue = float(metadata[sd][variable])

        mrpt_intersect_files = glob(sd + "/mrpt*binder*intersect*.dat")

        for fname in mrpt_intersect_files:
            fdict = parseHeader(getHeader(fname))
            L1 = fdict["L1"]
            L2 = fdict["L2"]
            if controlParameterName in fdict:
                # found an intersection point
                cp = fdict[controlParameterName]
                if controlParameterNameError in fdict:
                    # found an error estimate for the intersection point
                    cpError = fdict[controlParameterNameError]
                else:
                    cpError = 0.0                    
                mapToIntersectionResults.setdefault((L1,L2), {})[varvalue] = (cp, cpError)
            else:
                continue

    ## output collected results for each pair of (L1,L2)

    for L1,L2 in mapToIntersectionResults.keys():
        
        output_filename = prefix + "mrptbc-binderratio-intersect_" + variable + "-" + controlParameterName + \
                          "-L" + L1 + "L" + L2 + \
                          "_" + multivalString + ".values"

        with open(output_filename, 'w') as output_file:
            # prepend commonmetadata, add key = variable
            output_file.write(commonHeader)
            for k,v in zip(multivalKeys, tup):
                output_file.write("# %s = %s" % (k, v) + "\n")
            output_file.write("# L1 = %s" % L1 + "\n")            
            output_file.write("# L2 = %s" % L2 + "\n")
            output_file.write("# key = " + variable + "\n")
            output_file.write("## Binder ratio intersection\n")
            output_file.write("## %s %s %s\n" % (variable, controlParameterName, controlParameterNameError))

            for varvalue in sorted(mapToIntersectionResults[L1,L2].keys()):
                loc, err = mapToIntersectionResults[L1,L2][varvalue]
                if loc != "nan":
                    output_file.write('%s\t%s' % (varvalue, loc))
                    if err is not None:
                        output_file.write('\t%s' % err)
                else:
                    output_file.write('# %s -- no value' % varvalue)
                output_file.write('\n')
            output_file.close()
