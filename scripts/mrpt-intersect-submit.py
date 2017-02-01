#!/usr/bin/env python
import numpy as np
import os
import errno
import sys
import os.path
import re
import operator                 # itemgetter
import itertools
import subprocess
from subprocess import check_output
from glob import glob
import argparse
from scripthelpers import queryRunningJobs


# use mrpt to find intersection points for various quantites
#
# this script is for use on the cluster, submitting an individual
# cheops job for each requested intersection
#
# the results should then be gathered by a call to mrpt-collect-intersect.py --evalOnlyNew .







# template for jobs, strings preceded by % to be replaced before submitting
jobTemplate = r"""#!/bin/bash -l
#SBATCH --cpus-per-task=%cores
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1024
#SBATCH --time=01:00:00
#SBATCH --account=AG-Trebst
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mgerlac2@uni-koeln.de

module add intel/15.0 mkl/11.0

# This catches the output from an external Python script which 
# prints the number of seconds this job has left to run. 
export PBS_WALLTIME=$(slurm_seconds_left.py)

export WORKDIR=%maindir
cd $WORKDIR

# output start time
echo start time
date
echo

TIME1=$(date +%s)

env OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK %commandline
RETVAL=$?

TIME2=$(date +%s)

# output finish time
echo finish time
date
echo

exit "$RETVAL"
"""


def setup_job(jobname, output_directory, commandline, cores):
    job = jobTemplate.replace("%maindir", os.getcwd()) \
                     .replace("%cores", str(cores)) \
                     .replace("%commandline", commandline)

    jobfilename = "%s/%s.sh" % (output_directory, jobname)
    # %j is replaced by the job allocation number:
    outputfile = "%s/%s_output.%%j.log" % (output_directory, jobname) 
    errorfile  = "%s/%s_error.%%j.log" % (output_directory, jobname)

    qsubcommand = "sbatch --job-name=%s --output=%s --error=%s %s" % (
        jobname, outputfile, errorfile, jobfilename)

    with open(jobfilename, 'w') as jobfile:
        jobfile.write(job)
    # submit job
    try:
        stdout_and_stderr = check_output(qsubcommand, shell=True, stderr=subprocess.STDOUT)
        print stdout_and_stderr
    except subprocess.CalledProcessError, e:
        print "Command:", e.cmd
        print "Error:", e.output
        raise
    
    




def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


evalOnlyNew = False
# evalOnlyNew = True


# return list of lines starting with #, with the # removed, from file
# only beginning
def getHeader(f):
    header = []
    try:
        for line in open(f).xreadlines():
            if line[0] == "#":
                header.append(line[1:])
            else:
                break
        return header
    except:
        print "can't open", f


def parseHeader(header):
    d = {}
    # Regexp to match L = 256 and similar:
    p = re.compile(r"\s* ([^\s=]+) \s* = \s* (.+)", re.VERBOSE)
    try:
        for line in header:
            if line[0] == "#":
                continue
            m = p.match(line)
            try:
                d[m.group(1)] = m.group(2).strip()
            except:
                pass
        return d
    except:
        return None


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
parser.add_argument("--cores", type=int, default=4,
                    help="number of cpu cores to run each mrpt instance on, default: 4")

parser.add_argument("-f", nargs='*', type=str,
                    help="list of 4 column text files with entries: [L1 L2 cpMin cpMax] -- search intersection for these lattice sizes in that cp range; each file needs to come with a metadata block '# beta = ...' etc. uniquely specifying a group of dictionaries")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
default_options = "-b 1000 -j 20 --non-iterative -i 100000 --no-tau"
parser.add_argument("--options", type=str, default=default_options,
                    help="options to be passed to mrpt-find-intersect, default: '" + default_options + "'\n" +
                    "We will call mrpt-find-intersect separately for Binder cumulants and scaled susceptibilities [TODO: allow paramters to customize this behavior].")
parser.add_argument("--simindexjoined", dest='simindexjoined', action='store_true', help="only consider subdirectories with data joined from multiple simulations (by joinall.py, jointimeseries)")
parser.set_defaults(simindexjoined=False)
parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
parser.set_defaults(evalOnlyNew=False)

args = parser.parse_args()
cores = args.cores
options = args.options
prefix = args.prefix

excludedSubdirs = args.x if args.x else []
only_simindexjoined = args.simindexjoined
evalOnlyNew = args.evalOnlyNew
controlFiles = args.f

# print excludedSubdirs

# subdir candidates: only take directories in the top level
files_depth_1 = glob( "./*"  )
subdir_candidates = filter(lambda f: os.path.isdir(f), files_depth_1)
print "collected subdir_candidates"
sys.stdout.flush()

subdirs = [ f for f in subdir_candidates if ((f.startswith(prefix) or
                                              f.startswith("./" + prefix)) and
                                             (not (f in excludedSubdirs or
                                                   (f.startswith("./") and f[2:] in excludedSubdirs))))]
print "pruned subdir_candidates for prefix and excludedSubdirs"
sys.stdout.flush()

# remove subdirs not containing time series files
subdirs = [ sd for sd in subdirs if glob(sd + "/p*/*.series") ]
print "pruned subdirs not containing time series files"
# for subdir in subdirs:
#     print subdir
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
        if only_simindexjoined and "simindex" in header and header["simindex"] != "joined":
            continue
        else:
            infodata[sd] = header
subdirs = infodata.keys()
        


# for subdir in infodata.keys():
#     print subdir


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
exclude_from_metadata = ["L", "N", # in one group we want to have all possible lattice sizes
                         # compilation related keys:
                         "ARMA_VERSION", "BOOST_LIB_VERSION", "buildDate", "buildHost", "buildTime", "cppflags", "cxxflags", "gitBranch", "gitRevisionHash",
                         # if some measured quantities are not present in part of the data, it should not matter
                         "turnoffFermionMeasurements",
                         # unimportant for the data:
                         'saveConfigurationStreamInterval', 'dumpGreensFunction',
                         # rng seed will be different:
                         "rngSeed",
                         # control parameter values do not matter
                         "controlParameterCount", "controlParameterValues",
                         # target sweeps could have been modified for some
                         "sweeps", "thermalization",
                         # current state of simulation will be different
                         "sweepsDone", "sweepsDoneThermalization", "totalWallTimeSecs",
                         # dynamical data will be different
                         "globalShiftAccRatio", "wolffClusterUpdateAccRatio", "averageAcceptedWolffClusterSize",
                         "wolffClusterShiftUpdateAccRatio", "averageAcceptedWolffClusterSize",
                         # metadata from joinall.py
                         "join-total-samples", "join-subsample", "join-discard", "join-read"]
for sd, metadict in infodata.iteritems():
    L = metadict["L"]
    interesting_metadict = { k: metadict[k] for k in metadict if k not in exclude_from_metadata }

    #debug
    # if metadict["beta"] == "20" and metadict["L"] in ("14", "16"):
    #     print
    #     print interesting_metadict
    #     print

    groupable_list.append( (sd, L, sorted(interesting_metadict.iteritems())) )

groupable_list.sort(key = operator.itemgetter(2))

# turn into list of lists, one sub-list for each group

# [ [(L1, "subdir1"), (L2, subdir2), ...], ...]

L_subdir_groups = [ [ (entry[1], entry[0]) for entry in entries ]
                    for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(2)) ]

# for l_sd in L_subdir_groups:
#     print l_sd
# exit()


#evaluate subdir-groups
#######################


def find_intersection_for_subdirs(tuple_list):
    """if possible call mrpt-find-intersect via a SLURM job for subdirs in
       tuple_list, return output_directory if not skipped from
       evaluation
    
       tuple_list : [(L, subdir), ...]. return output directory if successful

    """
    
    print tuple_list

    if len(tuple_list) < 2:
        print "=> too few subdirectories, skip"
        return ""

    my_subdirs = [subdir for (L, subdir) in tuple_list]
    L_to_subdirs = {int(L): subdir for (L, subdir) in tuple_list}
    my_subdirs = L_to_subdirs.values()
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
    # intersections [the same control files for all quantities]
    L1_L2_cpMin_cpMax = np.loadtxt(my_cf, ndmin=2)

    output_directory = prefix + "mrpt-find-intersect"
    for key, value in my_cf_meta.items():
        output_directory += "_" + key + value
    mkdir_p(output_directory)
    print output_directory

    print "=> setting up evaluation"

    # generate info.dat file with common metadata
    L_infodat_files = ["%s/info.dat" % subdir for subdir in my_subdirs]
    L_infodat_meta  = [parseHeader(getHeader(f)) for f in L_infodat_files]
    combined_metadata = getCommonDictionary(L_infodat_meta)
    writeMetadictToFile("%s/info.dat" % output_directory, combined_metadata)

    runningJobs = queryRunningJobs("$USER").values()
    
    # submit mrpt-find-intersect calls
    for L1, L2, cpMin, cpMax in L1_L2_cpMin_cpMax:
        L1 = int(L1)
        L2 = int(L2)
        if evalOnlyNew and glob(output_directory + "/mrpt-*-intersect-l%dl%d.dat" % (L1,L2)):
            print "already evaluated: skip submitting new job: (%d, %d)" % (L1,L2)
            continue

        try:
            sd1 = L_to_subdirs[L1]
            sd2 = L_to_subdirs[L2]
        except KeyError:
            print "no data for (%d, %d), skip" % (L1,L2)
            continue
        info1 = sd1 + "/info.dat"
        info2 = sd2 + "/info.dat"
        zfile1 = sd1 + "/z.dat"
        zfile2 = sd2 + "/z.dat"        

        for quantity in ["BinderRatio", "ScaledKTSusceptibility"]:
            print "Going to find", quantity, "intersection for L1=", L1, ", L2=", L2

            jobname = output_directory + "_quantity_" + quantity + "_L1_" + str(L1) + "_L2_" + str(L2)

            if jobname in runningJobs:
                print "==> job already running"
                continue

            commandline = "mrpt-find-intersect " + options \
                          + " --quantity " + quantity \
                          + " --outputDirectory " + output_directory \
                          + " --info1 " + info1 + " --info2 " + info2 \
                          + " --loadz1 " + zfile1 + " --savez1 " + zfile1 \
                          + " --loadz2 " + zfile2 + " --savez2 " + zfile2 \
                          + " --cp-range %f %f" % (cpMin, cpMax)
            for sd in [sd1, sd2]:
                commandline += " %s/p*/associatedEnergy.series %s/p*/normMeanPhi.series" % (sd, sd)

            setup_job(jobname, output_directory, commandline, cores)
    return output_directory


for tuple_list in L_subdir_groups:
    find_intersection_for_subdirs(tuple_list)


