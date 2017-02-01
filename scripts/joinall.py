#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import os.path
import operator                 # itemgetter
import itertools
import subprocess
from glob import glob
import argparse
from scripthelpers import mkdir_p, getHeader, parseHeader, is_file_older_than_a_day, is_file_more_than_a_day_older_than_file

def getCommonDictionary(dictionaries):
    return dict(set.intersection(
        *(set(d.iteritems()) for d in dictionaries)))

def writeMetadictToFile(filename, meta):
    with open(filename, 'w') as f:
        for key in meta:
            f.write("# %s = %s\n" % (key, meta[key]))

def is_there_new_data(input_directories):
    result = False
    for sd in input_directories:
        file_to_check = os.path.join(sd, "info.dat")
        if not is_file_older_than_a_day(file_to_check):
            result = True
            break
    return result

def data_mucher_newer(input_directories, output_info_path):
    """check whether the input data is at most a day newer than the output
    data, i.e. whether there is some input data compared to which
    the output info path is at least a day older
    """
    result = False
    for sd in input_directories:
        file_to_check = os.path.join(sd, "info.dat")
        if is_file_more_than_a_day_older_than_file(output_info_path, file_to_check):
            # output info is more than a day older than this input info
            result = True
            break
    return result



# parse command line arguments and start program
################################################
parser = argparse.ArgumentParser()

parser.add_argument("--prefix", type=str, default="",
                    help="consider directories with this prefix to be joined")
parser.add_argument("--outputPrefix", type=str, default="",
                   help="optionally prepend this to the output directory names (this is hacked in, maybe flaky)")
parser.add_argument("-x", nargs='*', type=str,
                    help="list of subdirectories to exclude from the evaluation")
default_options = ""
parser.add_argument("--options", type=str, default=default_options,
                    help="extra options to be passed to jointimeseries, default: " + default_options)
parser.add_argument("--joinOnlyNew", dest='joinOnlyNew', action='store_true', help="skip subdirectory groups already joined before")
parser.add_argument("--no-joinOnlyNew", dest='joinOnlyNew', action='store_false', help="do not skip subdirectory groups already joined before")
parser.add_argument("--joinOnlyRecent", dest='joinOnlyRecent', action='store_true', help="skip subdirectory groups already joined before, if data is older than 24 hours")
parser.add_argument("--no-joinOnlyRecent", dest='joinOnlyRecent', action='store_false', help="do not skip subdirectory groups already joined before, if data is older than 24 hours")
parser.add_argument("--joinOnlyNecessary", dest='joinOnlyNecessary', action='store_true', help="skip subdirectory groups already joined before, if the joined data is more than 24 hours older than the input")
parser.add_argument("--no-joinOnlyNecessary", dest='joinOnlyNecessary', action='store_false', help="do not skip subdirectory groups already joined before, if the joined data is more than 24 hours older than the input")
parser.set_defaults(joinOnlyNew=False)
parser.set_defaults(joinOnlyRecent=False)
parser.set_defaults(joinOnlyNecessary=False)

args = parser.parse_args()
joinOnlyNew = args.joinOnlyNew
joinOnlyRecent = args.joinOnlyRecent
joinOnlyNecessary = args.joinOnlyNecessary
options = args.options
prefix = args.prefix
outputPrefix = args.outputPrefix
excludedSubdirs = args.x if args.x else []

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
           and ("simindex" in f and not "simindexjoined" in f)  # we want to have the simindex written into the directory name
           ]


#collect info.dat contents, maybe add simindex = 0
##################################################
infodata = {}
for sd in subdirs:
    header = parseHeader(getHeader(sd + "/info.dat"))
    if not header is None:
        infodata[sd] = header
        if not 'simindex' in infodata[sd]:
            infodata[sd]['simindex'] = '0'  # default simindex: 0
infodata = {k:v for k,v in infodata.iteritems() if v is not None}


# We need to group the subdirectories in groups that differ only by the simindex
################################################################################
groupable_list = []
# keys in the following links are removed from the metadata before sorting by the metadata
exclude_from_metadata = ["simindex", # in one group we want to have all simindex
                         # compilation related keys:
                         "ARMA_VERSION", "BOOST_LIB_VERSION", "buildDate", "buildHost", "buildTime", "cppflags", "cxxflags", "gitBranch", "gitRevisionHash",
                         # if some measured quantities are not present in part of the data, it should not matter
                         "turnoffFermionMeasurements",
                         # rng seed will be different:
                         "rngSeed",
                         # target sweeps could have been modified for some
                         "sweeps", "thermalization",
                         # current state of simulation will be different
                         "sweepsDone", "sweepsDoneThermalization", "totalWallTimeSecs",
                         # dynamical data will be different
                         "globalShiftAccRatio", "wolffClusterUpdateAccRatio", "averageAcceptedWolffClusterSize",
                         "wolffClusterShiftUpdateAccRatio", "averageAcceptedWolffClusterSize"]
for sd, metadict in infodata.iteritems():
    simindex = metadict["simindex"]
    interesting_metadict = { k: metadict[k] for k in metadict if k not in exclude_from_metadata }
    groupable_list.append( (sd, simindex, sorted(interesting_metadict.iteritems())) )

groupable_list.sort(key = operator.itemgetter(2))

# turn into list of lists, one sub-list for each group

simindex_subdir_groups = [ [ (entry[1], entry[0]) for entry in entries ]
                           for metalist, entries in itertools.groupby(groupable_list, operator.itemgetter(2)) ]
# for group in simindex_subdir_groups:
#     if len(group) != 4:
#         print( group )
# exit()


#join data from subdirectory-groups
###################################

def check_time_series_lengths(directory):
    """Check if all the time series files in @directory contain the same
    number of lines.  This is a consistency check."""
    # $ wc -l associatedEnergy.series 
    # 23227 associatedEnergy.series
    counts = []
    for f in glob(os.path.join(directory, "*.series")):
        output = subprocess.check_output(['wc', '-l', f])
        count = int(output.split()[0])
        assert count > 0
        counts.append(count)
    # check if all elements in counts are equal (or the list is empty)
    return len(set(counts)) <= 1

def join_simindex_subdirs(tuple_list):
    " tuple_list : [(simindex, subdir), ...]. return output directory if jointimeseries has been run "
    input_directories = [sd for simindex, sd in tuple_list]
    output_directory = tuple_list[0][1].replace("simindex"+tuple_list[0][0], "simindex"+"joined")
    if output_directory.startswith("./"):
        output_directory = "./" + outputPrefix + output_directory[2:]
    else:
        output_directory = outputPrefix + output_directory

    output_info = os.path.join(output_directory, "info.dat")
        
    if joinOnlyRecent and os.path.exists(output_info) and not is_there_new_data(input_directories):
        print( "already have joined data and newest data is older than 24 hours: skip", output_directory )
        return output_directory

    if joinOnlyNew and os.path.exists(output_info):
        print( "already have joined data: skip", output_directory )
        return output_directory

    if joinOnlyNecessary and os.path.exists(output_info) and not data_mucher_newer(input_directories,
                                                                                   output_info):
        print( "already have joined data, which is older or less than 24h newer than the data: skip", output_directory )
        return output_directory
    
    mkdir_p(output_directory)
    
    # join timeseries
    commandline = "jointimeseries " + options + " --outputDirectory " + output_directory
    for simindex, sd in tuple_list:
        commandline = commandline + " " + sd
    commandline = commandline + " ; exit 0"
    print( commandline )
    print( output_directory )
    pipes = subprocess.Popen(commandline, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE )
    std_out, std_err = pipes.communicate()
    print( std_out )
    if std_err:
        print( "Failure in jointimeseries", file=sys.stderr )
        print( "commandline was: " + commandline, file=sys.stderr )
        print( std_err, file=sys.stderr )
        print( file=sys.stderr )
    if not check_time_series_lengths(output_directory):
        print( "Time series length mismatch in joinall.py output directory " + output_directory, file=sys.stderr )
        print( file=sys.stderr )
    return output_directory

import multiprocessing as mp
# pool_size = min(4, os.getenv("SLURM_CPUS_PER_TASK"))
pool_size = os.getenv("SLURM_CPUS_PER_TASK")
if pool_size is None:
    pool_size = min(4, mp.cpu_count())
else:
    pool_size = int(pool_size)
pool = mp.Pool(processes=pool_size)

output_directories = pool.map(join_simindex_subdirs, simindex_subdir_groups, 1)
pool.close()
pool.join()



output_directories = [od for od in output_directories if od != '']

print( )
print( "joined data into directories" )
#print( output_directories )
print( )

# join info.dat files of sublevels into PT top folders
files_depth_1 = glob( "./*"  )
toplevel_subdir_candidates = filter(lambda f: os.path.isdir(f), files_depth_1)
toplevel_subdirs = [ f for f in toplevel_subdir_candidates if ((f.startswith(prefix) or
                                                                f.startswith("./" + prefix)) and
                                                               ("simindexjoined" in f) and
                                                               not glob(f + "/*.series")) ]
for tlsd in toplevel_subdirs:
    tlsd_contents = os.listdir(tlsd)
    my_subdirs = filter(lambda f: os.path.isdir(os.path.join(tlsd,f)), tlsd_contents)
    if not my_subdirs:
        # a directory without subdirectories: nothing useful here
        continue
    print( tlsd )
    
    infodata_list = [parseHeader(getHeader(os.path.join(tlsd,sd,"info.dat")))
                     for sd in my_subdirs]
    infodata_list = [d for d in infodata_list if d is not None]
    if len(infodata_list) == 0:
        print( "  no info.dat files" )
        continue
    # common info.dat // so aehnlich
    commoninfodata = dict(
        set.intersection( *(set(d.iteritems()) for d in infodata_list) )
    )
    with open( os.path.join(tlsd, "info.dat"), "w" ) as outf:
        for k,v in commoninfodata.items():
            outf.write('# %s = %s\n' % (k, v))
