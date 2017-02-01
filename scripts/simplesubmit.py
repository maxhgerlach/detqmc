#!/usr/bin/env python
from collections import OrderedDict
import itertools
import os
import os.path
import re
import subprocess
from subprocess import check_output
import argparse
from scripthelpers import mkdir_p, parseConf

# script to submit a job to SLURM, running an executable that is
# configured by commandline arguments (no resubmitting)
    
# template for jobs, strings preceded by % to be replaced before submitting
jobTemplate = r"""#!/bin/bash -l
#SBATCH --cpus-per-task=%cores
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --constraint=%constraint
#SBATCH --mem=%mem
#SBATCH --time=%walltime
#SBATCH --time-min=%mintime
#SBATCH --account=AG-Trebst
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mgerlac2@uni-koeln.de

module add intel/15.0 mkl/11.0

# This catches the output from an external Python script which 
# prints the number of seconds this job has left to run. 
export PBS_WALLTIME=$(slurm_seconds_left.py)

export WORKDIR=%maindir/%subdir
mkdir -p $WORKDIR
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


# these do not allow multiple values
possibleJobOptions = ['jobprefix', 'onlynew', 'cores', 'mem',
                      'walltime', 'mintime', 'maindir', 'subdir',
                      'constraint', 'partition']


#man sbatch says:
#
# Acceptable  time  formats  include  "minutes", "minutes:seconds", "hours:min-
# utes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
def walltimeInSeconds(walltimeString):
    if "-" in walltimeString:
        # day included: d-HH(:MM(:SS))
        daystr, timestr = walltimeString.split("-")
        colcount = timestr.count(":")
        if colcount == 2:
            # HH:MM:SS
            ftr = [3600,60,1]
        elif colcount == 1:
            # HH:MM
            ftr = [3600,60]
        elif colcount == 0:
            # HH
            ftr = [3600]
        else:
            raise Exception("invalid time specification")      
    else:
        # no day included: ((HH:))MM(:SS)
        daystr = "0"
        timestr = walltimeString
        colcount = timestr.count(":")
        if colcount == 2:
            # HH:MM:SS
            ftr = [3600,60,1]
        elif colcount == 1:
            # MM:SS
            ftr = [60,1]
        elif colcount == 0:
            # MM
            ftr = [60]
        else:
            raise Exception("invalid time specification")            
    seconds = int(daystr) * 24 * 3600
    seconds += sum([a*b for a,b in zip(ftr, map(int,timestr.split(':')))])
    return str(seconds)


# construct jobs from jobconffilename, only submit them if no job of that name is in list runningJobs
# ignore onlynew
def process(jobconffname, runningJobs, executable="detqmcsdw"):
    jobconf = parseConf(open(jobconffname, 'r'))

    options = jobconf.keys()
    onlynew = False

    # options to be passed to simulation instances
    simoptions = [opt for opt in options if opt not in possibleJobOptions]
    # the rest, to be replaced in the jobTemplate
    joboptions = [opt for opt in options if opt in possibleJobOptions]

    # filter out option value list of the form "option, [empty]"
    # ...
    optionsForJobName = []
    for opt in simoptions:
        values = jobconf[opt]
        if values[-1] == "":
            # the last value is empty
            # -> use this opt for the jobname, remove the empty entry from jobconf
            optionsForJobName.append(opt)
            del values[-1]

    if 'maindir' in joboptions:
        maindir = jobconf["maindir"][0]
    else:
        ##TODO: write to scratch
        # maindir =  os.path.relpath(os.getcwd(),
        #                       '/home/mgerlac2/heisenberg-kitaev/sim')
        maindir = os.getcwd()
    if 'jobprefix' in joboptions:
        jobprefix = jobconf["jobprefix"][0] + "_"
    else:
        jobprefix = ""
    if 'cores' in joboptions:
        cores = jobconf['cores'][0]
    else:
        cores = "1"
    if 'constraint' in joboptions:
        constraint = jobconf['constraint'][0]
    else:
        constraint = '""'
    if 'partition' in joboptions and jobconf['partition'][0] == 'devel':
        partitionOption = "--partition=devel"
    else:
        partitionOption = ""
    if 'mintime' in joboptions:
        mintime = jobconf['mintime'][0]
    else:
        mintime = jobconf['walltime'][0]

    # additionally: specify opt1:opt2 pairs with values val1_a:val2_a, val1_b:val2_b, ...
    pairoptions = [opt for opt in simoptions if re.match(r'.+:.+', opt)]
    optionpairs = [opt.split(':', 1) for opt in pairoptions]
    pairoptionvaluepairs = [[valpair.split(':', 1) for valpair in jobconf[opt]] for opt in pairoptions]

    nonpairoptions = [opt for opt in simoptions if not opt in pairoptions]
    nonpairvalues = [jobconf[opt] for opt in nonpairoptions]

    # create jobs for all possible combinations of mulitvaloption values...
    multivaloptions = [opt for opt in nonpairoptions if len(jobconf[opt]) > 1]
    
    #allsimvalues = [jobconf[opt] for opt in simoptions]

    # do not allow paired options that also occur as single/multival options:
    for (opt1, opt2) in optionpairs:
        for opt in (opt1, opt2):
            if opt in nonpairoptions:
                raise Exception("Error: %s occurs both as paired and as non-paired option" % opt)

    nonpairvalues_combinations = [values for values in itertools.product(*nonpairvalues)]

    # list of option->value dictionaries
    optvals = []
    for values in nonpairvalues_combinations:
        optval_shell = OrderedDict([(opt, val) for (opt, val) in zip(nonpairoptions, values)])
        if pairoptions == []:
            optvals.append(optval_shell)
        else:
            ## add all combinations pair-option values
            
            optvals_to_be_added = [optval_shell]

            for (opt1, opt2), valpairs in zip(optionpairs, pairoptionvaluepairs):
                unfinished_optvals = optvals_to_be_added
                optvals_to_be_added = []

                for u_o in unfinished_optvals:
                    for (val1, val2) in valpairs:
                        optval = u_o.copy()
                        optval[opt1] = val1
                        optval[opt2] = val2
                        optvals_to_be_added.append(optval)

            optvals.extend(optvals_to_be_added)
 
    # multivaloptions is used to create the job identifier. To make things more consistent,
    # from now on also include the pairoptions in multivaloptions:
    multivaloptions += [opt for optpair in optionpairs for opt in optpair]
    # also include those options explicitly set to be used for the jobname
    multivaloptions += [opt for opt in optionsForJobName]
    # Make unique!
    multivaloptions = list(set(multivaloptions))

    # iterate over all possible combinations of simulation values.
    for optval in optvals:
        commandline = executable + " " + " ".join(["--%s %s" % (opt, val) for opt,val in optval.items()])
        subdir = jobprefix + "_".join(["%s%s" % (opt, optval[opt]) for opt in multivaloptions])
        print subdir,

        if subdir in runningJobs:
            print "already running, skipping"
            continue

        print
        print commandline

        job = jobTemplate
        for jobopt in joboptions:
            job = job.replace("%" + jobopt, jobconf[jobopt][0])
        job = job.replace("%constraint", constraint)
        job = job.replace("%cores", cores)
        job = job.replace("%maindir", maindir)
        job = job.replace("%subdir", subdir)
        job = job.replace("%commandline", commandline)
        job = job.replace("%mintime", mintime)        
#       job = job.replace("%wtimeseconds", walltimeInSeconds(jobconf["walltime"][0]))

        jobname = subdir
        jobfilename = "%s/job.sh" % subdir
        outputfile = "%s/output.%%j.log" % subdir # %j is replaced by the job allocation number
        errorfile  = "%s/error.%%j.log" % subdir
        qsubcommand = "sbatch %s --job-name=%s --output=%s --error=%s %s" % (
            partitionOption, jobname, outputfile, errorfile, jobfilename)

        # put job script into a file
        mkdir_p(subdir)
        with open(jobfilename, 'w') as jobfile:
            jobfile.write(job)

        # submit job
        try:
            stdout_and_stderr = check_output(qsubcommand, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, e:
            print "Command:", e.cmd
            print "Error:", e.output
            raise

# return a dictionary {jobid -> jobname}.
# This is just so much easier than with Torque..
def queryRunningJobs(username):
    results = {}
    # TODO: does this need special case handling for array jobs?
    formatstring = r'"%.8i %j"' # jobid followed by jobname
    commandline = 'squeue --format=%s -u %s' % (formatstring, username)
    try:
        squeueOutput = check_output(commandline, shell=True)[1:] # discard header

        for line in squeueOutput.splitlines():
            id, name = line.strip().split()
            results[id] = name
    except subprocess.CalledProcessError as e:
        print e
        print 
    
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #important argument:
    parser.add_argument("--executable", nargs=1, type=str, help="specify simulation commandline without parameters")
    #positional arguments:
    parser.add_argument("fnames", nargs='*', help="job filenames to process (default: simulation.job)")
    args = parser.parse_args()

    executable = args.executable[0]  # list with only 1 entry

    if args.fnames:
        fnames = args.fnames
    else:
        fnames = ['simulation.job']

    for fn in fnames:
        if not os.path.exists(fn):
            raise IOError("No such file: " + fn)
        
    runningJobs = queryRunningJobs("$USER").values()

#    try:
    for fn in fnames:
        process(fn, runningJobs, executable=executable)
    # except Exception as e:
    #     print e
    #     print
    #     parser.print_help()
