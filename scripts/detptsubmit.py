#!/usr/bin/env python
from collections import OrderedDict
import itertools
import os
import os.path
import re
import subprocess
import argparse
import random
from slurm_seconds_left import walltimeInSeconds, secondsToTimeString
from subprocess import check_output
import numpy as np
from scripthelpers import mkdir_p, getHeader, parseConf

# SLURM adapted version of the various old (re)submit scripts, for the Intel compiled code

# only re-submit jobs which are no longer running and have not reached the number of targetSweeps;

# parse qstat -u $USER to find currently running jobs
# check the name of each with qstat -f
# if a name matches one of the subdirectories, skip that subdirectory
# for the remaining subdirectories:
#  look inside for info.dat, if it exists and sweepsDone < targetSweeps, skip the directory
#  else: resubmit that job!

# New: continuously resubmit the job from the job file to continue work
#  - put job id in error / output file name
#  - need to put job script into a file afterall to make this work

# Newer: support setting options via command line:
#  --sweeps          # new target sweeps
#  --saveInterval    # new interval in sweeps between saves
#  --executable      # simulation executable to run

# This version detptsubmit.py is adapted for the replica exchange /
# parallel tempering variation: detqmcptsdw
# This uses MPI for each simulation run.

# 2015-02-05: also support setting the executable via the job file,
# our commandline argument will overwrite that
#
# 2016-02-01: --excludeJobs support

# template for jobs, strings preceded by % to be replaced before submitting
jobTemplate_cheops = r"""#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --nodes=%nodes
#SBATCH --ntasks=%ntasks
#SBATCH --constraint=%constraint
#SBATCH --mem-per-cpu=%mem
#SBATCH --time=%walltime
#SBATCH --time-min=%mintime
#SBATCH --account=AG-Trebst
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mgerlac2@uni-koeln.de

module add %modules

# circumvent potential DAPL bug (Intel MPI)
export I_MPI_FABRICS="shm:ofa"

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

srun --msg-timeout=60 -K1 -n $SLURM_NTASKS %commandline
RETVAL=$?

TIME2=$(date +%s)

# output finish time
echo finish time
date
echo

# only resubmit if exit status 0 and enough time passed (e.g. simulation has not finished after 1 hour)
# an if we do not find the file "ABORT.%j" here or in the directory above

if [[ ("$TIME2"-"$TIME1" -gt 3600) && ( "$RETVAL" -eq 0 ) && \
      ( ! -e "ABORT.$SLURM_JOBID" ) && ( ! -e "../ABORT.$SLURM_JOBID" ) && \
      ( ! -e "ABORT.all" ) && ( ! -e "../ABORT.all" ) ]]; then
    # resubmit the same job
    echo resubmit
    %qresubcommand
    echo
fi



# submit a rsync job to save to /projects.
# But only do that if the exit state is zero => else we may have
# corrupt data and do not want to overwrite the valid data on
# /projects
if [[ ( "$RETVAL" -eq 0 ) ]]; then
   %qrsyncsubcommand
fi

exit "$RETVAL"
"""

jobTemplate_jureca = r"""#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=%nodes
#SBATCH --ntasks=%ntasks
#SBATCH --ntasks-per-node=24
#SBATCH --time=%walltime
#SBATCH --time-min=%mintime
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gerlach@thp.uni-koeln.de

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

srun -n $SLURM_NTASKS %commandline
RETVAL=$?

TIME2=$(date +%s)

# output finish time
echo finish time
date
echo

# only resubmit if exit status 0 and enough time passed (e.g. simulation has not finished after 1 hour)
# an if we do not find the file "ABORT.%j" here or in the directory above

if [[ ("$TIME2"-"$TIME1" -gt 3600) && ( "$RETVAL" -eq 0 ) && \
      ( ! -e "ABORT.$SLURM_JOBID" ) && ( ! -e "../ABORT.$SLURM_JOBID" ) && \
      ( ! -e "ABORT.all" ) && ( ! -e "../ABORT.all" ) ]]; then
    # resubmit the same job
    echo resubmit
    %qresubcommand
    echo
fi

exit "$RETVAL"
"""


rsyncJobTemplate = r"""#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=128mb
#SBATCH --time=2:00:00
#SBATCH --account=AG-Trebst
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mgerlac2@uni-koeln.de

# make sure target directory in /projects exists
MAINDIR=%maindir
SUBDIR=%subdir
BASENAME_MAINDIR=$(basename $MAINDIR)
TARGET_DIR=/projects/ag-trebst/mgerlac2/det/$BASENAME_MAINDIR/$SUBDIR
mkdir -p $TARGET_DIR

# rsync the subdir to /projects
cd $MAINDIR
rsync -arP $SUBDIR/ $TARGET_DIR
"""


# these do not allow multiple values
possibleJobOptions = ['jobprefix', 'onlynew', 'mem', 'walltime', 'mintime',
                      'maindir', 'constraint', 'nodes',
                      'ntasks', 'partition', 'executable', 'modules']

def get_build_host():
    fzj_path = "/etc/FZJ/systemname"
    if os.path.exists(fzj_path):
        with open(fzj_path, 'r') as f:
            build_host = f.read().strip()
    else:
        hostname_output = subprocess.check_output("hostname -f", shell=True)
        if "cheops" in hostname_output:
            build_host = "cheops"
    return build_host




# construct jobs from jobconffilename, only submit them if no job of
# that name is in list runningJobs or excludeJobs
#
# ignore onlynew
def process(jobconffname, runningJobs, executable=None,
            newSaveInterval=0, newSweeps=0, excludeJobs=None):
    if excludeJobs is None:
        excludeJobs = []
    
    build_host = get_build_host()
    
    jobconf = parseConf(open(jobconffname, 'r'))

    options = jobconf.keys()

    # options to be passed to simulation instances
    simoptions = [opt for opt in options if opt not in possibleJobOptions]
    # the rest, to be replaced in the jobTemplate
    joboptions = [opt for opt in options if opt in possibleJobOptions]

    # if rValues is not specified explicitly in the job file, check if
    # paramters rMin, rMax and rCount are given. Then we can compute
    # rValues here as a linspace.
    rValues_specified = False
    for opt in simoptions:
        if opt == "rValues":
            rValues_specified = True
            break
        else:
            # possibility of having rValues in an option pair
            m = re.match(r'(.+):(.+)', opt)
            if m:
                if m.group(1) == "rValues" or m.group(2) == "rValues":
                    rValues_specified = True
                    break
    if not rValues_specified:
        try:
            jobconf["rValues"] = []
            for rMin_str, rMax_str, rCount_str in zip(jobconf["rMin"], jobconf["rMax"], jobconf["rCount"]):
                rMin = float(rMin_str)
                rMax = float(rMax_str)
                rCount = int(rCount_str)
                rValues = np.linspace(rMin, rMax, rCount)
                jobconf["rValues"].append( " ".join(["%.3f" % f for f in rValues]) )
            simoptions.append("rValues")
            simoptions.remove("rMin")
            simoptions.remove("rMax")
            simoptions.remove("rCount")
            del jobconf["rMin"]
            del jobconf["rMax"]
            del jobconf["rCount"]
        except KeyError:
            raise Exception("No rValues specified")
            
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
    # if 'cores' in joboptions:
    #     cores = jobconf['cores'][0]
    # else:
    #     cores = "1"
    if 'nodes' in joboptions:
        nodes = jobconf['nodes'][0]
    else:
        nodes = "1"
    if 'ntasks' in joboptions:
        ntasks = jobconf['ntasks'][0]
    else:
        ntasks = "1"
    if 'constraint' in joboptions:
        constraint = jobconf['constraint'][0]
    else:
        constraint = '""'
    if 'partition' in joboptions and jobconf['partition'][0] == 'devel':
        partitionOption = "--partition=devel"
    else:
        partitionOption = ""
    if 'modules' in joboptions:
        modules = jobconf['modules'][0]
    else:
        modules = "intel/15.0 mkl/11.2 intelmpi/5.0.3"

    if not executable:
        if 'executable' in joboptions:
            executable = jobconf['executable'][0]
        else:
            raise Exception("No simulation executable specified!")

    if build_host == "jureca":
        if int(ntasks) % 24 != 0:
            raise Exception("Jureca: ntasks must be a multple of 24")
        if int(ntasks) != int(nodes) * 24:
            raise Exception("Jureca: ntasks must be 24 * nodes")

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

    nonpairvalues_combinations = [vals for vals in itertools.product(*nonpairvalues)]

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
        # potentially update save interval and sweeps
        if newSaveInterval != 0 and newSaveInterval != int(optval['saveInterval']):
            optval['saveInterval'] = str(newSaveInterval)
        if newSweeps > int(optval['sweeps']):
            optval['sweeps'] = str(newSweeps)

        commandline = executable + " " + " ".join(["--%s %s" % (opt, val) for opt,val in optval.items()])
        subdir = jobprefix + "_".join(["%s%s" % (opt, optval[opt]) for opt in multivaloptions])
        print subdir,

        if subdir in runningJobs:
            print "already running, skipping"
            continue
        if subdir in excludeJobs:
            print "excluding this jobname, skipping"
            continue

        subdirInfo = subdir + "/info.dat"
        if os.path.exists(subdirInfo):
            info = parseConf(getHeader(subdirInfo))
            sweeps = int(info["sweeps"][0]) # use newSweeps if set
            if newSweeps > sweeps:
                sweeps = newSweeps
            sweepsDone = int(info["sweepsDone"][0])
            thermalization = int(info["thermalization"][0])
            sweepsDoneThermalization = int(info["sweepsDoneThermalization"][0])
            if thermalization == sweepsDoneThermalization and sweeps == sweepsDone:
                print "job finished"
                continue
            else:
                print "will be continued",

        print

        print commandline

        if build_host == "cheops":
            # Add some variation to the requested wall times (add up too 30
            # extra minutes)
            extra_seconds = random.randint(0, 30*60)
        else:
            extra_seconds = 0
        walltime = secondsToTimeString(int(walltimeInSeconds(jobconf['walltime'][0])) + extra_seconds)
        if 'mintime' in joboptions:
            mintime  = secondsToTimeString(int(walltimeInSeconds(jobconf['mintime'][0])) + extra_seconds)
        else:
            mintime = walltime

        if build_host == "cheops":
            jobTemplate = jobTemplate_cheops
        elif build_host == "jureca":
            jobTemplate = jobTemplate_jureca

        job = jobTemplate
        job = job.replace("%mintime", mintime)
        job = job.replace("%walltime", walltime)
        for jobopt in joboptions:
            job = job.replace("%" + jobopt, jobconf[jobopt][0])
        job = job.replace("%constraint", constraint)
        # job = job.replace("%cores", cores)
        job = job.replace("%nodes", nodes)
        job = job.replace("%ntasks", ntasks)
        job = job.replace("%maindir", maindir)
        job = job.replace("%subdir", subdir)
        job = job.replace("%commandline", commandline)
        job = job.replace("%modules", modules)
        # job = job.replace("%wtimeseconds", walltimeInSeconds(jobconf["walltime"][0]))

        jobname = subdir
        jobfilename = "%s/job.sh" % subdir
        outputfile = "%s/output.%%j.log" % subdir # %j is replaced by the job allocation number
        errorfile  = "%s/error.%%j.log" % subdir
        qsubcommand = "sbatch %s --job-name=%s --output=%s --error=%s %s" % (
            partitionOption, jobname, outputfile, errorfile, jobfilename)

        # resubmitting the same job from within the job if it exited gracefully
        if build_host == "cheops":
            loginNode = "cheops0"
        elif build_host == "jureca":
            loginNode = "jrl02"
        resubmit_outputfile = "$WORKDIR/output.%j.log"       # %j is replaced by the job allocation number
        resubmit_errorfile  = "$WORKDIR/error.%j.log"
        qresubcommand = "ssh $USER@%s sbatch %s --job-name=%s --output=%s --error=%s $WORKDIR/job.sh" % (
            loginNode, partitionOption, jobname, resubmit_outputfile, resubmit_errorfile)
        job = job.replace("%qresubcommand", qresubcommand)

        if build_host == "cheops":
            # afterwards submit a job to copy to /projects.
            rsyncJob = rsyncJobTemplate.replace('%maindir', maindir).replace("%subdir", subdir)
            rsyncJobname = "rsync-" + jobname
            rsyncJobfilename = "%s/rsyncjob.sh" % subdir
            rsync_outputfile = "$WORKDIR/rsync-output.%j.log"       # %j is replaced by the job allocation number
            rsync_errorfile  = "$WORKDIR/rsync-error.%j.log"
            mkdir_p(subdir)
            with open(rsyncJobfilename, 'w') as rsyncJobfile:
                rsyncJobfile.write(rsyncJob)
            qrsyncsubcommand = "ssh $USER@cheops0 sbatch --job-name=%s --output=%s --error=%s $WORKDIR/rsyncjob.sh" % (
                rsyncJobname, rsync_outputfile, rsync_errorfile)
            job = job.replace("%qrsyncsubcommand", qrsyncsubcommand)        

        # put master job script into a file
        mkdir_p(subdir)
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


# return a dictionary {jobid -> jobname}.
# This is just so much easier than with Torque..
def queryRunningJobs(username, ignore_states=None):
    """return a dictionary {jobid -> jobname} for jobs currently in the queue
    
    If ignore_states is a list of job state shorthands like CG, then
    do not return jobs in these states even if they are in the queue.
    """
    if ignore_states is None:
        ignore_states = []
    for idx, s in enumerate(ignore_states):
        ignore_states[idx] = s.upper()

    results = {}
    # TODO: does this need special case handling for array jobs?
    formatstring = r'"%.8i %j %t"' # jobid, jobname, jobstate (short form)
    commandline = 'squeue --format=%s -u %s' % (formatstring, username)
    try:
        squeueOutput = check_output(commandline, shell=True)[1:] # discard header

        for line in squeueOutput.splitlines():
            jobid, jobname, jobstate = line.strip().split()
            if jobstate.upper() in ignore_states:
                continue
            else:
                results[jobid] = jobname
    except subprocess.CalledProcessError as e:
        print e
        print 
    
    return results


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # optional arguments:
    parser.add_argument("--sweeps", type=int, help="set new target sweeps")
    parser.add_argument("--saveInterval", type=int, help="set new saveInterval")
    parser.add_argument("--executable", type=str, help="specify simulation executable (otherwise take definition in job file)")
    parser.add_argument("--replaceCompleting", dest='replaceCompleting', action='store_true', help="resubmit jobs in SLURM state CG (completing), useful if they hang there for a while (off by default)")
    parser.add_argument("--no-replaceCompleting", dest='replaceCompleting', action='store_false', help="do not resubmit jobs in SLURM state CG (completing) -- this is the default")
    parser.set_defaults(replaceCompleting=False)
    # positional arguments:
    parser.add_argument("fnames", nargs='*', help="job filenames to process (default: simulation.job)")
    # extra optional argument
    parser.add_argument("--excludeJobs", nargs='*', help="list of job-names to exclude from (re-)submitting (e.g. because they have irreparable failures).  Add this after the list of job files.")
    args = parser.parse_args()

    if args.sweeps:
        newSweeps = args.sweeps
    else:
        newSweeps = 0
    if args.saveInterval:
        newSaveInterval = args.saveInterval
    else:
        newSaveInterval = 0
    if args.executable:
        executable = args.executable
    else:
        executable = None       # take default value from simulation.job
    if args.fnames:
        fnames = args.fnames
    else:
        fnames = ['simulation.job']

    for fn in fnames:
        if not os.path.exists(fn):
            raise IOError("No such file: " + fn)

    if args.excludeJobs:
        excludeJobs = args.excludeJobs
    else:
        excludeJobs = []

    replaceCompleting = args.replaceCompleting
    running_jobs_ignore_states = ["CG"] if replaceCompleting else []
    
    runningJobs = queryRunningJobs("$USER", running_jobs_ignore_states).values()

#    try:
    for fn in fnames:
        process(fn, runningJobs, executable=executable,
                newSaveInterval=newSaveInterval, newSweeps=newSweeps,
                excludeJobs=excludeJobs)
    # except Exception as e:
    #     print e
    #     print
    #     parser.print_help()
