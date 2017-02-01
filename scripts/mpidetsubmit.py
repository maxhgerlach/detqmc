#!/usr/bin/env python
import os
import os.path
import subprocess
import argparse
import random
from subprocess import check_output
from slurm_seconds_left import walltimeInSeconds, secondsToTimeString
from scripthelpers import parseConf

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



    
# template for jobs, strings preceded by % to be replaced before submitting
jobTemplate = r"""#!/bin/bash -l
#SBATCH --cpus-per-task=%cores
#SBATCH --nodes=%nodes
#SBATCH --ntasks=%ntasks
#SBATCH --constraint=%constraint
#SBATCH --mem-per-cpu=%mem
#SBATCH --time=%walltime
#SBATCH --time-min=%mintime
#SBATCH --account=AG-Trebst
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mgerlac2@uni-koeln.de

#module add intel/15.0 mkl/11.0 openmpi/1.6.5

# the MPI4PY-code currently only works with intelmpi
module rm openmpi
module add intel/15.0 mkl/11.0 intelmpi/5.0.2


# This catches the output from an external Python script which 
# prints the number of seconds this job has left to run. 
export PBS_WALLTIME=$(slurm_seconds_left.py)

cd %maindir

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

exit "$RETVAL"
"""


# these do not allow multiple values
possibleJobOptions = ['jobprefix', 'onlynew', 'cores', 'mem', 'walltime', 'mintime',
                      'maindir', 'constraint', 'nodes',
                      'ntasks', 'partition']


# process a jobconf file and produce the mpi job script that will run
# mpidetsched.py with paramters appropriate for the requested configuration 
# ignore onlynew
def process(jobconffname, executable="detqmcsdw", newSaveInterval=0, newSweeps=0):
    jobconf = parseConf(open(jobconffname, 'r'))

    options = jobconf.keys()

    # this script only needs to consider joboptions: the rest is
    # handled by mpidetsched.py
    joboptions = [opt for opt in options if opt in possibleJobOptions]

    # # filter out option value list of the form "option, [empty]"
    # # ...
    # optionsForJobName = []
    # for opt in simoptions:
    #     values = jobconf[opt]
    #     if values[-1] == "":
    #         # the last value is empty
    #         # -> use this opt for the jobname, remove the empty entry from jobconf
    #         optionsForJobName.append(opt)
    #         del values[-1]

    if 'maindir' in joboptions:
        maindir = jobconf["maindir"][0]
    else:
        maindir = os.getcwd()
    if 'jobprefix' in joboptions:
        jobprefix = jobconf["jobprefix"][0] + "_"
    else:
        jobprefix = ""
    # we will want to use either jobprefix or jobconffname for the
    # SLURM job name
    if 'cores' in joboptions:
        cores = jobconf['cores'][0]
    else:
        cores = "1"
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

    # Add some variation to the requested wall times (add up too 30
    # extra minutes)
    extra_seconds = random.randint(0, 30*60)
    walltime = secondsToTimeString(int(walltimeInSeconds(jobconf['walltime'][0])) + extra_seconds)
    if 'mintime' in joboptions:
        mintime  = secondsToTimeString(int(walltimeInSeconds(jobconf['mintime'][0])) + extra_seconds)
    else:
        mintime = walltime

    # TODO Check, ob es ueberhaupt noetig ist, zu resubmitten ?
    # Problem: aufwaendig, muessten schon hier alle Jobnamen erzeugen
    # und Unterverziechnisse checken

    # prepare job that will call mpidetsched.py
    commandline = "mpidetsched.py --sweeps %s --saveInterval %s --executable %s %s" % (
        newSweeps, newSaveInterval, executable, jobconffname)
    job = jobTemplate
    job = job.replace("%mintime", mintime)
    job = job.replace("%walltime", walltime)
    for jobopt in joboptions:
        job = job.replace("%" + jobopt, jobconf[jobopt][0])
    job = job.replace("%constraint", constraint)
    job = job.replace("%cores", cores)
    job = job.replace("%nodes", nodes)
    job = job.replace("%ntasks", ntasks)
    job = job.replace("%maindir", maindir)
    job = job.replace("%commandline", commandline)


    jobname = jobprefix[:-1] or jobconffname
    jobfilename = "%s/%s_job.sh" % (maindir, jobname)
    outputfile = "%s/output.%s.%%j.log" % (maindir, jobname) # %j is replaced by the job allocation number
    errorfile = "%s/error.%s.%%j.log" % (maindir, jobname)

    qsubcommand = "sbatch %s --job-name=%s --output=%s --error=%s %s" % (
        partitionOption, jobname, outputfile, errorfile, jobfilename)

    # resubmitting the same job from within the job if it exited gracefully
    resubmit_outputfile = "%s/output.%s.%%j.log" % (maindir, jobname) # %j is replaced by the job allocation number
    resubmit_errorfile = "%s/error.%s.%%j.log" % (maindir, jobname)
    qresubcommand = "ssh $USER@cheops0 sbatch %s --job-name=%s --output=%s --error=%s %s" % (
        partitionOption, jobname, resubmit_outputfile, resubmit_errorfile, jobfilename)
    job = job.replace("%qresubcommand", qresubcommand)

    # put job script into a file
    with open(jobfilename, 'w') as jobfile:
        jobfile.write(job)

    # submit job
    try:
        stdout_and_stderr = check_output(qsubcommand, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError, e:
        print "Command:", e.cmd
        print "Error:", e.output
        raise
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #optional arguments:
    parser.add_argument("--sweeps", type=int, help="set new target sweeps")
    parser.add_argument("--saveInterval", type=int, help="set new saveInterval")
    parser.add_argument("--executable", type=str, help="specify simulation executable (default: detqmcsdw)")
    #positional arguments:
    parser.add_argument("fnames", nargs='*', help="job filenames to process (default: simulation.job)")
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
        executable = "detqmcsdw"
    if args.fnames:
        fnames = args.fnames
    else:
        fnames = ['simulation.job']

    for fn in fnames:
        if not os.path.exists(fn):
            raise IOError("No such file: " + fn)

    for fn in fnames:
        process(fn, executable=executable, newSaveInterval=newSaveInterval, newSweeps=newSweeps)

