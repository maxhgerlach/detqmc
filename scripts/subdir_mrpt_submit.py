#!/usr/bin/env python
from glob import glob
import os
import os.path
import subprocess
from subprocess import check_output
import argparse
from scripthelpers import is_file_older_than_a_day, queryRunningJobs, is_file_more_than_a_day_older_than_file

# Run mrpt reweighting in each subdir, one job per subdir Note: One
# job takes like 30seconds or a few minutes, so this is not going to
# be to the admins' liking

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

export WORKDIR=%maindir/%subdir
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



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #important argument:
    default_options = "--cp-auto-range 0.001 -b 1000 -j 20 --non-iterative -i 100000 --no-tau"
    parser.add_argument("--options", nargs=1, type=str, default=[default_options], help="options to be passed to mrpt, default: " + default_options)
    parser.add_argument("--cores", nargs=1, type=int, default=[4], help="number of cpu cores to run each mrpt instance on, default: 4")
    parser.add_argument("--evalOnlyNew", dest='evalOnlyNew', action='store_true', help="skip already evaluated subdirectories")
    parser.add_argument("--no-evalOnlyNew", dest='evalOnlyNew', action='store_false', help="do not skip already evaluated subdirectories")
    parser.set_defaults(evalOnlyNew=False)
    parser.add_argument("--evalOnlyRecent", dest='evalOnlyRecent', action='store_true', help="skip already evaluated subdirectories, if data is older than 24 hours")
    parser.add_argument("--no-evalOnlyRecent", dest='evalOnlyRecent', action='store_false', help="do not skip already evaluated subdirectories, if data is older than 24 hours")
    parser.set_defaults(evalOnlyRecent=False)
    parser.add_argument("--evalOnlyNecessary", dest='evalOnlyNecessary', action='store_true', help="skip already evaluated subdirectories, if data is more than 24 hours than last evaluation result")
    parser.add_argument("--no-evalOnlyNecessary", dest='evalOnlyNecessary', action='store_false', help="do not skip already evaluated subdirectories, if data is more than 24 hours than last evaluation result")
    parser.set_defaults(evalOnlyNecessary=False)

    args = parser.parse_args()

    options = args.options[0]
    cores = args.cores[0]
    evalOnlyNew = args.evalOnlyNew
    evalOnlyRecent = args.evalOnlyRecent
    evalOnlyNecessary = args.evalOnlyNecessary

    # call in subdir
    commandline = "mrpt " + options + " --info ./info.dat --savez z.dat --loadz z.dat ./p*/associatedEnergy.series ./p*/normMeanPhi.series"

    # first level subdirs: contain data
    subdirs = [o for o in os.listdir(".") if os.path.isdir(o)]

    runningJobs = queryRunningJobs("$USER").values()
    
    for sd in subdirs:
        print sd,
        jobname = "mrpt_" + sd

        if jobname in runningJobs:
            print "==> already running"
            continue
        elif evalOnlyRecent and os.path.exists(sd + "/mrpt-dos.dat") and is_file_older_than_a_day(sd + "/info.dat"):
            print "=> already evaluated and data is older than 24 hours"
            continue
        elif evalOnlyNew and os.path.exists(sd + "/mrpt-dos.dat"):
            print "=> already evaluated"
            continue
        elif evalOnlyNecessary and os.path.exists(sd + "/mrpt-dos.dat") and not is_file_more_than_a_day_older_than_file(sd + "/mrpt-dos.dat", sd + "/info.dat"):
            print "=> already evaluated and data is not at least a day younger than evaluation result: skip evaluation, but take into account"
            continue
        elif not glob(sd + "/p*/*.series"):
            print "=> no time series found: skip"
            continue

        print "==> submit mrpt job"

        job = jobTemplate.replace("%subdir", sd) \
                         .replace("%maindir", os.getcwd()) \
                         .replace("%cores", str(cores)) \
                         .replace("%commandline", commandline)        

        jobfilename = "%s/mrpt_job.sh" % sd
        outputfile = "%s/mrpt_output.%%j.log" % sd # %j is replaced by the job allocation number
        errorfile  = "%s/mrpt_error.%%j.log" % sd
        qsubcommand = "sbatch --job-name=%s --output=%s --error=%s %s" % (
            jobname, outputfile, errorfile, jobfilename)

        with open(jobfilename, 'w') as jobfile:
            jobfile.write(job)
        # submit job
        try:
            stdout_and_stderr = check_output(qsubcommand, shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError, e:
            print "Command:", e.cmd
            print "Error:", e.output
            raise

