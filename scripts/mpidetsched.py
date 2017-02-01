#!/usr/bin/env python

# This is the script to be run via the SLURM batch system 

from collections import OrderedDict
import itertools
import os
import os.path
import errno
import re
import subprocess
import argparse
from mpi4py import MPI
import time
from scripthelpers import mkdir_p, getHeader, parseConf

possibleJobOptions = ['jobprefix', 'onlynew', 'cores', 'mem',
                      'walltime', 'mintime', 'maindir', 'constraint', 'nodes',
                      'ntasks', 'partition']   # these do not allow multiple values


# This is a helper for findWhatToDo
# remove the tasks that don't need to be run (because they are finished)
def filterOutDone(taskCommandlineDictionary, newSweeps=0):
    filteredDictionary = OrderedDict()
    for subdir, commandline in taskCommandlineDictionary.iteritems():
        print subdir,
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
            else:
                print "will be continued"
                filteredDictionary.update({subdir: commandline})
        else:
            print
            filteredDictionary.update({subdir: commandline})            
    return filteredDictionary

    

# setup a list of tasks to be computed (task name = subdirectory) with
# the respective command lines, return as a dictionary
def findWhatToDo(jobconffname, executable="detqmcsdw", newSaveInterval=0, newSweeps=0):
    jobconf = parseConf(open(jobconffname, 'r'))

    options = jobconf.keys()
    onlynew = False

    # options to be passed to simulation instances
    simoptions = [opt for opt in options if opt not in possibleJobOptions]

    if 'jobprefix' in options:
        jobprefix = jobconf["jobprefix"][0] + "_"
    else:
        jobprefix = ""
    
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

    # additionally: specify opt1:opt2 pairs with values val1_a:val2_a, val1_b:val2_b, ...
    pairoptions = [opt for opt in simoptions if re.match(r'.+:.+', opt)]
    optionpairs = [opt.split(':', 1) for opt in pairoptions]
    pairoptionvaluepairs = [[valpair.split(':', 1) for valpair in jobconf[opt]] for opt in pairoptions]

    nonpairoptions = [opt for opt in simoptions if not opt in pairoptions]
    nonpairvalues = [jobconf[opt] for opt in nonpairoptions]

    # create jobs for all possible combinations of mulitvaloption values...
    multivaloptions = [opt for opt in nonpairoptions if len(jobconf[opt]) > 1]
    
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

    # iterate over all possible combinations of simulation values, construct an ordered dictionary
    taskCommandlineDictionary = OrderedDict()
    for optval in optvals:
        # potentially update save interval and sweeps
        if newSaveInterval != 0 and newSaveInterval != int(optval['saveInterval']):
            optval['saveInterval'] = str(newSaveInterval)
        if newSweeps > int(optval['sweeps']):
            optval['sweeps'] = str(newSweeps)

        commandline = executable + " " + " ".join(["--%s %s" % (opt, val) for opt,val in optval.items()])
        subdir = jobprefix + "_".join(["%s%s" % (opt, optval[opt]) for opt in multivaloptions])

        taskCommandlineDictionary[subdir] = commandline

    return filterOutDone(taskCommandlineDictionary, newSweeps)


def walltimeInSeconds(walltimeString):
    ftr = [3600,60,1]
    return str(sum([a*b for a,b in zip(ftr, map(int,walltimeString.split(':')))]))

# Allow atleast 35 minutes of remaining runtime
safety_walltime = 35*60
#safety_walltime = 10

# Constants
TAG_SUBDIR = 1
TAG_COMMANDLINE = 2
TAG_WALLTIME_SECONDS = 3
TAG_STATUS_MSG = 4

MSG_MANAGER_HAS_WORK = 1
MSG_MANAGER_HAS_NO_WORK = 2
MSG_WORKER_FINISHED = 3
MSG_WORKER_ERROR = 4


# manage a list of task, send tasks to workers via MPI
def manager(taskCommandlineDictionary):
    still_todo = taskCommandlineDictionary
    num_tasks = len(still_todo)
    
    comm = MPI.COMM_WORLD
    jobid = os.getenv("SLURM_JOBID", default="nojobid")
    num_workers = comm.Get_size() - 1
    print "Job ID:", jobid
    print "Processes: 1 master and %d workers" % num_workers
    print "Tasks to do:", num_tasks
    active_workers = range(1, num_workers+1)
    worker_to_subdir = {}

    granted_wtime = int(os.getenv("PBS_WALLTIME", default=4294967295)) # PBS_WALLTIME is to be set 
    start_time = MPI.Wtime()    # current time in seconds
    def elapsed_time():
        return MPI.Wtime() - start_time
    def remaining_time():
        return int(granted_wtime - elapsed_time())
    print "Remaining wall time:", remaining_time()
    
    def send_task_from_still_todo(destination):
        comm.send(MSG_MANAGER_HAS_WORK, tag=TAG_STATUS_MSG, dest=destination)
        subdir, commandline = still_todo.popitem(last=False)
        comm.send(subdir, tag=TAG_SUBDIR, dest=destination)
        comm.send(commandline, tag=TAG_COMMANDLINE, dest=destination)
        comm.send(remaining_time(), tag=TAG_WALLTIME_SECONDS, dest=destination)
        worker_to_subdir[destination] = subdir
        print "Worker %d processes task %s" % (destination, subdir)

    # distribute initial work
    for i in range(1, min(num_tasks+1, num_workers+1)):
        send_task_from_still_todo(destination=i)
        
    # if necessary, tell some workers that there is no work for them
    # and dismiss them
    for i in range(num_tasks+1, num_workers+1):
        comm.send(MSG_MANAGER_HAS_NO_WORK, tag=TAG_STATUS_MSG, dest=i)
        active_workers.remove(i)

    abort_job = False
    abort_filename = "ABORT." + jobid

    failed_tasks = 0
    
    def give_new_task_if_we_have_any(destination):
        "if our list of tasks is not empty and the flag abort_job is not set, send out a new task"
        if len(still_todo) > 0 and not abort_job:
            # send new work
            send_task_from_still_todo(destination=i)
        else:
            # no more work or time is up, dismiss the worker
            comm.send(MSG_MANAGER_HAS_NO_WORK, tag=TAG_STATUS_MSG, dest=i)
            active_workers.remove(i)

    while len(active_workers) > 0:
        if not abort_job and remaining_time() < safety_walltime:
            print "Remaining wall time is less than %d seconds, stopping distribution of tasks\n" % safety_walltime
            abort_job = True
        elif not abort_job and os.path.exists(abort_filename):
            print "Found file %s, stopping distribution of tasks\n" % abort_filename
            abort_job = True
        
        # Check for idle workers
        for i in active_workers:
            if comm.Iprobe(source=i, tag=TAG_STATUS_MSG):
                msg = comm.recv(source=i, tag=TAG_STATUS_MSG)
                if msg == MSG_WORKER_FINISHED:
                    print "SUCCESS: Worker %d, task %s" % (i, worker_to_subdir[i])
                    give_new_task_if_we_have_any(destination=i)
                elif msg == MSG_WORKER_ERROR:
                    print "FAILURE: Worker %d, task %s" % (i, worker_to_subdir[i])
                    failed_tasks += 1
                    give_new_task_if_we_have_any(destination=i)                    
                else:
                    raise Exception("Manager reveived an invalid status message from worker %d" % i)
                print "Remaining wall time:", remaining_time(), ", remaining tasks:", len(still_todo)
                print
        # save some cpu
        time.sleep(0.1)
    print "End of job.\nRemaining wall time: %d\nTasks not started: %d" % (remaining_time(), len(still_todo))
    print "Failed tasks: %d" % failed_tasks
    return

# receive tasks via MPI and complete them
def worker():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    jobid = os.getenv("SLURM_JOBID", default="nojobid")

    # Receive work, as long as there is any.
    # When finished: ask manager for more work
    while True:
        msg = comm.recv(source=0, tag=TAG_STATUS_MSG)
        if msg == MSG_MANAGER_HAS_WORK:
            subdir = comm.recv(source=0, tag=TAG_SUBDIR)
            commandline = comm.recv(source=0, tag=TAG_COMMANDLINE)
            granted_wtime = comm.recv(source=0, tag=TAG_WALLTIME_SECONDS)

            # setup actual simulation run:
            mkdir_p(subdir)
            os.chdir(subdir)
            if jobid == "nojobid":
                stamp = "t"+str(int(time.time()))
            else:
                stamp = jobid
            fname_outlog = "output.%s.log" % stamp
            fname_errlog = "error.%s.log" % stamp
            if os.path.exists(fname_outlog):
                raise Exception("Output log file %s/%s already exists" % (subdir, fname_outlog))
            if os.path.exists(fname_errlog):
                raise Exception("Error log file %s/%s already exists" % (subdir, fname_errlog))
            with open(fname_outlog, 'w') as outlog, open(fname_errlog, 'w') as errlog:
                outlog.writelines(["Host name: " + MPI.Get_processor_name() + "\n"])
                outlog.writelines(["Command line:\n", commandline + "\n"])
                outlog.flush()
                # make sure the simulation process knows how much walltime is left
                sub_env = os.environ.copy()
                sub_env["PBS_WALLTIME"] = str(granted_wtime)
                sub_env["OMP_NUM_THREADS"] = str(1)
                # run the command
                returncode = subprocess.call(commandline, shell=True, stdout=outlog, stderr=errlog,
                                             env=sub_env)
                # report back to the manager
                if returncode == 0:
                    comm.send(MSG_WORKER_FINISHED, tag=TAG_STATUS_MSG, dest=0)
                else:
                    # some kind of failure, detaisl should be in the error log
                    comm.send(MSG_WORKER_ERROR, tag=TAG_STATUS_MSG, dest=0)
            os.chdir("..")
            
        elif msg == MSG_MANAGER_HAS_NO_WORK:
            # nothing to do, end
            return
        else:
            raise Exception("Worker %d received an invalid status message" % rank)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    #optional arguments:
    parser.add_argument("--sweeps", type=int, help="set new target sweeps")
    parser.add_argument("--saveInterval", type=int, help="set new saveInterval")
    parser.add_argument("--executable", type=str, help="specify simulation executable (default: detqmcsdw)")
    #positional arguments:
    parser.add_argument("fname", nargs='?', help="job filename")
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
    if args.fname:
        fname = args.fname
    else:
        fname = "simulation.job"

    # MPI action
    comm = MPI.COMM_WORLD
    if comm.Get_rank() == 0:
        taskCommandlineDictionary = findWhatToDo(fname, executable, newSaveInterval, newSweeps)
        manager(taskCommandlineDictionary)
    else:
        worker()
