#!/usr/bin/env python

# print a summary of the jobs corresponding to the subdirectories of
# the current directory

import os, os.path
import time
from collections import namedtuple
from slurm_seconds_left import secondsToTimeStringWithDays
from scripthelpers import getHeader, parseConf, queryJobsInState


def queryRunningJobs(username):
    return queryJobsInState(username, 'R')

def queryPendingJobs(username):
    return queryJobsInState(username, 'PD')    

def getSubdirs():
    return [f for f in os.listdir('.') if os.path.isdir(f)]

def file_age_in_seconds(pathname):
    return time.time() - os.path.getmtime(pathname)


if __name__ == '__main__':
    subdir_jobs = getSubdirs()
    
    subdirs_with_info = [s for s in subdir_jobs if os.path.exists(s + "/info.dat")]

    interesting_fields = ['sweeps', 'sweepsDone', 'thermalization', 'sweepsDoneThermalization', 'totalWallTimeSecs']

    Jobinfo = namedtuple('Jobinfo', interesting_fields)
    subdir_jobinfo = {}
    subdir_replica_count = {}
    finished_jobs = []
    incomplete_jobs = []
    for s in subdirs_with_info:
        info = parseConf(getHeader(s + "/info.dat"))

        try:
            jinfo = Jobinfo(*[int(info[f][0]) for f in interesting_fields])
        except KeyError:
            # info without entries like sweepsDone, don't show
            # (probably is simindexjoined)
            continue

        subdir_jobinfo[s] = jinfo

        if jinfo.sweepsDone == jinfo.sweeps and jinfo.thermalization == jinfo.sweepsDoneThermalization:
            finished_jobs.append(s)
        else:
            incomplete_jobs.append(s)

        subdir_replica_count[s] = 1
        if "controlParameterValues" in info:
            subdir_replica_count[s] = len(info["controlParameterValues"][0].split())
        
    
    running_jobs = queryRunningJobs('$USER').values()

    dead_never_started_jobs = [s for s in subdir_jobs if s not in subdirs_with_info
                               and s not in running_jobs]
    dead_incomplete_jobs = [s for s in incomplete_jobs if s not in running_jobs]
    running_incomplete_jobs = [s for s in running_jobs if s in incomplete_jobs]
    running_nodata_jobs = [s for s in running_jobs if s in subdir_jobs and s not in subdirs_with_info]

    titles = ['Name', 'Replicas',  'Thermalization    ', 'Sweeps          ', 'totalWallTime    ', 'infoAge' ]
    def printTitles(fieldwidth1):
        print "{0:<{1}}".format(titles[0], fieldwidth1), 
        for x in titles[1:]:
            print x,
        print
    def printJobList(category, joblist):
        if len(joblist) == 0:
            return
        print category
        fieldwidth = max( [len(j) for j in joblist] )        
        printTitles(fieldwidth)
        for j in joblist:
            # Name
            print "{0:<{1}}".format(j, fieldwidth),
            try:
                info = subdir_jobinfo[j]

                # Number of replicas
                print "{0}".format(subdir_replica_count[j]).ljust(len(titles[1])),
            
                # Thermalization, Sweeps
                for field1, field2, title in [("sweepsDoneThermalization", "thermalization", titles[2]),
                                              ("sweepsDone", "sweeps", titles[3])]:
                    print "{0} / {1}".format(getattr(info, field1), getattr(info, field2)).ljust(len(title)),
                # walltime
                print secondsToTimeStringWithDays(info.totalWallTimeSecs).ljust(len(titles[4])),
                # age of info.dat
                print secondsToTimeStringWithDays(file_age_in_seconds(os.path.join(j, "info.dat"))).ljust(len(titles[5])),

                # if j in subdir_jobinfo:
                #     for field in titles[1:]:
                #         print "{0:>{1}}".format(getattr(subdir_jobinfo[j], field), len(field)),
            except KeyError:
                pass
            print
        print
            
    printJobList("DEAD, NEVER STARTED", dead_never_started_jobs)
    printJobList("DEAD, INCOMPLETE", dead_incomplete_jobs)
    printJobList("FINISHED", finished_jobs)    
    printJobList("RUNNING, NO DATA", running_nodata_jobs)
    printJobList("RUNNING, INCOMPLETE", running_incomplete_jobs)

