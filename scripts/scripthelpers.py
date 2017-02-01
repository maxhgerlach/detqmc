import datetime
import os
import sys
import re
import errno
import subprocess
from collections import OrderedDict

# functions used by parallelevalcollectresults.py, detptsubmit.py etc.



# return a dictionary {jobid -> jobname}.
# This is just so much easier than with Torque..
def queryRunningJobs(username):
    results = {}
    # TODO: does this need special case handling for array jobs?
    formatstring = r'"%.8i %j"' # jobid followed by jobname
    commandline = 'squeue --format=%s -u %s' % (formatstring, username)
    try:
        squeueOutput = subprocess.check_output(commandline, shell=True)[1:] # discard header

        for line in squeueOutput.splitlines():
            id, name = line.strip().split()
            results[id] = name
    except subprocess.CalledProcessError as e:
        print >>sys.stderr, e
        print >>sys.stderr 
    
    return results

# return a dictionary {jobid -> jobname}.
def queryJobsInState(username, desired_state):
    formatstring = r'"%.8i %.2t %j"' # jobid followed by state and jobname
    commandline = 'squeue --format=%s -u %s' % (formatstring, username)
    try:
        squeueOutput = subprocess.check_output(commandline, shell=True)[1:] # discard header
    except:
        squeueOutput = ""
    results = {}
    for line in squeueOutput.splitlines():
        id, state, name = line.strip().split()
        if state == desired_state:
            results[id] = name
    return results
    

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


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
        if len(header) == 0:
            print >>sys.stderr, "file is empty:", f
            return None
        else:
            return header
    except:
        print >>sys.stderr, "can't open", f
        return None
        
def parseHeader(header):
    if header == None:
        return None
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

# take a string with configuration settings confOption = confValue1, confValue2, ... \n
# return a ordered dictionary mapping confOption -> [List of confValues]
# If the list ends in a "," not followed by another confValue, this confOption 
# will be used to construct the jobname regardless of it having only one possible
# value or not.
def parseConf(header):
    d = OrderedDict()
    # Regexp to match "key = anything":
    p = re.compile(r"\s* ([^\s=]+) \s* = (.+)", re.VERBOSE)
    for line in header:
        if line[0] == "#":
            continue
        m = p.match(line)
        try:
            d[m.group(1)] = [s.strip() for s in m.group(2).split(',')]
        except:
            pass
    return d



def is_file_older_than_days(filepath, days):
    days_ago = datetime.datetime.now() - datetime.timedelta(days=days)
    filetime = datetime.datetime.fromtimestamp(os.path.getctime(filepath))
    return filetime < days_ago

def is_file_older_than_a_day(filepath):
    return is_file_older_than_days(filepath, 1)

def is_file_more_than_a_day_older_than_file(filepath1, filepath2):
    "check that file1 is sufficiently older than file2"
    filetime1 = datetime.datetime.fromtimestamp(os.path.getctime(filepath1))
    filetime2 = datetime.datetime.fromtimestamp(os.path.getctime(filepath2))
    filetimedelta = filetime2 - filetime1   # how much later in time was file2 changed than file1
    return filetimedelta > datetime.timedelta(days = 1)
    
