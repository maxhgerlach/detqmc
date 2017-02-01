#!/usr/bin/env python

import subprocess

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

def secondsToTimeString(total_seconds):
    "convert $total_seconds to HH:MM:SS format"
    ts = int(total_seconds)
    hours = ts / (60*60)
    ts -= hours*60*60
    minutes = ts / 60
    ts -= minutes * 60
    seconds = ts
    return "%d:%02d:%02d" % (hours, minutes, seconds)

def secondsToTimeStringWithDays(total_seconds):
    "convert $total_seconds to [Days d]HH:MM:SS format"
    ts = int(total_seconds)
    days = ts / (24*60*60)
    ts -= days*(24*60*60)
    hours = ts / (60*60)
    ts -= hours*60*60
    minutes = ts / 60
    ts -= minutes * 60
    seconds = ts
    if days == 0:
        return "%d:%02d:%02d" % (hours, minutes, seconds)
    else:
        return "%d-%02d:%02d:%02d" % (days, hours, minutes, seconds)

    
if __name__ == "__main__":
    time_str = subprocess.check_output(r"squeue -h -j $SLURM_JOBID -o %L", shell=True)
    print walltimeInSeconds(time_str)
