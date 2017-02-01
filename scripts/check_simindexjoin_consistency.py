#!/usr/bin/env python
import subprocess
import os
import sys
from glob import glob

def check_time_series_lengths(directory):
    """Check if all the time series files in @directory contain the same
    number of lines.  This is a consistency check."""
    # $ wc -l associatedEnergy.series 
    # 23227 associatedEnergy.series
    counts = []
    for f in glob(os.path.join(directory, "*.series")):
        output = subprocess.check_output(['wc', '-l', f])
        count = int(output.split()[0])
        counts.append(count)
    # check if all elements in counts are equal (or the list is empty)
    return len(set(counts)) <= 1

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print "Usage: check_simindexjoin_consistency.py directory"
        print "Will print directory, if *.series files in directory have differing lengths"
        exit()

    directory = sys.argv[1]

    if not check_time_series_lengths(directory):
        print directory
        
