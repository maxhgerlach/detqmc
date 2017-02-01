#!/usr/bin/env python
from collections import defaultdict
import os
import os.path
import re
import subprocess
from subprocess import check_output
from scripthelpers import mkdir_p, getHeader, parseConf


if __name__ == "__main__":
    # used for sdw-o2-lambda3 directly on a Jureca log-in node

    # need to deal with multiple simindexes; results put onto simindexjoined
    evalOnlyNew = True
    options = "-j 20"
    sweeps_discard = 1500
    source_toplevel = os.path.abspath(".")
    dest_toplevel   = source_toplevel

    # find source directories, where we keep our stored system configurations    
    source_directories = [d for d in os.listdir(source_toplevel) if
                          os.path.isdir(os.path.join(source_toplevel, d)) and
                          'simindex' in d and
                          not 'simindexjoined' in d]
    source_directories.sort()

    ## group those source directories that only differ by their simindex? entry
    groups = defaultdict(list)
    target_directories = {}
    for d in source_directories:
        # key: remove simindex..[_] part from d
        key = re.sub(r'simindex\d+_?', '', d)
        groups[key].append(d)
        if key not in target_directories:
            target_directories[key] = re.sub(r'simindex\d+', 'simindexjoined', d)
    
    ## for each group find the p*_r* data directories; also read out
    ## the property "saveConfigurationStreamInterval", which we will
    ## use to adjust sweeps_discard
    data_directories = {}
    saveConfigurationStreamInterval = {}
    for key in groups.keys():
        example_source_dir = os.path.join(source_toplevel, groups[key][0])
        data_directories[key] = [d for d in os.listdir(example_source_dir) if
                                 d.startswith("p") and
                                 os.path.isdir(os.path.join(example_source_dir, d))]
        try:
            meta = parseConf(getHeader(os.path.join(example_source_dir, "info.dat")))
            # scsi = int(meta['saveConfigurationStreamInterval'])
            scsi = int(meta['saveConfigurationStreamInterval'][0])
        except KeyError:
            scsi = 1
        saveConfigurationStreamInterval[key] = scsi

    ## handle key after key
    for key in groups.keys():
        print key,
        print "==> consider running sdwcorr"

        input_directories  = [os.path.join(source_toplevel, d) for d in groups[key]]
        output_directory = os.path.join(dest_toplevel, target_directories[key])
        
        my_discard = sweeps_discard / saveConfigurationStreamInterval[key]
        my_options = options + " -d %d" % my_discard

        commandlines = ""
        for dd in data_directories[key]:
            ind = " ".join([os.path.join(d, dd) for d in input_directories])
            outd = os.path.join(output_directory, dd)
            if evalOnlyNew and os.path.exists(os.path.join(outd, "corr_ft.npz")):
                print "already evaluated", dd, " .. skip"
                continue
            commandlines = commandlines + ("sdwcorr --input %s --output %s %s \n" % (ind, outd, my_options))

        if commandlines == "":
            print "==> nothing to do"
            print
            continue

        # just run command after command and print its (little) output directly to stdout
        
        mkdir_p(output_directory)
        os.chdir(output_directory)   # this would be %maindir if submitted as a job

        for cl in commandlines.splitlines():
            try:
                stdout_and_stderr = check_output(cl, shell=True, stderr=subprocess.STDOUT)
                print stdout_and_stderr
            except subprocess.CalledProcessError, e:
                print "Command:", e.cmd
                print "Error:", e.output
                raise

        print
        

