import os
import sys

#Cluster Settings
NODES_MIN = 1
NODES_MAX = 4
PARTITION = 'normal'
DAYS = 3
HOURS = 0
MEMORY = 500 # in MB

def make_sbatch(jobline, fname, fname_out, jobname, nNodes_min, nNodes_max, partition, days, hours, nodeMemory): #Make batch file to execute job
    sbatch = "#SBATCH"
    fout = open(fname, 'w+')
    fout.write("#!/bin/sh\n")

    minutes = 0
    seconds = 0

    fout.write(sbatch + " --job-name=" + jobname +"\n")
    fout.write(sbatch + " --partition=" + partition + "\n")
    fout.write(sbatch + " --time=" +str(days) + "-" + str(hours) + ":" + str(minutes) + ":" + str(seconds) + " # days-hours:minutes:seconds\n")
    if nNodes_min == nNodes_max:
        fout.write(sbatch + " --nodes=" + str(nNodes_min) + "\n")
    else:
        fout.write(sbatch + " --nodes=" + str(nNodes_min) + "-" + str(nNodes_max) + "\n")
    fout.write(sbatch + " --mem-per-cpu=" + str(nodeMemory) + " # in MB\n")
    fout.write(sbatch + " -o " + str(fname_out) + "\n")
    fout.write(jobline)

    makeprogram = "chmod u+x " + fname
    os.system(makeprogram)
