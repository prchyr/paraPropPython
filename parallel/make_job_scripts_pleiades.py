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
    """
        jobline :
        fname :
        fname_out :
        jobname :
        nNodes_min
        nNodes_max
        partion :
        days :
        hours :
        nodeMemory :
    """
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

#Input Arguments
fname_in = sys.argv[1] #List of jobs to be executed
path2jobs = sys.argv[2] #Directory to save bash scripts to

if __name__ == "__main__":
    fin = open(fname_in, "r+")
    if os.path.isdir(path2jobs) == False:
        os.mkdir(path2jobs)
    for jobline in fin:
        cols = jobline.split()
        src_depth = cols[3]
        jobname = "src" + src_depth
        sbatch_file = path2jobs + "/src" + src_depth + ".sh"
        out_file = path2jobs + "/src" + src_depth + ".out"
        make_sbatch(jobline, sbatch_file, out_file, jobname, NODES_MIN, NODES_MAX, PARTITION, DAYS, HOURS, MEMORY)
    fin.close()
