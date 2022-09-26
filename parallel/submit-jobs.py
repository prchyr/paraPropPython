from os import system, listdir
from os.path import isfile, join
from sys import argv

"""
This program submits a list of jobs to a cluster (i.e. pleiades)
Each job is defined using a bash file stored in a directory
-> this directory and the contained bash scripts are defined beforehand by running pleiades-job-sbatch.py

"""

path2jobs = argv[1]

def submitjob(path, jobfile):
    sbatch = "sbatch"
    command = sbatch + " " + path + "/" + jobfile
    system(command)

if __name__ == "__main__":
    shfiles = [f for f in listdir(path2jobs) if isfile(join(path2jobs, f))]
    for f in shfiles:
        submitjob(path2jobs, f)
