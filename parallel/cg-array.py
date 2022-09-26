import os
import sys
import numpy as np
import subprocess
import time

A_arr = np.arange(0.1, 1.0, 0.05)
N = len(A_arr)
nRandom = 100
fname_config = 'config_aletsch.txt'
fname_nprofile = 'share/colle-gnifetti.txt'
directory = 'colle-gnifetti'

fname_prefix = directory + '_A='
fname_out = directory + '.h5'
fname_scan = directory + '-scan.txt'
fscan = open(fname_scan, 'w')
for i in range(N):
    print(i, A_arr[i])
    fname_matrix = directory + '_matrix.h5'

    command_genSim = 'python genSim.py ' + fname_config + ' ' + fname_nprofile + ' ' + fname_out + ' ' + str(A_arr[i]) + ' ' + str(nRandom)
    print(command_genSim)
    os.system(command_genSim)

    fname_joblist = directory + '/' + fname_prefix + str(A_arr[i]) + '_joblist.txt'
    command_create_scripts = 'python make_job_scripts_pleiades.py ' + fname_joblist + ' ' + directory
    print(command_create_scripts)
    os.system(command_create_scripts)

    command_submit = 'python submit-jobs.py ' + directory
    print(command_submit)
    os.system(command_submit)

    fscan.write(directory + '/' + fname_matrix + '\n')

    numjobs = int(subprocess.check_output('squeue | grep "kyriacou" | wc -l', shell=True))
    while numjobs > 0:
        time.sleep(10)
        print('waiting, N-jobs left: ', numjobs)