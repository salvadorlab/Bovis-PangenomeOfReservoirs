#Noah A. Legall
#Salvador Lab
#Description: Run for loop that will submit a submission job for each accession in an accession list
#Last Update: May 4th 2020
import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in
import time # for time stamps

#1. create automatically the submission script for qsub
print("[{}] beginning the accession list download".format(time.strftime('%a %H:%M:%S')))
qsub_script = open("download_fastq.sh","w")
qsub_script.write(
"""
#!/bin/bash
#PBS -q batch
#PBS -N download_fastq
#PBS -l nodes=1:ppn=5
#PBS -l mem=10gb
#PBS -l walltime=2:00:00
#PBS -M noahausxsapelo2xdump@gmail.com
#PBS -m abe
#PBS -o /scratch/noahaus/noahaus_out
#PBS -e /scratch/noahaus/noahaus_err
#PBS -j oe

ACC=${accession}

cd $PBS_O_WORKDIR
ml SRA-Toolkit/2.9.1-centos_linux64

fastq-dump --split-3 $ACC
"""
)
qsub_script.close()

#2. run a for loop that will submit every job with a different accession
acc_file = sys.argv[1].strip()
acc_list = []
acc = open(acc_file,'r')

for line in acc:
    acc_list.append(line.strip())

for accession in acc_list:
    os.system("qsub -v \"accession={}\" download_fastq.sh".format(accession))
    print("[{}] accession {} has been submitted".format(time.strftime('%a %H:%M:%S'),accession))

os.remove("download_fastq.sh")
print("[{}] all accessions are submitted to the cluster".format(time.strftime('%a %H:%M:%S')))
