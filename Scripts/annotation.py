#Noah A. Legall
#Salvador Lab
#Description: Run for loop that will annotate
#Last Update: June 2nd 2020
import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in
import glob # grabs files by name and puts them in a list
import re # we can do regular expression features with this
import time # for time stamps

#0. functions for script
logger = lambda message: "[{}] {}".format(time.strftime('%a %H:%M:%S'),message)

#1. create automatically the submission script for qsub
print(logger("creating the annotation submission script"))

qsub_script = open("annotation.sh","w")
qsub_script.write(
"""
#!/bin/bash
#PBS -q batch
#PBS -N annotate
#PBS -l nodes=1:ppn=5
#PBS -l mem=25gb
#PBS -l walltime=2:00:00
#PBS -M noahausxsapelo2xdump@gmail.com
#PBS -m abe
#PBS -o /scratch/noahaus/noahaus_out
#PBS -e /scratch/noahaus/noahaus_err
#PBS -j oe


ml prokka/1.14.6-foss-2018a

ASSEM=${assembly}
OUT=${out}
PRE=${prefix}

prokka --proteins mbovis_annot.gb --outdir $OUT --prefix $PRE $ASSEM
"""
)
qsub_script.close()

#2. run a for loop that will submit every job with a different group of fastq reads
assembly = sorted(glob.glob("*.scaffold.fa"))
# something new i'm using. lambdas are anonymous functions that don't need a formal name. basically quick and dirty function creation
# saves LOC if the function is relatively simple. here, I'm creating a list of output directory names.
annot_dir = list(map(lambda assemb: re.sub('.scaffold.fa','_annot',assemb), assembly))
prefix = list(map(lambda assemb: re.sub('.scaffold.fa','',assemb), assembly))

#qsub -v reference=/path/to/reference.fa bash.sh
for i in range(len(assembly)):
    os.system("qsub -v \"assembly={},out={},prefix={}\" annotation.sh".format(assembly[i],annot_dir[i],prefix[i]))
    print(logger("{} is being annotated".format(assembly[i])))

os.remove("annotation.sh")
print(logger("All annotations submitted to the cluster"))
