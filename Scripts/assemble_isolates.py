#Noah A. Legall
#Salvador Lab
#Description: Run for loop that will submit a submission job for each assembly
#Last Update: May 10th 2020
import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in
import glob # grabs files by name and puts them in a list
import re # we can do regular expression features with this
import time # for time stamps

#0. functions for script
logger = lambda message: "[{}] {}".format(time.strftime('%a %H:%M:%S'),message)

#1. create automatically the submission script for qsub
print(logger("creating the assembly submission script"))

qsub_script = open("assemble_isolates.sh","w")
qsub_script.write(
"""
#!/bin/bash
#PBS -q batch
#PBS -N assemble_isolates
#PBS -l nodes=1:ppn=5
#PBS -l mem=25gb
#PBS -l walltime=2:00:00
#PBS -M noahausxsapelo2xdump@gmail.com
#PBS -m abe
#PBS -o /scratch/noahaus/noahaus_out
#PBS -e /scratch/noahaus/noahaus_err
#PBS -j oe

R1=${r1}
R2=${r2}
OUT=${out}
SCAF=${scaffold}

cd $PBS_O_WORKDIR
ml add spades/3.12.0-k_245

spades.py --only-assembler --careful -1 $R1 -2 $R2 -o $OUT
cd $OUT
mv scaffolds.fasta $SCAF
"""
)
qsub_script.close()

#2. run a for loop that will submit every job with a different group of fastq reads
r1 = sorted(glob.glob('*trimmed.R1.fastq'))
r2 = sorted(glob.glob('*trimmed.R2.fastq'))
# something new i'm using. lambdas are anonymous functions that don't need a formal name. basically quick and dirty function creation
# saves LOC if the function is relatively simple. here, I'm creating a list of output directory names.
scaffold_file = list(map(lambda trim: re.sub('.paired.trimmed.R1.fastq','.scaffold.fa',trim), r1))
assemble_dir = list(map(lambda trim: re.sub('.paired.trimmed.R1.fastq','_assembly',trim), r1))

#qsub -v reference=/path/to/reference.fa bash.sh
for i in range(len(r1)):
    os.system("qsub -v \"r1={},r2={},out={},scaffold={}\" assemble_isolates.sh".format(r1[i],r2[i],assemble_dir[i],scaffold_file[i]))
    print(logger("DeBrujin Graph Assembly occuring now between {} and {}. Output resides in {} directory".format(r1[i],r2[i],assemble_dir[i])))

os.remove("assemble_isolates.sh")
print(logger("All isolates submitted to the cluster"))
