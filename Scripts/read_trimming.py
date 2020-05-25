#Noah A. Legall
#Salvador Lab
#Description: Trim all the isolates using trimmomatic
#Last Update: May 11th 2020
import sys # use to access arguments
import os # use in order to call commands from the terminal script is called in
import glob # grabs files by name and puts them in a list
import re # we can do regular expression features with this
import time # for time stamps

#0. functions for script
logger = lambda message: "[{}] {}".format(time.strftime('%a %H:%M:%S'),message)

#1. create automatically the submission script for qsub
print(logger("creating the read trimming submission script"))
qsub_script = open("read_trimming.sh","w")
qsub_script.write(
"""
#!/bin/bash
#PBS -q batch
#PBS -N read_trimming
#PBS -l nodes=1:ppn=10
#PBS -l mem=25gb
#PBS -l walltime=2:00:00
#PBS -M noahausxsapelo2xdump@gmail.com
#PBS -m abe
#PBS -o /scratch/noahaus/noahaus_out
#PBS -e /scratch/noahaus/noahaus_err
#PBS -j oe

R1=${r1}
R2=${r2}
R1_PAIRED_OUT=${r1_paired_out}
R2_PAIRED_OUT=${r2_paired_out}
R1_UNPAIRED_OUT=${r1_unpaired_out}
R2_UNPAIRED_OUT=${r2_unpaired_out}

cd $PBS_O_WORKDIR
module load Trimmomatic/0.36-Java-1.8.0_144

java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 8 $R1 $R2 $R1_PAIRED_OUT $R1_UNPAIRED_OUT \
$R2_PAIRED_OUT $R2_UNPAIRED_OUT \
ILLUMINACLIP:/usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

rm $R1_UNPAIRED_OUT $R2_UNPAIRED_OUT
"""
)
qsub_script.close()

#2. run a for loop that will submit every job with a different group of fastq reads
r1 = sorted(glob.glob('*_1.fastq'))
r2 = sorted(glob.glob('*_2.fastq'))
# something new i'm using. lambdas are anonymous functions that don't need a formal name. basically quick and dirty function creation
# saves LOC if the function is relatively simple. here, I'm creating a list of output directory names.
r1_paired = list(map(lambda raw: re.sub('_1.fastq','.paired.trimmed.R1.fastq',raw), r1))
r1_unpaired = list(map(lambda raw: re.sub('_1.fastq','.unpaied.trimed.R1.fastq',raw),r1))
r2_paired = list(map(lambda raw: re.sub('_2.fastq','.paired.trimmed.R2.fastq',raw), r2))
r2_unpaired = list(map(lambda raw: re.sub('_2.fastq','.unpaied.trimed.R2.fastq',raw),r2))

#qsub -v reference=/path/to/reference.fa bash.sh
for i in range(len(r1)):
    os.system("qsub -v \"r1={},r2={},r1_paired_out={},r1_unpaired_out={},r2_paired_out={},r2_unpaired_out={}\" read_trimming.sh".format(r1[i],r2[i],r1_paired[i],r1_unpaired[i],r2_paired[i],r2_unpaired[i]))
    print(logger("Read trimming performed on {} {}".format(r1[i],r2[i])))

os.remove("read_trimming.sh")
print(logger("All isolates submitted to the cluster"))
