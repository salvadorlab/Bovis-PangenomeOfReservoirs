#Noah A. Legall
#Salvador Lab
#Description: Generate Statistics for Genome Assemblies using QUAST
#Last Update: May 14th 2020

import glob # for gathering similarly named files
import os # interact with the OS that this is running in
import re # regular expression utilization
import functools # use for the reduce function.
import time # for time stamps

logger = lambda message: "[{}] {}".format(time.strftime('%a %H:%M:%S'),message)

print(logger("beginning the generation of statistics"))
qsub_script = open("assembly_stats.sh","w")
qsub_script.write(
"""
#!/bin/bash
#PBS -q batch
#PBS -N assembly_stats
#PBS -l nodes=1:ppn=20
#PBS -l mem=20gb
#PBS -l walltime=10:00:00
#PBS -M noahausxsapelo2xdump@gmail.com
#PBS -m abe
#PBS -o /scratch/noahaus/noahaus_out
#PBS -e /scratch/noahaus/noahaus_err
#PBS -j oe

ISOLATES=${isolates}
REF=${ref}

cd $PBS_O_WORKDIR
ml QUAST/5.0.2-foss-2018a-Python-2.7.14
quast.py -t 15 -r $REF -o assembly_stats $ISOLATES
"""
)
qsub_script.close()

###QUAST
#create a string of contigs to grab stats from
scaffold = sorted(glob.glob("*.scaffold.fa"))
# reduce() uses a sequential function on a entire list. here I use it to aggregate strings.
scaffold_string = reduce(lambda final, current: final + "{} ".format(current),scaffold)
print(logger("list of isolates to input: {}".format(scaffold_string)))

os.system("qsub -v \"isolates={},ref=./NC_002945v4.fasta\" assembly_stats.sh")
os.remove('assembly_stats.sh')
print(logger("all isolates submitted for statistic generation."))
