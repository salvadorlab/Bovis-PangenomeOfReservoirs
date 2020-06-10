#Noah A. Legall
#Salvador Lab
#Description: Generate Statistics for Genome Assemblies using QUAST
#Last Update: May 14th 2020

import glob # for gathering similarly named files
import os # interact with the OS that this is running in
import re # regular expression utilization
from functools import reduce # use for the reduce function.
import time #for time stamps

logger = lambda message: "[{}] {}".format(time.strftime('%a %H:%M:%S'),message)

print(logger("beginning the generation of statistics"))
###QUAST
#create a string of contigs to grab stats from
scaffold = sorted(glob.glob("*.scaffold.fa"))
# reduce() uses a sequential function on a entire list. here I use it to aggregate strings.
scaffold_string = reduce(lambda final, current: final + " {}".format(current),scaffold)
print(logger("list of isolates to input: {}".format(scaffold_string)))
os.system("quast.py -t 15 -r ./NC_002945v4.fasta -o assembly_stats {}".format(scaffold_string))
print(logger("all isolates submitted for statistic generation."))
