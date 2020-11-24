# python script to perform a sliding window of tajima's D for an alignment
import dendropy
from dendropy.calculate import popgenstat # for calculating tajima's D
import argparse
from Bio import AlignIO #to parse the FASTA alignment
from Bio import SeqIO
from multiprocessing import Pool

#arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input',dest='input',help='the FASTA multi-sequence alignment input')
#parser.add_argument('-o','--output',dest='output',help='prefix of the csv file output')
parser.add_argument('-w','--window',dest='window',help='the window size to calculate tajima\'s D')
parser.add_argument('-s','--step',dest='step',help='the step size to instruct where to place the next window')
parser.add_argument('-t','--threads',dest='threads',help='number of cores to enable parallelization')

args = parser.parse_args()

#functions
def tajima_run(region):
    start = region.split(",")[0]
    end = region.split(",")[1]
    with open("temp.fa","w") as handle:
        SeqIO.write(align[:,start:end],handle,"fasta")
        seqs = dendropy.DnaCharacterMatrix.get(path="temp.fa",schema="fasta")
        try:
            tajimaD = dendropy.calculate.popgenstat.tajimas_d(seqs)
            print("{},{},{}".format(start,end,tajimaD))
        except:
            print("{},{},NA".format(start,end))


def implementation(x):
    p = Pool(args.threads)
    p.map(tajima_run,x)

#script
print("reading in alignment")
align = AlignIO.read(args.input, "fasta")
align_len = len(align[1])
print("alignment read. beginning tajima's D analysis")

regions = []
for i in range(0,align_len-(args.window-1),arg.step):
    regions.append("{},{}".format(i,i+args.window))
implementation(regions)
