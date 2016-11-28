


import sys
import argparse
from collections import defaultdict
import re
import itertools

# experimentally infered scores for MM locations, should be kept fixed
scoreMatrix=[0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583]


# mean pairwise distance between mismatches
# if only 1 MM then this would be 19, this is because to cancel the 
# effect of the term in the formula

def mean_d(mm):
	ss=0
	cc=0
	if len(mm)<2:
		return 19
	for x in range(len(mm)-1):
		ss+=abs(mm[x+1]-mm[x])
	return float(ss)/(len(mm)-1)


# if no MM, then the score is 100
# this is to shift the score matrix, basically for guide sequences smaller tahn 20bp, for example for 18bp it is shifted to left 2 units	
def hit_score(mm,shift=0):
	ss1=1
	for e in mm:
		ss1=ss1*(1-scoreMatrix[e-1+shift])
	ss2=float(1)/(((float(19-mean_d(mm))/19)*4)+1)
	ss3=float(1)/(len(mm)*len(mm))
	return ss1*ss2*ss3*100

	
def _f(a):
	
	b=map(lambda x:x.split(','),a)
	b=reduce(lambda x,y:x+y,b)
	b=list(set(b))
	b=map(int,b)
	b=sorted(b)
	
	return b



parser = argparse.ArgumentParser(description='')
parser.add_argument('-i','--infile',dest='infile',help='prefix of fasta file of sequences')
parser.add_argument('-o','--out',dest='out',help='output file prefix')
args = parser.parse_args()


filename=args.infile
outname=args.out

d=defaultdict(list)
s={}

with open(filename) as inf:
	with open(outname,'w') as outf:
		for no,line in enumerate(inf.readlines()):
			linev=line.split()
			if int(linev[9])==0:
				score=100
			else:
				shift=20-len(linev[1])
				score=hit_score(map(int,linev[8].split(",")),shift)
			outf.write('%s\t%.11f\n' %('\t'.join(linev),score))






