


import sys
import argparse


parser = argparse.ArgumentParser(description='')
parser.add_argument('-i','--infile',dest='infile',help='prefix of fasta file of sequences')
parser.add_argument('-o','--out',dest='out',help='output file prefix')
args = parser.parse_args()


filename=args.infile
outname=args.out


with open(outname,'w') as outf:
	with open(filename) as inf:
		for line in inf:
			linev=line.split()
			seqID=linev[0]
			chrom=linev[2]
			seq=linev[4]
			pos=linev[3]
			strand=linev[1]
			#outf.write('%s\t%d\t%d\t%s\n' %(chrom,int(loc),int(loc)+len(seq)+3,seq))
			if (strand=='+'):
				outf.write('%s\t%d\t%d\t%s\t%s + \n' %(chrom,int(pos),int(pos)+len(seq)+3,seqID,0))
			elif (strand=='-'):
				outf.write('%s\t%d\t%d\t%s\t%s - \n' %(chrom,int(pos)-3,int(pos)+len(seq),seqID,0))
