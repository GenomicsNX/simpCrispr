

import sys
import argparse


# given a dict of sequences output a file with all possible mutations 
def getFakeMutations(seq):

	bb=[]
	tag=[]

	a=list(seq.upper())
	for i in range(len(a)):
		for e in ['A','T','G','C']:
			s0=a[i]
			if(a[i]!=e):
				bb.append(list2str(replace(a[:],i,e)))
				tag.append('pos%d%s>%s' %(i,s0,e))
	#cout(bb,10)

	return [bb,tag]

# ---------------------------------------------------------------

def replace(a,i,e):
	b=a
	b[i]=e
	return b

# ---------------------------------------------------------------

def list2str(l,delim=''):
	return delim.join([str(x) for x in l])

# ---------------------------------------------------------------


parser = argparse.ArgumentParser(description='')
parser.add_argument('-i','--infile',dest='infile',help='prefix of fasta file of sequences')
parser.add_argument('-o','--out',dest='out',help='output file prefix')
args = parser.parse_args()


filename=args.infile
outname=args.out

cc=0
with open(outname,'w') as outf:
	with open(filename) as inf:
		for line in inf:
			data=line.split()				
			seqID=data[0]
			seq=data[1]
			[fakeseq,tag]=getFakeMutations(seq)
			for no2,e in enumerate(fakeseq):
				#outf.write('>%s\n' % (seqID))
				#outf.write('%s-%s\n' % (seqID,tag[no2]))
				outf.write('>%s-seq%d\n%s\n' % (seqID,cc,e.upper()))
				#outf.write('%s\n' %(e.upper()))
				cc=cc+1




