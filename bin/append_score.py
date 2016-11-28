
import sys
import argparse
import pandas as pd
import numpy as np


def get_score(hit_scores):
	try:
		sum=0
		for e in hit_scores:
			sum=sum+e
		return (100/(100+sum))*100
	except:
		pass
		#print "error",hit_scores

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i','--infile',dest='infile',help='prefix of fasta file of sequences')
parser.add_argument('-g','--guide',dest='guide',help='prefix of fasta file of sequences')
parser.add_argument('-o','--out',dest='out',help='prefix of fasta file of sequences')


args = parser.parse_args()
filename=args.infile
guide=args.guide
outname=args.out


gg=dict()
with open(guide) as fin:
	for line in fin:
		linev=line.split()
		gg[linev[0]]={'seq':linev[1]}


df=pd.read_csv(filename,delim_whitespace=True,header=None)

with open(outname,'w') as fout:
	for k in gg:
		try:
			vv=df[(df[0]==k)]
			seq=gg[k]['seq']
			ontarget=vv[ (vv[2]=="on-target") & (vv[9]==0)]
			offtarget=vv[vv[9]>0]	
			chrom=ontarget.iloc[0,3]
			st=ontarget.iloc[0,4]
			pos=ontarget.iloc[0,5]
			pam=ontarget.iloc[0,7]
			# get the OT where score <100 , score=100 means that on-target
			total_score=int(get_score(np.array(offtarget[10])))
			n_OT_sites=len(offtarget)							
			#fout.write('%s\t%s\t%s\t%s\t%s\n' %(k,seq,pam,int(score(np.array(vv[9]))),len(np.array(vv[9]))))
			fout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(k,seq,pam,chrom,st,pos,n_OT_sites,total_score))
		except:
			#print sys.exc_info()
			#print seq,chrom,st,pos,pam
			pass













