#!/usr/bin/env python


import argparse
import re



def getComplementary(s,rever=True): 
	s=s.upper() 
	a=list(s)
	b=[] 
	for e in a: 
		if(e=='G'): b.append('C') 
		elif (e=='C'):b.append('G') 
		elif (e=='A'):b.append('T') 
		elif (e=='T'):b.append('A') 
	if rever: 
		return ''.join(reversed(b)) 
	else: 
		return ''.join(b) 


def write_to_file(bb,outf,pam,maxl=20): 
	
	global counter
	fw_strand=[]
	rv_strand=[]
	start=0


	if pam=="NGG":
		findPos(bb,fw_strand,'GG')
		findPos(bb,rv_strand,'CC')
	elif pam=="NAG":
		findPos(bb,fw_strand,'AG')
		findPos(bb,rv_strand,'CT')

		
	for p in fw_strand:
		if p<maxl: continue
		if start+p-maxl-1<0:continue
		#outf.write('>%s%.3d\n%s\t%s\t%s\n' % ("Guide",counter,bb[p-maxl-1:p-1],bb[p-1:p+2],'+'))
		outf.write('>%s%.3d\n%s\n' % ("sgRna",counter,bb[p-maxl-1:p-1]))
		counter+=1

	for p in rv_strand:
		if int(p+3+maxl)>len(bb):continue
		#outf.write('>%s%.3d\n%s\t%s\t%s\n' % ("Guide",counter,getComplementary(bb[p+3:p+maxl+3]),getComplementary(bb[p:p+3]),'-'))
		outf.write('>%s%.3d\n%s\n' % ("sgRna",counter,getComplementary(bb[p+3:p+maxl+3])))
		counter+=1



def findPos(bb,ll,ss):
	p=0
	while True:
		p=bb.find(ss,p)
		if(p==-1):break
		ll.append(p)
		p+=1
	




parser = argparse.ArgumentParser(description='')
parser.add_argument('-f','--fasta',dest='fasta',help='Fasta file of sequence')
parser.add_argument('-l','--maxl',dest='maxl',default=20,type=int,help='Length of Guiding Sequences [Default:20]')
parser.add_argument('-o','--out',dest='out',help='Output Name')
parser.add_argument('-p','--pam',dest='pam',help='PAM Sequence NGG or NGA [Default: Both]')



args = parser.parse_args()
maxl=args.maxl
fasta_file=args.fasta
out=args.out
counter=1

sequences=[]

for line in open(fasta_file).readlines():
	if re.search('^>',line):
		continue
	
	s_temp=line.strip().upper()
	if len(s_temp)<maxl:
		continue

	sequences.append(s_temp)



#if len(labels)!=len(sequences):
#labels=['']*len(sequences)
#seq_string=''.join(ss)




# make guides text file
with open('%s' %out,'w') as outf:
	for seq_string in sequences:
		if args.pam!="both":
			write_to_file(seq_string,outf,args.pam,maxl)
		else:
			write_to_file(seq_string,outf,'NGG',maxl)
			write_to_file(seq_string,outf,'NAG',maxl)




