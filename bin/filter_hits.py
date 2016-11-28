

import sys
import argparse
import re




def rev(s):
	if s=='A':return 'T'
	if s=='C':return 'G'
	if s=='G':return 'C'
	if s=='T':return 'A'
	if s=='a':return 't'
	if s=='c':return 'g'
	if s=='g':return 'c'
	if s=='t':return 'a'
	return 'N'
# -----------------------------------------------------


def complement(ss):
	ss_new=[]
	for e in ss[::-1]:
		ss_new.append(rev(e))
	return ''.join(ss_new)
# -----------------------------------------------------



# extract sequence and PAM from the string
def _f(ss,strand,l=20):
	if strand=="-":
		pam=complement(ss[0:3])
		seq=complement(ss[3:len(ss)])
	else:
		pam=ss[l:l+3]
		seq=ss[0:l]
	return seq,pam
# -----------------------------------------------------


def get_MM(ss1,ss2):
	pp=1
	mm=[]
	for x,y in zip(ss1.upper(),ss2.upper()):
		if x!=y:
			mm.append(pp)
		pp+=1
	return ','.join(map(str,mm))

# -----------------------------------------------------




parser = argparse.ArgumentParser(description='')
parser.add_argument('-i','--infile',dest='infile',help='prefix of fasta file of sequences')
parser.add_argument('-o','--out',dest='out',help='output file prefix')
parser.add_argument('-g','--guide',dest='guide',help='output file prefix')


args = parser.parse_args()


filename=args.infile
outname=args.out
guidefile=args.guide


# load guide sequences
gg=dict()
with open(guidefile) as inf:
	for line in inf:
		linev=line.split()			
		#guide_id=linev[0].split("-")[0]
		guide_id=linev[0]
		#pam=linev[3]
		guide_seq=linev[1]
		gg[guide_id]=guide_seq



dd=[]
cc=0
with open(filename) as inf:
	try:
		with open(outname,'w') as outf:
			for line in inf:

				linev=line.split()
				chrom=linev[2]
				pos=int(linev[3])
				strand=linev[1]
				guide_id=linev[0].split("-")[0]
				guide_seq=gg[guide_id]
				seq,pam=_f(linev[9],strand,len(guide_seq))

				if not (pam[1:3].upper()=='GG' or pam[1:3].upper()=='AG'):
					continue
				key='%s:%s:%s' %(chrom,strand,pos)
				if key in dd:
					continue # it is already in the list
				else:
					dd.append(key)
			
				mm_loc=get_MM(guide_seq,seq)
				if len(mm_loc)==0:
					mm_loc='--'
					n_mm=0
				else:
					n_mm=len(mm_loc.split(","))
				if strand=="-":
					pos=pos-3
				outf.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(cc,guide_id,guide_seq,seq,pam,mm_loc,n_mm,chrom,strand,pos))
				cc+=1
	except:
		print sys.exc_info()
		print linev
