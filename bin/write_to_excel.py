


import sys
import argparse
import pandas as pd
import numpy as np


filename=sys.argv[1]
out=sys.argv[2]
n_top_primer=sys.argv[3]
dn_primer=sys.argv[4]
up_primer=sys.argv[5]



writer=pd.ExcelWriter('%s.xlsx' %out,engine='xlsxwriter')


try:

	df=pd.read_csv('%s.sorted.guides' %out,delim_whitespace=True,header=None)
	col_names=['sgRNA ID','Guide Sequence (5->3)','PAM','Chromosome','Strand','Position (0-based)','# of Offtargets','Score']
	col_widths=[20,25,5,15,10,20,20,10]
	df.columns=col_names
	df.sort_values(by=['Score'],ascending=False,inplace=True)
	df.to_excel(writer,sheet_name='sgRNA',index=False,header=True,startrow=0,startcol=1)

	i=1
	for w in col_widths: 
		writer.sheets['sgRNA'].set_column(i,i,w)
		i+=1

except:
	pass

try:
	df=pd.read_csv('%s.failed.guides' %out,delim_whitespace=True,header=None)
	df.to_excel(writer,sheet_name='Removed sgRNA')
	col_names=['sgRNA ID','Guide Sequence (5->3)','PAM','Chromosome','Strand','Position (0-based)']
	col_widths=[20,25,10,15,10,20]
	df.columns=col_names
	df.to_excel(writer,sheet_name='Removed sgRNA',index=False,header=True,startrow=0,startcol=1)

	i=1
	for w in col_widths: 
		writer.sheets['Removed sgRNA'].set_column(i,i,w)
		i+=1
except:
	pass


try:
	
	df=pd.read_csv('%s.sorted.offtargets'%out,delim_whitespace=True,header=None)
	col_names=['sgRNA ID','Guide Sequence (5->3)','Offtarget ID','Chromosome','Strand','Position (0-based)','Offtarget Sequence (5->3)',
				'PAM','Mismatch Positions','# of Mismatches','Offtarget Score']
	col_widths=[20,25,13,13,10,20,25,5,20,20,15]
	df.columns=col_names
	df.to_excel(writer,sheet_name='Offtargets',index=False,header=True,startrow=0,startcol=1)


	i=1
	for w in col_widths: 
		writer.sheets['Offtargets'].set_column(i,i,w)
		i+=1
except:
	pass


try:
	
	df=pd.read_csv('%s.offtargets.exons'%out,delim_whitespace=True,header=None)
	col_names=['sgRNA ID','Guide Sequence (5->3)','Offtarget ID','Chromosome','Strand','Position (0-based)','Offtarget Sequence (5->3)',
				'PAM','Mismatch Positions','# of Mismatches','Offtarget Score','Annotation']
	col_widths=[20,25,13,13,10,20,25,5,20,20,15,20]
	df.columns=col_names
	df.to_excel(writer,sheet_name='Offtargets-annotated',index=False,header=True,startrow=0,startcol=1)


	i=1
	for w in col_widths: 
		writer.sheets['Offtargets-annotated'].set_column(i,i,w)
		i+=1
except:
	pass


try:

	df=pd.read_csv('%s.top%s.primers'%(out,n_top_primer),delim_whitespace=True,header=None)
	col_names=['Offtarget ','Flanking Sequence (%s downstream, %s upstream)' %(dn_primer,up_primer)]
	col_widths=[100,100]
	df.columns=col_names
	df.to_excel(writer,sheet_name='Primer Regions',index=False,header=True,startrow=0,startcol=1)

	i=1
	for w in col_widths: 
		writer.sheets['Primer Regions'].set_column(i,i,w)
		i+=1

except:
	pass



writer.save()














