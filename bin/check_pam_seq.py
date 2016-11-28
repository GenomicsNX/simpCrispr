#!/usr/bin/env python


# check if PAM sequence is correct , NGG or NAG
# also remove duplicates, i.e if there are more than one sgRNA then remove that

import sys
import pandas as pd

filename=sys.argv[1]
outprefix=sys.argv[2]


df=pd.read_csv(filename,delim_whitespace=True,header=None)
df['pam']=df.apply(lambda s: s[2][1:3]=="GG" or s[2][1:3]=="AG",axis=1)

# only consider the correct PAM sequences
df_new=df[df['pam']]
set_dup=set(df_new[df_new.duplicated([0])][0])
df['dup']=df.apply(lambda s: s[0] not in set_dup,axis=1)


vv_passed=df[df['pam'] & df['dup']]
vv_failed=df[~df['dup']]

show_cols=[0,1,2,3,4,5]

vv_passed[show_cols].to_csv('%s.guides' %outprefix,sep="\t",index=False,header=False)
vv_failed[show_cols].to_csv('%s.failed.guides' %outprefix,sep="\t",index=False,header=False)




