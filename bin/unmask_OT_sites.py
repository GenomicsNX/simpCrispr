

# replace lower case letter to upper case

import sys
import argparse
import re


with open(sys.argv[1]) as fin:
	with open(sys.argv[2],'w') as fout:
		for line in fin:
			vv=line.split()
			vv[6]=vv[6].upper()
			vv[7]=vv[7].upper()
			fout.write('\t'.join(vv)+'\n')





