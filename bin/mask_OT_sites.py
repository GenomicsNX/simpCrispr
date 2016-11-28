

# just filter the Offtarget sites having lower case letter in any where
# lower case letter in reference means it is lower complexity and repetitive regions

import sys
import argparse
import re
import itertools


with open(sys.argv[1]) as fin:
	with open(sys.argv[2],'w') as fout:
		for line in fin:
			vv=line.split()
			if not vv[6].isupper():
				continue
			if not vv[7].isupper():
				continue
			fout.write(line)




