#!/usr/bin/env python

"""
Collect sequence-similarity based prediction results for each function

Copyright (c) 2021 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import sys
import os
import os.path
import re
import argparse
import math
import numpy as np
import logging


try:
	from fugassem import utilities
	from fugassem import config
except:
	sys.exit ("fugassem is not install!")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Collect sequence-similarity based prediction results for each function
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--input",
	                    help = '[REQUIRED] gene by function sequence-similarity matrix',
	                    required = True)
	parser.add_argument('-t', "--trim",
	                    help = '[OPTIONAL] trim off the suffix of each function, e.g. "__seqSimilarity"',
						default = None)
	parser.add_argument('-o', "--output",
	                    help = '[REQUIRED] output file name',
	                    required = True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect seqSimilarity info
#==============================================================
def collect_seq_info (raw_file, trim, outfile):
	open_file = open(raw_file, "r")
	open_out = open(outfile, "w")
	info = open_file.readline().strip().split("\t")
	mytitle = info[0]
	for item in info:
		if info.index(item) == 0:
			continue
		if trim:
			item = re.sub(trim, "", item)
			mytitle = mytitle + "\t" + item
	open_out.write(mytitle + "\n")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		open_out.write (line + "\n")
	open_file.close()
	open_out.close()


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start collect_seqSimilarity_results process")

	collect_seq_info (values.input, values.trim, values.output)
	
	config.logger.info ("Successfully finish collect_seqSimilarity_results")

# end: main

if __name__ == '__main__':
	main()
