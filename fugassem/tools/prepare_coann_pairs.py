#!/usr/bin/env python

"""
Prepare co-ann paris for annotation matrix

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
import logging
import math

try:
	from fugassem import config,utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find FUGAsseM python package." +
		" Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Prepare co-ann paris for annotation matrix
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--input",
	                    help = '[REQUIRED] input annotation file',
	                    required = True)
	parser.add_argument('-t', "--header",
	                    help = '[OPTIONAL] whether include header or not in the input annotation file [Default: no]',
						choices= ['yes', 'no'],
						default = 'no')
	parser.add_argument('-o', "--output",
	                    help = 'output refined file',
	                    required = True)
	values = parser.parse_args()

	return values
# get_args



def collect_anno_info (raw_file, header):
	flag = 0
	anns = {}
	funcs = {}
	families = {}
	titles = {}
	for line in utilities.gzip_bzip2_biom_open_readlines(raw_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag == 0:
			flag = 1
			if header == "yes":
				for item in info:
					titles[info.index(item)] = item
				continue
		info = line.split("\t")
		myid = info[0]
		myf = info[-1]
		families[myid] = ""
		if not myid in anns:
			anns[myid] = {}
		anns[myid][myf] = 1
		if not myf in funcs:
			funcs[myf] = {}
		funcs[myf][myid] = ""
	# foreach line

	return families, funcs, anns


def prep_coann_pair (families, funcs, anns, outfile):
	# get co-annotation info
	coann = {}
	for myid1 in families.keys():
		if myid1 in anns:
			for myf in anns[myid1].keys():
				myv = anns[myid1][myf]
				if myf in funcs:
					for myid2 in funcs[myf].keys():
						if myid1 != myid2:
							mynew = myid1 + "\t" + myid2
							#if not mynew in coann:
							#	coann[mynew] = {}
							#coann[mynew][myf] = myv
							coann[mynew] = myv
					# foreach annotated gene for a given function
			# foreach annotated item for a given gene
		# if gene is annotated
	# foreach gene

	# write co-annotation into file
	open_out = open(outfile, "w")
	for myid1 in sorted(families.keys()):
		for myid2 in sorted(families.keys()):
			if myid1 == myid2:
				continue
			else:
				mynew = myid1 + "\t" + myid2
				if mynew in coann:
					myv = coann[mynew]
				else:
					#myv = 0
					continue
			mystr = mynew + "\t" + str(myv)
			open_out.write(mystr + "\n")
	open_out.close()
	
	if os.path.isfile(outfile):
		os.system("gzip " + outfile)

# prep_coann_pair


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start prepare_coann_pairs process")

	families, funcs, anns = collect_anno_info (values.input, values.header)
	prep_coann_pair (families, funcs, anns, values.output)

	config.logger.info ("Successfully finish prepare_coann_pairs process")

# end: main

if __name__ == '__main__':
	main()
