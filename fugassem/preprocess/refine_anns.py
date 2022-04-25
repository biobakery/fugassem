#!/usr/bin/env python

"""
Refine the annotation info using new list info

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
import argparse
import subprocess
import tempfile
import re
import numpy as np
import math
import statistics
import logging


try:
	from fugassem import utilities
	from fugassem import config
except:
	sys.exit ("fugassem is not installed!")


def parse_arguments():
	"""
	Parse the arguments from the user
	"""

	parser = argparse.ArgumentParser(
		description = "Refine the annotation info using new list info\n",
		prog = "refine_anns")
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] known annotations of features\n",
		required = True)
	parser.add_argument(
		"-l", "--list",
		help = "[REQUIRED] a white list of functions\n",
		required = True)
	parser.add_argument(
		"-t", "--type",
		help = "[REQUIRED] type of function\n",
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] output refined file\n",
		required = True)

	return parser.parse_args()


def collect_function_list (func_file, func_type):
	"""
	Collect the white list of functions
	Input: the file name of the function list
	Output: myfuncs = {func1, func2, ...}
	"""

	config.logger.info ('collect_function_list')

	myfuncs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (func_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		info2 = info[0].split(":")
		mytype = func_type
		if info2[0] == "GO":
			myid = ":".join(info2[0:2])
			info[0] = re.sub("\s+\[BP\]", "[BP]", info[0])
			info[0] = re.sub("\s+\[MF\]", "[MF]", info[0])
			info[0] = re.sub("\s+\[CC\]", "[CC]", info[0])
			if re.search("\[BP\]", info[0]):
				mytype = "BP"
			if re.search("\[MF\]", info[0]):
				mytype = "MF"
			if re.search("\[CC\]", info[0]):
				mytype = "CC"
		else:
			myid = info2[0]
		if len(info) > 1:
			myindex = 1
			while myindex < len(info):
				myv = info[myindex]
				if not myv in myfuncs:
					myfuncs[myv] = {}
				myfuncs[myv][myid] = myv + "\t" + info[0] + "__" + mytype
				myindex = myindex + 1

	return myfuncs


def collect_annotation (ann_file):
	"""
	Collect the annotations of features
	Input:
		ann_file - the file name of annotations: from feature to function
		funcs - white list of functions
	Output: anns = {func1:[f1, f2, f3], func2:[f1], ...}
	"""

	config.logger.info ('collect_annotation')

	anns = {}
	for line in utilities.gzip_bzip2_biom_open_readlines(ann_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		info2 = info[-1].split(":")
		if info2[0] == "GO":
			myid_new = ":".join(info2[0:2])
		else:
			myid_new = info2[0]
		if not myid in anns:
			anns[myid] = {}
		anns[myid][myid_new] = line

	return anns

def refine_annotation (raw_ann, new_ann, outfile):
	"""
	Refine annotations
	Input:
		raw_ann - raw annotations
		new_ann - new annotations
	Output: output annotation
	"""

	config.logger.debug('refine_annotation')
	
	outfile1 = re.sub(".tsv$", ".simple.tsv", outfile)
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	
	"""
	for myid in sorted(raw_ann.keys()):
		for myf in sorted(raw_ann[myid].keys()):
			open_out.write(raw_ann[myid][myf] + "\n")
			tmp = myf.split("__")[0]
			open_out1.write(myid + "\t" + tmp + "\n")
		if myid in new_ann:
			for myf in sorted(new_ann[myid].keys()):
				if myf in raw_ann[myid]:
					continue
				else:
					open_out.write(new_ann[myid][myf] + "\n")
					tmp = myf.split("__")[0]
					open_out1.write(myid + "\t" + tmp + "\n")
	"""

	for myid in sorted(new_ann.keys()):
		#if myid in raw_ann:
		#	continue
		for myf in sorted(new_ann[myid].keys()):
			open_out.write(new_ann[myid][myf] + "\n")
			tmp = myf.split("__")[0]
			open_out1.write(myid + "\t" + tmp + "\n")
	open_out.close()
	open_out1.close()


def main():
	args_value = parse_arguments()

	config.logger.info ("Start refine_anns process......")

	new_ann = collect_function_list (args_value.list, args_value.type)
	raw_ann = collect_annotation (args_value.annotation)
	refine_annotation (raw_ann, new_ann, args_value.output)

	config.logger.info("Successfully finished refine_anns process!")


if __name__ == '__main__':
	main()
