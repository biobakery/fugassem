#!/usr/bin/env python

"""
Prepare sequence similarity information as evidence for prediction

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
		description = "Prepare sequence similarity information as evidence for prediction\n",
		prog = "prepare_seqSimilarity")
	parser.add_argument(
		"-a", "--family",
		help = "[REQUIRED] protein family file\n",
		required = True)
	parser.add_argument(
		"-t", "--type",
		help = "[OPTIONAL] type of sequence similarity info [Default: homology]\n",
		default = "homology")
	parser.add_argument(
		"-s", "--similarity",
		help = "[REQUIRED] sequence similarity info file\n",
		required = True)
	parser.add_argument(
		"-f", "--function",
		help = "[REQUIRED] original function annotation file\n",
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] name of output file\n",
		required = True)

	return parser.parse_args()


def collect_families (family_file):
	"""
	Collect protein families
	Input: the file name of protein family
	Output: families = {cluster1, cluster2, ...}
	"""
	families = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (family_file):
		line = line.strip()
		if not len(line):
			continue
		flag = flag + 1
		if flag == 1:
			continue
		info = line.split("\t")
		families[info[0]] = ""

	return families


def collect_similarity_unit (simi_file):
	"""
	Collect sequence similarity based on grouped units
	Input: the file name of sequence similarity
	Output: pair = {{c1:c2}, ...}
	"""
	units = {}
	maps = {}
	pairs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (simi_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myu = info[1]
		if not myid in maps:
			maps[myid] = {}
		maps[myid][myu] = ""
		if not myu in units:
			units[myu] = {}
		units[myu][myid] = ""

	for myid1 in maps.keys():
		for myu in maps[myid1].keys():
			if myu in units:
				for myid2 in units[myu].keys():
					if myid1 != myid2:
						if not myid1 in pairs:
							pairs[myid1] = {}
						if not myid2 in pairs[myid1]:
							pairs[myid1][myid2] = 1
						else:
							pairs[myid1][myid2] += 1
	maps = {}
	units = {}

	return pairs


def collect_function (func_file):
	"""
	Collect function maps
	Input: the file name of function annotation
	Output: func = {{f1: c1, c2}, ...}
	"""
	funcs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (func_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myf = info[1]
		if not myid in funcs:
			funcs[myid] = {}
		funcs[myid][myf] = ""

	return funcs


def build_function_matrix (families, pairs, funcs, outfile):
	"""
	Build function association matrix based on sequence similarity
	Input: families, pairs, funcs
	Output: output function association info file
	"""

	anns = {}
	myfuncs = {}
	for myid in families.keys():
		#if myid in funcs:
		#	if not myid in anns:
		#		anns[myid] = {}
		#	for myf in funcs[myid].keys():
		#		myfuncs[myf] = ""
		#		anns[myid][myf] = 1
		if myid in pairs:
			for myid2 in pairs[myid].keys():
				if myid == myid2:
					continue
				myv = float(pairs[myid][myid2])
				if myid2 in funcs:
					for myf in funcs[myid2].keys():
						if not myid in anns:
							anns[myid] = {}
						if not myf in anns[myid]:
							myfuncs[myf] = ""
							anns[myid][myf] = myv
						if float(anns[myid][myf]) < float(myv):
							anns[myid][myf] = myv
			# foreach pair
		# if exits pair
	# foreach protein family

	# write info into file
	mytitle = "ID"
	for myf in sorted(myfuncs.keys()):
		mytitle = mytitle + "\t" + myf + "__seqSimilarity"
	open_out = open(outfile, "w")
	open_out.write(mytitle + "\n")
	for myid in sorted(families.keys()):
		mystr = myid
		for myf in sorted(myfuncs.keys()):
			if myid in anns:
				if myf in anns[myid]:
					mystr = mystr + "\t" + str(anns[myid][myf])
				else:
					mystr = mystr + "\t0"
			else:
				mystr = mystr + "\t0"
		open_out.write(mystr + "\n")
	open_out.close()

# function build_function_matrix


def main():
	args_value = parse_arguments()

	config.logger.info ("Start prepare_seqSimilarity process......")

	## get info ##
	families = collect_families (args_value.family)
	funcs = collect_function (args_value.function)

	#units = ["pfam", "uniref50", "contig", "DDI", "homology"]
	paris = {}
	#if args_value.type in units:
	pairs = collect_similarity_unit (args_value.similarity)
	build_function_matrix(families, pairs, funcs, args_value.output)

	config.logger.info("Successfully finished prepare_seqSimilarity process!")

if __name__ == '__main__':
	main()
