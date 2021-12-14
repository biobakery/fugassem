#!/usr/bin/env python

"""
Prepare annotation used for prediction

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
		description = "Prepare annotation information of features\n",
		prog = "prepare_annotation")
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] metawibele finalized annotation file\n",
		required = True)
	parser.add_argument(
		"-m", "--mapping",
		help = "[OPTIONAL] function name mapping file [Default: None]\n",
		default = None)
	parser.add_argument(
		"-t", "--type",
		help = "[REQUIRED] annotation type\n",
		choices = ["UniRef90_GO", "UniRef90_GO_BP", "UniRef90_GO_CC", "UniRef90_GO_MF", "UniRef90_COG", "UniRef90_eggNOG", "UniRef90_KEGG-KOs",
		           "InterProScan_PfamDomain", "Denovo_transmembrane", "Denovo_signaling", "DOMINE_interaction"],
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] output file\n",
		required = True)

	return parser.parse_args()


def collect_name_map (map_file):
	"""
	Collect function names from the mapping file
	Input: the file name of the function names
	Output: funcs = {func1, func2, ...}
	"""
	names = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		names[info[0]] = info[1]

	return names


def collect_function_info (ann_file, func_type, names, outfile):
	"""
	Collect the white list of functions
	Input: the file name of the function list
	Output: output feature2function info file
	"""

	outfile1 = re.sub(".tsv", ".simple.tsv", outfile)
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	titles = {}
	open_file = open(ann_file, "r")
	line = open_file.readline().strip()
	info = line.split("\t")
	for item in info:
		titles[item] = info.index(item)
	hits = {}
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mytype = info[titles["feature"]]
		mytype = re.sub("\(", "_", mytype)
		mytype = re.sub("\)", "", mytype)
		myann = info[titles["annotation"]]
		mycat = info[titles["category"]]
		if mytype == "quality":
			if myann == "good":
				hits[myid] = ""
			continue
		if mytype != func_type:
			continue
		if myann == "NA":
			continue
		if not myid in hits:
			continue
		tmp = myann.split(";")
		for item in tmp:
			if re.search("GO\:", item):
				mym = re.search("(GO\:[^\]]+)", item)
				item = mym.group(1)
			myname = item
			open_out1.write(myid + "\t" + myname + "\n")
			if item in names:
				myname = myname + ":" + names[item]
			myname = myname + "__" + func_type
			open_out.write(myid + "\t" + myname + "\n")

	open_out.close()
	open_out1.close()

# function collect_function_info


def main():
	args_value = parse_arguments()

	config.logger.info ("Start prepare_annotation process......")

	## get info ##
	if args_value.mapping:
		names = collect_name_map (args_value.mapping)
	else:
		names = {}
	collect_function_info (args_value.annotation, args_value.type, names, args_value.output)

	config.logger.info("Successfully finished prepare_annotation process!")

if __name__ == '__main__':
	main()
