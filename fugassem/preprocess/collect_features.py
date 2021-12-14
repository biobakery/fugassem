#!/usr/bin/env python

"""
Prepare features for machine learning

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
	os.sys.stderr ("fugassem is not install!")


def parse_arguments():
	"""
	Parse the arguments from the user
	"""

	parser = argparse.ArgumentParser(
		description = "Prepare feature information of ML\n",
		prog = "collect_features")
	parser.add_argument(
		"-r", "--raw",
		help = "[REQUIRED] raw feature file\n",
		required = True)
	parser.add_argument(
		"-n", "--new",
		help = "[REQUIRED] new feature file for combining\n",
		required = True)
	parser.add_argument(
		"-t", "--header",
		help = "[OPTIONAL] whether the new-feature file has header or not [Default: None] \n",
		default = None)
	parser.add_argument(
		"-f", "--type",
		help = "[REQUIRED] data type of new feature\n",
		choices = ["matrix", "homology", "vector", "function"],
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] output file\n",
		required = True)

	return parser.parse_args()


def collect_feature_matrix (new_file, header):
	"""
	Collect new features
	Input: the file name of new features for genes
	Output: funcs = {func1, func2, ...}
	"""
	features = {}
	titles = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (new_file):
		line = line.strip()
		if not len(line):
			continue
		if header:
			if header != "no":
				if flag == 0:
					flag = 1
					info = line.split("\t")
					for i in info:
						titles[info.index(i)] = re.sub(":", "_", i)
					continue
		flag = 1
		info = line.split("\t")
		myid = info[0]
		myindex = 1
		while myindex < len(info):
			myf = titles[myindex]
			myv = info[myindex]
			if not myf in features:
				features[myf] = {}
			features[myf][myid] = myv
			myindex = myindex + 1

	return features

def collect_feature_hom (new_file, header):
	"""
	Collect new features
	Input: the file name of new features for genes
	Output: funcs = {func1, func2, ...}
	"""
	features = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (new_file):
		line = line.strip()
		if not len(line):
			continue
		if header:
			if header != "no":
				if flag == 0:
					flag = 1
					continue
		flag = 1
		info = line.split("\t")
		info[1] = re.sub(":", "_", info[1])
		if not info[1] in features:
			features[info[1]] = {}
		features[info[1]][info[0]] = 1

	return features


def collect_feature_func (new_file, header):
	"""
	Collect new features
	Input: the file name of new features for genes
	Output: funcs = {func1, func2, ...}
	"""
	features = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (new_file):
		line = line.strip()
		if not len(line):
			continue
		if header:
			if header != "no":
				if flag == 0:
					flag = 1
					continue
		flag = 1
		info = line.split("\t")
		info[1] = re.sub(":", "_", info[1])
		if not info[1] in features:
			features[info[1]] = {}
		features[info[1]][info[0]] = "yes"

	return features


def combine_features (raw_file, features, func_type, outfile):
	"""
	Combine new features to the raw file
	Input: the file name of raw features, a dictionary of new features
	Output: write out the combined features into file
	"""

	meta = ["function", "vector"]
	open_out = open(outfile, "w")
	open_file = open(raw_file, "r")
	line = open_file.readline().strip()
	info = line.split("\t")
	items = {}
	for i in info:
		items[i] = ""
	names = sorted(features.keys())
	add_boundary = 0
	if func_type in meta:
		if len(info) > 1:
			if "Feature" in items:
				open_out.write(info[0] + "\t" + "\t".join(names) + "\t" + "\t".join(info[1:len(info)])+ "\n")
			else:
				add_boundary = 1
				open_out.write(info[0] + "\t" + "\t".join(names) + "\t" + "Feature\t" + "\t".join(info[1:len(info)])+ "\n")
		else:
			add_boundary = 1
			open_out.write(info[0] + "\t" + "\t".join(names) + "\t" + "Feature\n")
	else:
		if "Feature" in items:
			open_out.write(line + "\t" + "\t".join(names) + "\n")
		else:
			add_boundary = 1
			if len(info) > 1:
				open_out.write(info[0] + "\t" + "Feature\t" + "\t".join(info[1:len(info)]) + "\t" + "\t".join(names) + "\n")
			else:
				open_out.write(info[0] + "\t" + "Feature\t" + "\t".join(names) + "\n")
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		mystr = ""
		if func_type in meta:
			for i in names:
				if myid in features[i]:
					myv = features[i][myid]
				else:
					if func_type == "function":
						myv = "no"
					else:
						myv = str(0)
				if mystr == "":
					mystr = myv
				else:
					mystr = mystr + "\t" + myv
			if len(info) > 1:
				if add_boundary == 1:
					myline = myid + "\t" + mystr + "\tSTART" + "\t" + "\t".join(info[1:len(info)])
				else:
					myline = myid + "\t" + mystr + "\t" + "\t".join(info[1:len(info)])
			else:
				myline = myid + "\t" + mystr + "\tSTART"
		else:
			for i in names:
				if myid in features[i]:
					myv = str(features[i][myid])
				else:
					#if func_type == "homology":
					myv = str(0)
				if mystr == "":
					mystr = myv
				else:
					mystr = mystr + "\t" + myv
			if add_boundary == 1:
				if len(info) > 1:
					myline = myid + "\tSTART\t" + "\t".join(info[1:len(info)]) + "\t" + mystr
				else:
					myline = myid + "\tSTART\t" + mystr
			else:
				myline = line + "\t" + mystr
		open_out.write(myline + "\n")
	open_out.close()

# function combine_features


def main():
	values = parse_arguments()

	config.logger.info ("Start collect_features process......")

	## get info ##
	if values.type == "function":
		features = collect_feature_func (values.new, values.header)
	if values.type == "homology":
		features = collect_feature_hom (values.new, values.header)
	if values.type == "matrix" or values.type == "vector":	
		features = collect_feature_matrix (values.new, values.header)
	combine_features(values.raw, features, values.type, values.output)


	config.logger.info("Successfully finished collect_features process!")

if __name__ == '__main__':
	main()
