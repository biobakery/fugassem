#!/usr/bin/env python

"""
Prepare features for next layer of machine learning process

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
import argparse
import subprocess
import tempfile
import re
import logging
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import namedtuple
import math


try:
	from fugassem import utilities
	from fugassem import config
except:
	sys.exit ("fugassem is not install!")


def parse_arguments():
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "Prepare features for next layer of machine learning process\n",
		prog = "prepare_nextlayer_ml")
	parser.add_argument(
		"-r", "--raw",
		help = "[REQUIRED] input file including raw features used for previous ML\n",
		required = True)
	parser.add_argument(
		"-l", "--list",
		help = "[REQUIRED] input list file of previous ML output folders\n",
		required = True)
	parser.add_argument(
		'-e', "--extension",
	    help = '[OPTIONAL] extension of ML result file',
	    choices = [".xval.roc.tsv"],
	    default = ".xval.roc.tsv")
	parser.add_argument(
		'-a', "--label",
	    help = '[OPTIONAL] select predicted label name',
	    default = "yes")
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def collect_raw_feature (infile):
	"""
	Collect raw features
	Input:
		infile - feature file
	Output: function list and gene list
	"""

	config.logger.info('collect_raw_feature')

	funcs = {}
	titles = {}
	title = ""
	head_t = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (infile):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if head_t == 0:
			head_t = 1
			title = line
			for i in info:
				titles[i] = info.index(i)
			continue
		if re.search("^Feature\t", line):
			break
		myid = info[0]
		if re.search("__", myid):
			continue
		funcs[myid] = line

	return title, funcs


def collect_prediction_score (pred_file, data_type, prediction, predict_label):
	"""
	Collect prediction score from previous ML process
	Input:
		pred_file - previous learning process
		data_type - data type for this function prediction
		prediction - dictionary including all predictions for all functions
		predict_label - select predicted label
	Output: selected scores
	"""

	config.logger.info('collect_prediction_score')

	# collect prediction
	if not data_type in prediction:
		prediction[data_type] = {}
	open_file = open(pred_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myindex = 2
		while myindex < len(info):
			myv = info[myindex]
			if myindex % 2 == 0:
				if myv == predict_label:
					try:
						myw = info[myindex+1]
						prediction[data_type][myid] = myw
					except:
						config.logger.info("Weight value doesn't exist for this label: " + myv)
			myindex = myindex + 1
	open_file.close()


def prep_features4next (list_file, funcs, gene_header, extension, predict_label, outfile):
	"""
	Prepare features for next round of ML process
	Input:
		list_file - input list file of previous ML output folders
		funcs - a dictionary of function info
		gene_header - header of genes
		outfile - file name for output
	Output: write down collected features
	"""

	config.logger.info('prep_features4next')

	# collect predicted scores from previous ML process
	scores = {}
	open_list = open(list_file, "r")
	for mypath in open_list:
		mypath = mypath.strip()
		if not len(mypath):
			continue
		filelist = utilities.find_files(mypath, extension, None)
		for myfile in filelist:
			if not os.path.isfile(myfile):
				config.logger.info("Error! The prediction file doesn't exist: " + myfile)
				continue
			base = os.path.basename(myfile)
			data_type = re.sub(extension, "", base)
			if re.search("_results", mypath):
				mym = re.search("([^\/]+)_results", mypath)
				data_type = data_type + "__" + mym.group(1)
			collect_prediction_score(myfile, data_type, scores, predict_label)
		# foreach function
	# foreach data type
	open_list.close()

	# output features
	titles = {}
	info = gene_header.split("\t")
	genes = info[1:len(info)]
	for i in info:
		titles[i] = info.index(i)
	open_out = open(outfile, "w")
	open_out.write(gene_header + "\n")
	for myf in sorted(funcs.keys()):
		open_out.write(funcs[myf] + "\n")
	for myt in sorted(scores.keys()):
		mystr = myt
		for gene in genes:
			if gene	in scores[myt]:
				mystr = mystr + "\t" + scores[myt][gene]
			else:
				mystr = mystr + "\t0"
		open_out.write(mystr + "\n")
	mystr = "Feature\t" + "START\t" * len(genes)
	mystr = re.sub("\t$", "", mystr)
	open_out.write(mystr + "\n")
	open_out.close()


def main():
	args_value = parse_arguments()

	config.logger.info ("Start prepare_nextlayer_ml process ......")

	## Collect function info ##
	title, funcs = collect_raw_feature (args_value.raw)

	## Prepare features for next round of ML process ##
	prep_features4next (args_value.list, funcs, title, args_value.extension, args_value.label, args_value.output)

	config.logger.info ("Successfully finished prepare_nextlayer_ml process!")

if __name__ == '__main__':
	main()
