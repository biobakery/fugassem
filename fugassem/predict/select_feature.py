#!/usr/bin/env python

"""
Perform feature selection for machine learning

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
	from fugassem.common.corrmp import NaNCorrMp
except:
	sys.exit ("corrmp is not install!")
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
		description = "Perform feature selection for machine learning\n",
		prog = "select_feature")
	parser.add_argument(
		"-i", "--input",
		help = "[REQUIRED] input feature table for ML\n",
		required = True)
	parser.add_argument(
		"-r", "--redundancy",
		help = "[OPTIONAL] maximum correlation between features for filtering more redundant features [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-m", "--method",
		help = "[OPTIONAL] correlation methods, [ Default: Pearson ]\n",
		choices = ["Pearson", "Spearman", "Kendall"],
		default = "Pearson")
	parser.add_argument(
		"-c", "--core",
		help = "[OPTIONAL] number of threads, [ Default: 1 ]\n",
		default = 1)
	parser.add_argument(
		"-n", "--information",
		help = "[OPTIONAL] minimum informative fraction of features for filtering less informative features [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-p", "--importance",
		help = "[OPTIONAL] minimum prediction importance of features from 1st round prediction for filtering less importance features [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-l", "--implist",
		help = "[OPTIONAL] input list file including the prediction importance files. It should be customized when '--importance' para is used [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def collect_feature (infile):
	"""
	Collect features
	Input:
		infile - feature file
	Output: function list and feature list
	"""

	config.logger.info('collect_feature')

	funcs = []
	features = {}
	open_file = open(infile, "r")
	title = open_file.readline().strip()
	flag = 0
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		if re.search("^Feature\t", line):
			funcs.append(line)
			flag = 1
			continue
		if flag == 0:
			funcs.append(line)
		if flag == 1:
			info = line.split("\t")
			features[info[0]] = "\t".join(info[1:len(info)])
	open_file.close()

	return title, funcs, features


def flt_redundant_feature (redu_level, features, genes, corr_method, cores):
	"""
	Filter redundant features
	Input:
		redu_level - maximum correlation value
		feature - a dictionary of feature info
		gene - gene arrays
		corr_method - correlation method
		cores - number of threads
	Output: selected features
	"""

	config.logger.info('flt_redundant_feature')

	# build data frame
	evidence_row = sorted(features.keys())
	try:
		evidence_table_row = namedtuple("evidence_table_row", evidence_row, rename = False)
	except:
		evidence_table_row = namedtuple("evidence_table_row", evidence_row, verbose = False, rename = False)
	feature_table = pd.DataFrame(index = genes, columns = evidence_table_row._fields)

	for myf in evidence_row:
		myvalue = features[myf]
		feature_table[myf] = myvalue.split("\t")

	# correlation
	if corr_method == "Pearson":
		corr, p_value, se = NaNCorrMp.calculate_with_pvalue_se (feature_table, n_jobs=cores, zero="lenient", transform=False) 
	else:
		corr_method = corr_method.lower()
		corr = abunds_table.corr (method = corr_method)
		p_value = None
		se = None

	# select features
	mycorr = corr.to_dict() # {'col1': {'row1': 1, 'row2': 2}, 'col2': {'row1': 0.5, 'row2': 0.75}}
	flt = {}
	for myf1 in mycorr.keys():
		for myf2 in mycorr[myf1].keys():
			if myf2 == myf1:
				continue
			myv = mycorr[myf1][myf2]
			try:
				myv = float(myv)
				if math.isnan(myv):
					continue
				if myv >= float(redu_level):
					flt[myf2] = ""
			except:
				continue
	
	refined_feature = {}
	for myf in features.keys():
		myvalue = features[myf]
		if myf in flt:
			continue
		refined_feature[myf] = myvalue

	return refined_feature


def flt_nonformative_feature (info_level, features):
	"""
	Filter non-informative features
	Input:
		info_level - minimum informative fraction of a given feature
		feature - a dictionary of feature info
	Output: selected features
	"""

	config.logger.info('flt_nonformative_feature')

	# select features
	refined_feature = {}
	for myf in features.keys():
		myvalue = features[myf]
		info = myvalue.split("\t")
		num = len(info)
		mynum = 0
		for myv in info:
			try:
				myv = float(myv)
				if math.isnan(myv) or myv <= 0:
					continue
				mynum = mynum + 1
			except:
				continue
		myper = mynum * 1.0 / num
		if myper <= float(info_level):
			continue
		refined_feature[myf] = myvalue

	return refined_feature


def flt_nonimportant_feature (imp_level, implist, features):
	"""
	Filter non-inportant features
	Input:
		imp_level - minimum prediction importance of a given feature
		implist - a list of file including the prediction importance file for each function
		feature - a dictionary of feature info
	Output: selected features
	"""

	config.logger.info('flt_nonimportant_feature')

	# collect less important features
	hit = {}
	open_list = open(implist, "r")
	for myfile in open_list:
		myfile = myfile.strip()
		if not len(myfile):
			continue
		if not os.path.isfile(myfile):
			config.logger.info("Error! The prediction file doesn't exist: " + myfile)
			continue
		open_file = open(myfile, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			myf = info[0]
			myv = info[-1]
			try:
				myv = float(myv)
				#if math.isnan(myv):
				#	flt[myf] = ""
				#if myv <= float(imp_level):
				#	flt[myf] = ""
				if myv > float(imp_level):
					hit[myf] = ""
			except:
				#flt[myf] = ""
				continue
		open_file.close()
	open_list.close()

	# select features
	refined_feature = {}
	for myf in features.keys():
		myvalue = features[myf]
		if myf in hit:
			refined_feature[myf] = myvalue

	return refined_feature


def write_feature (title, funcs, features, outfile):
	"""
	Write selected features info output file
	Input: functions and selected features
	Output: output into file
	"""

	config.logger.debug('write_feature')

	out_file = open(outfile, 'w')
	out_file.write(title + "\n")
	for item in funcs:
		out_file.write(item + "\n")
	for myf in sorted(features.keys()):
		out_file.write(myf + "\t" + features[myf] + "\n")
	out_file.close()


def main():
	args_value = parse_arguments()

	config.logger.info ("Start select_feature process ......")

	## Collect info ##
	title, funcs, feature = collect_feature (args_value.input)

	## Select feature ##
	if args_value.information:
		try:
			info_level = float(args_value.information)
			feature = flt_nonformative_feature (info_level, feature)
		except ValueError:
			sys.exit("Error! Please use valid value for minimum information fraction by '--information' parameter")
	if args_value.redundancy:
		try:
			cores = int(args_value.core)
			redu_level = float(args_value.redundancy)
			genes = title.split("\t")
			genes = genes[1:len(genes)]
			feature = flt_redundant_feature (redu_level, feature, genes, args_value.method, cores)
		except ValueError:
			sys.exit("Error! Please use valid value for maximum redundant level by '--redundancy' parameter")
	if args_value.importance:
		try:
			imp_level = float(args_value.importance)
			if not args_value.implist:
				sys.exit("Error! Please list file including the prediction importance files by '--implist' parameter")
			feature = flt_nonimportant_feature (imp_level, args_value.implist, feature)
		except ValueError:
			sys.exit("Error! Please use valid value for minimum prediction importance by '--importance' parameter")

	# write select features
	write_feature (title, funcs, feature, args_value.output)

	config.logger.info ("Successfully finished select_feature process!")

if __name__ == '__main__':
	main()
