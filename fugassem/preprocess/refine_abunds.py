#!/usr/bin/env python

"""
Refine abundance table based on prevalence and covariate info (e.g. dna-abundance, taxa-abundance)

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
from collections import namedtuple

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
		description = "Refine abundance table based on prevalence and covariate info\n",
		prog = "refine_abunds")
	parser.add_argument(
		"-i", "--input",
		help = "[REQUIRED] input abundance table used for correlation\n",
		required = True)
	parser.add_argument(
		"-r", "--raw",
		help = "[REQUIRED] raw abundance table without transform and smooth\n",
		required = True)
	parser.add_argument(
		"-c", "--covariate",
		help = "[OPTIONAL] covariate abundance table used for technical-zero filtering [ Default: None ] \n",
		default = None)
	parser.add_argument(
		"-m", "--method",
		help = "[OPTIONAL] pre-filtering approach [ Default: \"none\" ] \n",
		choices = ["lenient", "semi-strict", "strict", "none", "None"],
		default = "none")
	parser.add_argument(
		"-a", "--abund",
		help = "[OPTIONAL] the minimum abundance for each feature [ Default: 0 ] \n",
		default = 0)
	parser.add_argument(
		"-d", "--detected",
		help = "[OPTIONAL] the minimum detected value for each covariate feature [ Default: 0 ] \n",
		default = 0)
	parser.add_argument(
		"-p", "--prev",
		help = "[OPTIONAL] the minimum prevalence per feature [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-q", "--coverage",
		help = "[OPTIONAL] the minimum detection coverage per sample [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def filter_non_detected_sample (header1, header2, header3, abunds, abunds_raw, abunds_cov, min_abund, min_detect, flt_method):
	"""
	Exclude samples that are likely to have technical zeros
	Input:
		header - header of abundance table
		abunds - a dictionary of abundance table used for correlation
		abunds_raw - a dictionary of abundance table without transform and smooth
		abunds_cov - a dictionary of abundance table of covariate
		min_abund - minimum detected abundance for each feature
		min_detect - minimum detected abundance for each covariate
		flt_method - pre-filtering approach
	Output: abunds_new = {Cluster_XYZ: abunds, ...}
	"""

	config.logger.info ('filter_non_detected_sample')

	## collect info
	raw_info = {}
	cov_info = {}
	titles1 = {}
	titles2 = {}
	titles3 = {}
	info = header1.split("\t")
	for item in info:
		titles1[info.index(item)] = item
	info = header2.split("\t")
	for item in info:
		titles2[info.index(item)] = item
	info = header3.split("\t")
	for item in info:
		titles3[info.index(item)] = item
	if len(abunds_cov.keys()) > 0:
		for mykey in abunds_cov:
			myline = abunds_cov[mykey]
			info = myline.split("\t")
			myid = info[0].split(utilities.c_strat_delim)[0]
			myindex = 1
			while myindex < len(info):
				mys = titles3[myindex]
				myv = info[myindex]
				flag = 0
				if myv != "NA" and myv != "NaN" and myv != "nan":
					myv = float(myv)
					if myv > float(min_detect):
						myv = 1
						flag = 1
				if flag == 0:
					myv = 0
				if not myid in cov_info:
					cov_info[myid] = {}
				cov_info[myid][mys] = myv
				myindex = myindex + 1
		# foreach covarite

	for mykey in abunds_raw:
		myline = abunds_raw[mykey]
		info = myline.split("\t")
		myid = info[0]
		myindex = 1
		while myindex < len(info):
			mys = titles2[myindex]
			myv = info[myindex]
			flag = 0
			if myv != "NA" and myv != "NaN" and myv != "nan":
				myv = float(myv)
				if myv > float(min_abund):
					myv = 1
					flag = 1
			if flag == 0:
				myv = 0
			if not myid in raw_info:
				raw_info[myid] = {}
			raw_info[myid][mys] = myv
			myindex = myindex + 1
	# foreach feature

	## filter zeros
	refined_feature = {}
	for mykey in raw_info:
		myid = mykey.split(utilities.c_strat_delim)
		mynum = 0
		mynum_cov = 0
		for mys in raw_info[mykey].keys():
			myv = raw_info[mykey][mys]
			if myv == 1:
				mynum = mynum + 1
			myv_cov = "NA"
			if len(cov_info.keys()) > 0:
				if myid[-1] in cov_info:
					if mys in cov_info[myid[-1]]:
						myv_cov = cov_info[myid[-1]][mys]
				if myv_cov == "NA":
					if myid[0] in cov_info:
						if mys in cov_info[myid[0]]:
							myv_cov = cov_info[myid[0]][mys]
			if myv_cov != "NA":
				if myv_cov == 1:
					mynum_cov = mynum_cov + 1
				if flt_method == "semi-strict":
					if myv == 1 or myv_cov == 1: # filter 0/0
						if not myid[0] in refined_feature:
							refined_feature[myid[0]] = {}
						refined_feature[myid[0]][mys] = ""
				if flt_method == "strict":
					if myv == 1 and myv_cov == 1: # filter 0/0, 1/0, 0/1
						if not myid[0] in refined_feature:
							refined_feature[myid[0]] = {}
						refined_feature[myid[0]][mys] = ""
		if flt_method == "lenient":
			if mynum > 0 and mynum_cov > 0:
				if not myid[0] in refined_feature:
					refined_feature[myid[0]] = {}
				for mys in raw_info[mykey].keys():
					refined_feature[myid[0]][mys] = ""
		if flt_method == "none" or flt_method == "None" or not flt_method:	
			if not myid[0] in refined_feature:
				refined_feature[myid[0]] = {}
			for mys in raw_info[mykey].keys():
				refined_feature[myid[0]][mys] = ""


	## filter abundance table
	abunds_new = {}
	for mykey in abunds:
		myline = abunds[mykey]
		info = myline.split("\t")
		myid = info[0].split(utilities.c_strat_delim)[0]
		mystr = []
		myindex = 1
		while myindex < len(info):
			mys = titles1[myindex]
			myv = info[myindex]
			flag = 0
			if myid in refined_feature:
				if mys in refined_feature[myid]:
					flag = 1
			if flag == 0:
				myv = "NaN"
			mystr.append(myv)
			myindex = myindex + 1
		abunds_new[info[0]] = info[0] + "\t" + "\t".join(mystr)

	return abunds_new


def filter_low_prevalent_feature (abunds, abunds_raw, min_abund, min_prev):
	"""
	Filter features with prevalence less than the cutoff
	Input: a dictionary of abundance info
	Output: abunds_new = {Cluster_XYZ: abunds, ...}
	"""

	config.logger.info ('filter_low_prevalent_feature')

	abunds_new = {}
	for mykey in abunds_raw:
		myline = abunds_raw[mykey]
		info = myline.split("\t")
		myabund = info[1:len(info)]
		mynum = 0
		for i in myabund:
			if i != "NA" and i != "NaN" and i != "nan":
				myv = float(i)
				if myv > float(min_abund):
					mynum += 1
		myprev = mynum / len(myabund)
		if myprev < float (min_prev):
			continue
		myid = mykey.split(utilities.c_strat_delim)
		if mykey in abunds:
			line = abunds[mykey]
			info = line.split("\t")
			abunds_new[info[0]] = line
		else:
			if myid[0] in abunds:
				line = abunds[myid[0]]
				info = line.split("\t")
				abunds_new[info[0]] = line

	return abunds_new


def filter_low_coverage_sample (header1, header2, abunds, abunds_raw, min_cov, min_abund):
	"""
	Filter samples with detection coverage less than the cutoff
	Input: a dictionary of abundance info
	Output: abunds_new = {Cluster_XYZ: abunds, ...}
	"""

	config.logger.info ('filter_low_coverage_sample')

	# select samples
	samples = []
	features = {}
	cov_info = {}
	titles2 = {}
	samples_raw = header2.split("\t")
	for item in samples_raw:
		titles2[samples_raw.index(item)] = item
	for mykey in abunds_raw:
		myline = abunds_raw[mykey]
		info = myline.split("\t")
		myid = info[0]
		features[myid] = ""
		myindex = 1
		while myindex < len(info):
			mys = titles2[myindex]
			myv = info[myindex]
			if myv != "NA" and myv != "NaN" and myv != "nan":
				myv = float(myv)
				if myv > float(min_abund):
					if not mys in cov_info:
						cov_info[mys] = {}
					cov_info[mys][myid] = ""
			myindex = myindex + 1
	# foreach feature
	mytotal = len(features.keys())
	for mys in samples_raw:
		if mys in cov_info:
			myper = len(cov_info[mys].keys()) * 1.0 / mytotal
			if myper >= float(min_cov):
				samples.append(mys)
	header_new = "ID\t" + "\t".join(samples)

	# collect new abunds info
	titles1 = {}
	info = header1.split("\t")
	for item in info:
		titles1[item] = info.index(item)
	abunds_new = {}
	for mykey in abunds:
		myline = abunds[mykey]
		info = myline.split("\t")
		mystr = info[0]
		for mys in samples:
			if mys in titles1:
				myindex = titles1[mys]
				myv = info[myindex]
				mystr = mystr + "\t" + myv
			else:
				mystr = mystr + "\tNaN"
		abunds_new[info[0]] = mystr

	return abunds_new, header_new


def output_abunds (abunds, header, outfile):
	"""
	Output refined abundance table
	Input: refined abundance dictionary
	Output: output into file
	"""

	config.logger.info ('output_abunds')

	out_file = open(outfile, 'w')
	out_file.write(header + "\n")
	for myid in sorted(abunds.keys()):
		myline = abunds[myid]
		info = myline.split("\t")
		myid = info[0].split(utilities.c_strat_delim)
		info[0] = myid[0]
		myline = "\t".join(info)
		out_file.write(myline + "\n")


def main():
	args_value = parse_arguments()
    
	if args_value.prev:
		if args_value.prev == "None" or args_value.prev == "none":
			args_value.prev = None
	if args_value.coverage:
		if args_value.coverage == "None" or args_value.coverage == "none":
			args_value.coverage = None

	config.logger.info ("Start refine_abunds process ......")

	## Get info ##
	abunds, header1 = utilities.read_data_from_file (args_value.input)
	abunds_raw, header2 = utilities.read_data_from_file(args_value.raw)
	if args_value.prev:
		abunds = filter_low_prevalent_feature (abunds, abunds_raw, args_value.abund, args_value.prev)
	if args_value.method:
		if args_value.covariate and args_value.covariate != "None" and args_value.covariate != "none":
			abunds_cov, header3 = utilities.read_data_from_file(args_value.covariate)
		else:
			config.logger.info ("WARNING! No covariate file is provied, so skip filtering genes based on covariate abundance.")
			abunds_cov = {}
			header3 = ""
		abunds = filter_non_detected_sample(header1, header2, header3, abunds, abunds_raw, abunds_cov, args_value.abund, args_value.detected, args_value.method)
	if args_value.coverage:
		if float(args_value.coverage) > 1 or float(args_value.coverage) < 0:
			sys.exit ("Error! Please provide valide detection coverage per sample scaling in [0, 1] via --coverage.")
		abunds, header1 = filter_low_coverage_sample(header1, header2, abunds, abunds_raw, args_value.coverage, args_value.abund)

	## output info ##
	output_abunds (abunds, header1, args_value.output)

	config.logger.info ("Successfully finish refine_abunds process!")

if __name__ == '__main__':
	main()
