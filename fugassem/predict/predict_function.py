#!/usr/bin/env python

"""
Predict function based on scoring information

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
		description = "Predict functions based on a scoring system\n",
		prog = "predict_function")
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] known annotations of features\n",
		required = True)
	parser.add_argument(
		"-l", "--list",
		help = "[OPTIONAL] a white list of functions [Default: None]\n",
		default = None)
	parser.add_argument(
		"-b", "--score",
		help = "[REQUIRED] scoring file on co-variation or importance prediction\n",
		required = True)
	parser.add_argument(
		"-m", "--method",
		help = "[REQUIRED] the method for calculating score\n",
		choices = ["max", "mean", "ML"],
		required = True)
	parser.add_argument(
		"-f", "--abs-flag",
		help = "[OPTIONAL] use absolute correlation value or not [Default: no]\n",
		choices = ["yes", "no"],
		default = "no")
	parser.add_argument(
		"-t", "--threshold",
		help = "[OPTIONAL] specify the threshold score for prediction [0-1] [Default: None]\n",
		default = None)
	parser.add_argument(
		"-s", "--start",
		help = "[OPTIONAL] specify the minimum threshold score for testing [0-1] [Default: 0]\n",
		default = 0)
	parser.add_argument(
		"-e", "--end",
		help = "[OPTIONAL] specify the maximum threshold score for testing [0-1] [Default: 1]\n",
		default = 1)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] output file\n",
		required = True)

	return parser.parse_args()


def collect_function_list (func_file):
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
		myid_new = info[0].split(":")
		if myid_new[0] == "GO":
			info[0] = ":".join(myid_new[0:2])
		else:
			info[0] = myid_new[0]
		myfuncs[info[0]] = line

	return myfuncs


def collect_annotation (ann_file, myfuncs):
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
		myid = info[-1]
		myid_new = myid.split(":")
		if (len(myid_new) > 1):
			myid_new = ":".join(myid_new[0:2])
		else:
			myid_new = myid_new[0]
		myid_new2 = myid.split("__")[0]
		myid = re.sub("__", "\t", myid)
		if len(myfuncs.keys()) > 0:
			if not myid_new in myfuncs and not myid_new2 in myfuncs:
					continue
		myf = info[0]
		if not myid in anns:
			anns[myid] = {}
		anns[myid][myf] = ""

	return anns

def collect_score (score_file):
	"""
	Collect scoring information of features
	Input: the file name of scoring file, e.g. co-variation file or importance file
	Output: funcs = {f1->f2: 0,01, ...}
	"""

	config.logger.info ('collect_score')

	score = {}
	features = {}
	titles = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (score_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag == 0:
			flag = 1
			for i in info:
				titles[info.index(i)] = i
			continue
		myid = info[0]
		features[myid] = ""
		myindex = 1
		while myindex < len(info):
			myf = titles[myindex]
			myv = info[myindex]
			if myv == "":
				myv = float("NaN")
			try:
				myv = float(myv)
			except:
				myv = float("NaN")
			try:
				if float(myv) > 1:
					myv = 1
			except:
				myv = float("NaN")
			if myid != myf:
				if not math.isnan(myv):
					score[myid + "\t" + myf] = float(myv)
			myindex = myindex + 1

	return features, score

def assign_hit (anns, score):
	"""
	Assign hit for each feature-function pair
	Input:
		anns - annotation info of features
		features - scores of features
	Output: a data frame for feature-function scoring
	"""

	# select functions
	anns_new = {}
	funcs = {}
	ann_map = {}
	for item in score.keys():
		myid, mya = item.split("\t")
		funcs[mya] = ""
	for item in anns.keys():
		myf, myc = item.split("\t")
		myid_new = myf.split(":")
		if (len(myid_new) > 1):
			myid_new = ":".join(myid_new[0:2])
		else:
			myid_new = myid_new[0]
		myid_new2 = myf.split("__")[0]
		myid2 = "NA"
		if myid_new in funcs:
			myid2 = myid_new
		if myid_new2 in funcs:
			myid2 = myid_new2
		if  myid2 != "NA":
			ann_map[myid2] = item
			if not myid2 in anns_new:
				anns_new[myid2] = {}
			for i in anns[item].keys():
				anns_new[myid2][i] = item

	score_new = {}
	for item in score.keys():
		myid, mya = item.split("\t")
		myv = score[item]
		if mya in anns_new:
			myann = ann_map[mya]
			hit = 0
			if myid in anns_new[mya]:
				hit = 1
			score_new[myid + "\t" + myann] = str(myv) + "\t" + str(hit)

	return score_new


def refine_score (anns, features, cov, score_method, abs_flag):
	"""
	Assign score for each feature-function pair
	Input:
		anns - annotation info of features
		features - list of features
		cov - pairwise covariation info of feature
		score_method - approach for scoring
	Output: a data frame for feature-function scoring
	"""

	config.logger.info ('assign_score')

	scores = {}
	for myf1 in features:
		for mya in anns:
			covs = []
			hit = 0
			if myf1 in anns[mya]:
				hit = 1
			for myf2 in anns[mya]:
				myid = myf1 + "\t" + myf2
				if myid in cov:
					myv = cov[myid]
					if abs_flag == "yes":
						covs.append(abs(myv))
					else:
						covs.append(myv)
				#else:
					#config.logger.info ("Warning! No pairwise covariation info for " + myf1 + " and " + myf2)

			myscore = float("NaN")
			if score_method == "max":
				try:
					if len(covs) > 0:
						covs_new = [float(x) for x in covs if math.isnan(x) == False]
						myscore = max(covs_new)
				except:
					config.logger.info("Warning! Failed in calculating the maximum correlation value: " + myf1)
					myscore = float("NaN")
			if score_method == "mean":
				try:
					if len(covs) > 0:
						covs_new = [float(x) for x in covs if math.isnan(x) == False]
						myscore = sum(covs_new) / len(covs)
				except:
					config.logger.info("Warning! Failed in calculating the mean correlation value: " + myf1)
					myscore = float("NaN")
			if math.isnan(myscore):
				continue
			
			scores[myf1 + "\t" + mya] = str(myscore) + "\t" + str(hit)

	return scores


def perform_prediction (scores, pred_cutoff, start, end):
	"""
	Predict function based on the cutoff of scoring
	Input:
		scores = {f1: func1, ...}
	Output: a data frame for feature-function prediction
	"""

	config.logger.info('perform_prediction')

	prediction = {}
	cutoffs = np.linspace(float(start), float(end), num = 10).tolist()
	for item in scores:
		myscore, myhit = scores[item].split("\t")
		myscore = float(myscore)
		mypred = 0
		mycut = "NA"
		if pred_cutoff:
			mycut = float(pred_cutoff)
			if myscore >= pred_cutoff:
				mypred = 1
			if not item in prediction:
				prediction[item] = []
			prediction[item].append(scores[item] + "\t" + str(mycut) + "\t" + str(mypred))
		else:
			for i in cutoffs:
				mycut = i
				mypred = 0
				if myscore >= float(i):
					mypred = 1
				else:
					mypred = 0
				if not item in prediction:
					prediction[item] = []
				prediction[item].append(scores[item] + "\t" + str(mycut) + "\t" + str(mypred))

	return prediction


def write_prediction (scores, prediction, outfile):
	"""
	Output predicted scores
	Input:
		scores - scoring for each pair of feature and function
		prediction - function prediction based on different cutoff
	Output: output prediction
	"""

	config.logger.debug('write_prediction')
	
	raw_anns = {}
	open_out = open(outfile, "w")
	open_out.write("feature\tfunc\tcategory\tscore\traw_ann\n")
	for i in sorted(scores.keys()):
		if int(scores[i].split("\t")[-1]) == 1:
			term = i.split("\t")
			if not term[1] in raw_anns:
				raw_anns[term[1]] = {}
			raw_anns[term[1]][term[0]] = float(scores[i].split("\t")[0])
		open_out.write(i + "\t" + scores[i] + "\n")
	open_out.close()
	
	"""
	outfile1 = re.sub(".tsv", ".cutoffs.tsv", outfile)
	open_out = open(outfile1, "w")
	open_out.write("feature\tfunc\tcategory\tscore\traw_ann\tcutoff\tnew_ann\n")
	for i in sorted(prediction.keys()):
		for item in prediction[i]:
			open_out.write(i + "\t" + item + "\n")
	open_out.close()
	"""

	# output select functions
	order_func = {}
	order_num = {}
	for i in raw_anns.keys():
		mynum = len(raw_anns[i].keys())
		if not mynum in order_num:
			order_num[mynum] = {}
		order_num[mynum][i] = ""
		myscore = []
		for j in raw_anns[i].keys():
			myscore.append(raw_anns[i][j])
		myscore_new = [x for x in myscore if math.isnan(x) == False]
		#myscore = statistics.median(myscore)
		mys = sum(myscore_new) / len(myscore)
		mys = float(mys)
		if math.isnan(mys):
			continue
		if not mys in order_func:
			order_func[mys] = {}
		order_func[mys][i] = ""

	outfile2 = re.sub(".tsv", ".pred_func_num.txt", outfile)
	open_out = open(outfile2, "w")
	for mynum in sorted(order_num.keys(), key=int, reverse = True):
		for term in sorted(order_num[mynum].keys()):
			open_out.write(term + "\t" + str(mynum) + "\n")
	open_out.close()
	
	outfile2 = re.sub(".tsv", ".pred_func_score.txt", outfile)
	open_out = open(outfile2, "w")
	for mys in sorted(order_func.keys(), key=float, reverse = True):
		for term in sorted(order_func[mys].keys()):
			open_out.write(term + "\t" + str(mys) + "\n")
	open_out.close()


def main():
	args_value = parse_arguments()

	config.logger.info ("Start predict_function process......")

	## get info ##
	if args_value.list:
		myfuncs = collect_function_list (args_value.list)
	else:
		myfuncs = {}
	anns = collect_annotation (args_value.annotation, myfuncs)
	features, score = collect_score (args_value.score)

	## perform prediction ##
	if args_value.method == "ML":
		scores = assign_hit (anns, score)
	else:
		scores = refine_score (anns, features, score, args_value.method, args_value.abs_flag)
	prediction = perform_prediction (scores, args_value.threshold, args_value.start, args_value.end)

	## output info ##
	write_prediction (scores, prediction, args_value.output)

	config.logger.info("Successfully finished predict_function process!")


if __name__ == '__main__':
	main()
