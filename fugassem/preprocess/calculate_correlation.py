#!/usr/bin/env python

"""
Calculate correlation based on abundance table

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
import scipy.stats
import multiprocessing as mp
from collections import namedtuple
import statistics
import math

try:
	from fugassem.common import corrmp
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
		description = "Calculate pairwise correlation based on abundance information\n",
		prog = "calculate_correlation")
	parser.add_argument(
		"-i", "--input",
		help = "[REQUIRED] input abundance table used for correlation\n",
		required = True)
	parser.add_argument(
		"-m", "--method",
		help = "[OPTIONAL] correlation methods [ Default: Pearson ]\n",
		choices = ["Pearson", "Pearson_adv", "Spearman", "Spearman_adv", "Kendall"],
		default = "Pearson")
	parser.add_argument(
		"-z", "--zero",
		help = "[OPTIONAL] method to handle zeros for correlation calculation [ Default: none ]:\n"
			   "none: do not apply filtering\n"
			   "lenient: ignore features that are zeros across all samples\n"
			   "strict: ignore the samples where one of the paired features is zero\n",
		choices = ["none", "lenient", "strict"],
		default = "none")
	parser.add_argument(
		"-t", "--transform",
		help = "[OPTIONAL] transform into log-scale when doing correlation [ Default: True ]\n",
		default = True)
	parser.add_argument(
		"-n", "--merge",
		help = "[OPTIONAL] merge method for advanced correlation [ Default: arithmetic_mean ]\n",
		choices = ["arithmetic_mean", "harmonic_mean", "median", "max", "min", "sum"],
		default = "arithmetic_mean")
	parser.add_argument(
		"-r", "--normalization",
		help = "[OPTIONAL] normalize correlation matrix [ Default: None ]\n",
		choices = [None, "ranking", "zscoring"],
		default = None)
	parser.add_argument(
		"-q", "--rank-method",
		help = "[OPTIONAL] operation method for ranking [ Default: single ]\n",
		choices = ["single", "table"],
		default = "single")
	parser.add_argument(
		"-c", "--core",
		help = "[OPTIONAL] number of threads, [ Default: 1 ]\n",
		default = 1)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def merging (values, metric):
	myv = float("nan")
	if metric == "arithmetic_mean":
		myv = statistics.mean(values)
	if metric == "harmonic_mean":
		values_n = []
		for i in values:
			if i <= 0:
				continue
			values_n.append(i)
		if len(values_n) > 0:
			myv = statistics.harmonic_mean(values_n)
	if metric == "max":
		myv = max(values)
	if metric == "min":
		myv = max(values)
	if metric == "median":
		myv = statistics.median(values)
	if metric == "sum":
		myv = sum(values)
	return myv


def merge_corr (corr1, corr2, merge_method):
	"""
	Merge pairwise correlation between features from different approaches
	Input:
		corr1 - correlation matrix based one method
		corr2 - correlation matrix based another method
		merge_method - method for merging
	Output: correlation matrix
	"""
	Xcolumns = corr1.columns.values.tolist()
	Xrows = Xcolumns
	corr = pd.DataFrame(index = Xrows, columns = Xcolumns)
	for i in Xcolumns:
		for j in Xrows:
			myv = float("nan")
			myv1 = float(corr1[i][j])
			try:
				myv2 = float(corr2[i][j])
			except:
				myv2 = float("nan")
			values = []
			if not math.isnan(myv1):
				values.append(myv1)
			if not math.isnan(myv2):
				values.append(myv2)
			myv = merging (values, merge_method)
			corr[i][j] = myv
	
	return corr


def zscoring (df):
	"""
	Calcualte z-scores of a correlation matrix
	Input:
		df - correlation matrix
	Output: correlation matrix
	"""

	corr = scipy.stats.zscore(df, nan_policy='omit')

	return corr


def ranking (df, rank_method):
	"""
	Rank correlation matrix
	Input:
		df - correlation matrix
		rank_method - method for ranking
	Output: correlation matrix
	"""

	Xcolumns = df.columns.values.tolist()
	if rank_method == "single":
		corr = df.rank(axis=0, method='average', numeric_only=True, na_option='keep', ascending=True, pct=True)
		corr = corr / corr.max()
	if rank_method == "table":
		numpy_matrix = df.to_numpy()
		flat_matrix = numpy_matrix.flatten()
		flat_series = pd.Series(flat_matrix)
		flat_series = flat_series.rank(axis=0, method='average', numeric_only=True, na_option='keep', ascending=True, pct=True)
		flat_series = flat_series / flat_series.max()
		ranks = flat_series.values
		rank_matrix = ranks.reshape(numpy_matrix.shape)
		corr = pd.DataFrame(rank_matrix, index=Xcolumns, columns=Xcolumns)

	return corr


def pairwise_correlation (abunds, samples, corr_method, merge_method, zero_method, trans_method, cores):
	"""
	Calculate pairwise correlation between features
	Input:
		abunds - a dictionary of abundance info
		samples - sample arrays
		corr_method - correlation method
		trans_method - whether to transform into log-scale
		cores - number of threads
	Output: correlation matrix
	"""

	config.logger.info('pairwise_correlation')

	# build data frame
	evidence_row = sorted(abunds.keys())
	evidence_table_row = namedtuple("evidence_table_row", evidence_row, rename = False)
	abunds_table = pd.DataFrame(index = samples, columns = evidence_table_row._fields)

	for item in evidence_row:
		myvalue = abunds[item]
		abunds_table[item] = myvalue

	# correlation
	if corr_method == "Pearson":
		corr, p_value, se = NaNCorrMp.calculate_with_pvalue_se (abunds_table, n_jobs=cores, zero=zero_method, transform=trans_method)
	elif corr_method == "Pearson_adv":
		JC_matrix = corrmp.jaccard (abunds_table)
		corr, p_value, se = NaNCorrMp.calculate_with_pvalue_se (abunds_table, n_jobs=cores, zero=zero_method, transform=trans_method)
		if isinstance(corr, pd.DataFrame):
			corr = corr.replace({None: np.nan})
			corr = corr.replace({"": np.nan})
		corr = merge_corr (JC_matrix, corr, merge_method)
	elif corr_method == "Spearman":
		corr, p_value = NaNCorrMp.calculate_with_pvalue (abunds_table, n_jobs=cores, zero=zero_method, transform=trans_method, corr="spearman")
		se = None
	elif corr_method == "Spearman_adv":
		JC_matrix = corrmp.jaccard (abunds_table)
		corr, p_value = NaNCorrMp.calculate_with_pvalue (abunds_table, n_jobs=cores, zero=zero_method, transform=trans_method, corr="spearman")
		if isinstance(corr, pd.DataFrame):
			corr = corr.replace({None: np.nan})
			corr = corr.replace({"": np.nan})
		corr = merge_corr (JC_matrix, corr, merge_method)
		se = None
	else:
		corr_method = corr_method.lower()
		corr = abunds_table.corr (method = corr_method)
		p_value = None
		se = None
	
	if isinstance(corr, pd.DataFrame):
		corr = corr.replace({None: np.nan})
		corr = corr.replace({'': np.nan})
	if isinstance(p_value, pd.DataFrame):
		p_value = p_value.replace({None: np.nan})
		p_value = p_value.replace({'': np.nan})
	if isinstance(se, pd.DataFrame):
		se = se.replace({None: np.nan})
		se = se.replace({'': np.nan})

	return corr, p_value, se


def output_correlation (corr, p_value, se, outfile):
	"""
	Print correlation matrix
	Input: corr matrix and p-value matrix
	Output: output into file
	"""

	config.logger.info('output_correlation')

	keys = corr.columns.values.tolist()
	out_file = open(outfile, 'w')
	out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
	out_file.close()
	corr.to_csv(outfile, mode = 'a', sep = '\t', header = False, na_rep = float("NaN"))

	outfile1 = re.sub(".tsv", ".pvalues.tsv", outfile)
	if p_value is not None:
		keys = p_value.columns.values.tolist()
		out_file = open(outfile1, 'w')
		out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
		out_file.close()
		p_value.to_csv(outfile1, mode = 'a', sep = '\t', header = False, na_rep = float("NaN"))

	outfile2 = re.sub(".tsv", ".stderr.tsv", outfile)
	if se is not None:
		keys = se.columns.values.tolist()
		out_file = open(outfile2, 'w')
		out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
		out_file.close()
		se.to_csv(outfile2, mode = 'a', sep = '\t', header = False, na_rep = float("NaN"))


def main():
	args_value = parse_arguments()

	config.logger.info ("Start calculate_correlation process ......")

	## Get info ##
	abunds_raw, header = utilities.read_data_from_file (args_value.input, False, "yes", "yes")
	samples = header.split("\t")[1:]
	abunds = {}
	for myid in abunds_raw.keys():
		info = abunds_raw[myid].split("\t")
		myabund = info[1:len(info)]
		myabund = [float(x) for x in myabund]
		abunds[myid] = myabund

	## do correlation ##
	try:
		cores = int(args_value.core)
	except ValueError:
		config.logger.info ("Error! Please use valid number of threads!")
		cores = 1
	if args_value.zero == "none":
		args_value.transform = False
	if args_value.transform == "True":
		args_value.transform = True
	if args_value.transform == "False":
		args_value.transform = False
	corr, p_value, se = pairwise_correlation (abunds, samples, args_value.method, args_value.merge, args_value.zero, args_value.transform, cores)

	## rank correlation ##
	if args_value.normalization:
		if args_value.normalization == "zscoring":
			corr = zscoring (corr)
		if args_value.normalization == "ranking":
			corr = ranking (corr, args_value.rank_method)

	## write outputs
	output_correlation (corr, p_value, se, args_value.output)


	config.logger.info ("Successfully finished calculate_correlation process!")

if __name__ == '__main__':
	main()
