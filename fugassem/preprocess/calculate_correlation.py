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
import multiprocessing as mp
from collections import namedtuple


try:
	from fugassem.common.nancorrmp import NaNCorrMp
except:
	sys.exit ("nancorrmp is not install!")
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
		help = "[OPTIONAL] correlation methods, [ Default: Pearson ]\n",
		choices = ["Pearson", "Spearman", "Kendall", "Pearson_SE"],
		default = "Pearson")
	parser.add_argument(
		"-c", "--core",
		help = "[OPTIONAL] number of threads, [ Default: 1 ]\n",
		default = 1)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def pairwise_correlation (abunds, samples, corr_method, cores):
	"""
	Calculate pairwise correlation between features
	Input:
		abunds - a dictionary of abundance info
		samples - sample arrays
		corr_method - correlation method
		cores - number of threads
	Output: correlation matrix
	"""

	config.logger.info('pairwise_correlation')

	# build data frame
	evidence_row = sorted(abunds.keys())
	try:
		evidence_table_row = namedtuple("evidence_table_row", evidence_row, verbose = False, rename = False)
	except:
		evidence_table_row = namedtuple("evidence_table_row", evidence_row, rename = False)
	abunds_table = pd.DataFrame(index = samples, columns = evidence_table_row._fields)

	for item in evidence_row:
		myvalue = abunds[item]
		abunds_table[item] = myvalue

	# correlation
	if corr_method == "Pearson":
		corr, p_value = NaNCorrMp.calculate_with_p_value (abunds_table, n_jobs = cores)
		se = None
	else:
		if corr_method == "Pearson_SE":
			corr, p_value, se = NaNCorrMp.calculate_with_pvalue_se (abunds_table, n_jobs=cores)
		else:
			corr_method = corr_method.lower()
			corr = abunds_table.corr (method = corr_method)
			p_value = None
			se = None

	return corr, p_value, se


def print_correlation (corr, p_value, se, outfile):
	"""
	Print correlation matrix
	Input: corr matrix and p-value matrix
	Output: output into file
	"""

	config.logger.debug('print_correlation')

	keys = corr.columns.values.tolist()
	out_file = open(outfile, 'w')
	out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
	out_file.close()
	corr.to_csv(outfile, mode = 'a', sep = '\t', header = False)

	outfile1 = re.sub(".tsv", ".pvalues.tsv", outfile)
	if p_value is not None:
		keys = p_value.columns.values.tolist()
		out_file = open(outfile1, 'w')
		out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
		out_file.close()
		p_value.to_csv(outfile1, mode = 'a', sep = '\t', header = False)

	outfile2 = re.sub(".tsv", ".stderr.tsv", outfile)
	if se is not None:
		keys = se.columns.values.tolist()
		out_file = open(outfile2, 'w')
		out_file.write("ID" + "\t" + "\t".join(keys) + "\n")
		out_file.close()
		se.to_csv(outfile2, mode = 'a', sep = '\t', header = False)

def main():
	args_value = parse_arguments()

	config.logger.info ("Start calculate_correlation process ......")

	## Get info ##
	abunds_raw, header = utilities.read_data_from_file (args_value.input)
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
	corr, p_value, se = pairwise_correlation (abunds, samples, args_value.method, cores)
	print_correlation (corr, p_value, se, args_value.output)

	config.logger.info ("Successfully finished calculate_correlation process!")

if __name__ == '__main__':
	main()
