#!/usr/bin/env python

"""
Impute NaN in data matrix

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
import re
import argparse
import logging
import math
import pandas as pd

try:
	from fugassem import config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find FUGAsseM python package." +
		" Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Impute NaN in data matrix
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--input",
	                    help = '[REQUIRED] input raw matrix file',
	                    required = True)
	parser.add_argument('-m', "--method",
	                    help = '[REQUIRED] methods for imputation',
	                    choices = ["zero", "mean", "none"],
						default = "zero")
	parser.add_argument('-n', "--nonna",
	                    help = '[OPTIONAL] minimum fraction of non-NA per feature [ Default: 0]',
						default = 0)
	parser.add_argument('-o', "--output",
	                    help = 'output refined file',
	                    required = True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# Impute NAs
#==============================================================
def flt_nan (raw_file, non_na):
	coexp = pd.read_csv(raw_file, sep="\t", header=0, index_col=0)
	total = len(coexp.columns.values.tolist())
	non_na = int(total * float(non_na))
	coexp = coexp.dropna(axis=1, thresh=non_na)

	return coexp


def imputation (method, coexp, outfile):
	keys = coexp.columns.values.tolist()
	open_out = open(outfile, "w")
	open_out.write("ID" + "\t" + "\t".join(keys) + "\n")
	for index, row in coexp.iterrows():
		rows = [str(x) for x in row.tolist()]
		line = row.name + "\t" + "\t".join(rows)
		info = line.split("\t")
		values = []
		myindex = 1
		while myindex < len(info):
			myv = info[myindex]
			try:
				myv = float(myv)
			except:
				if myv == "NA":
					myv = float("NaN")
				else:
					config.logger.info("Warning! Invalid value: " + myv)
					myv = float("NaN")
			if math.isnan(myv):
				if method == "zero":
					myv = 0
			else:
				values.append(myv)
			info[myindex] = str(myv)
			myindex = myindex + 1
		if method == "mean":
			mymean = sum(values) / len(info)
			myindex = 1
			while myindex < len(info):
				if math.isnan(float(info[myindex])):
					info[myindex] = str(mymean)
				myindex = myindex + 1
		line = "\t".join(info)
		open_out.write(line + "\n")
	# foreach line
	open_out.close()

# collect_info


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start impute_nan process")

	coexp = flt_nan (values.input, values.nonna)
	imputation (values.method, coexp, values.output)
	
	config.logger.info ("Successfully finish impute_nan process")

# end: main

if __name__ == '__main__':
	main()
