#!/usr/bin/env python

"""
Collect ML prediction results for each function

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
import math
import numpy as np
import logging


try:
	from fugassem import utilities
	from fugassem import config
except:
	sys.exit ("fugassem is not install!")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Collect ML info for each function
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-p', "--path",
	                    help = '[REQUIRED] directary of importance result files',
	                    required = True)
	parser.add_argument('-e', "--extension",
	                    help = '[REQUIRED] extension of ML result file',
	                    choices = [".importance.tsv", ".xval.roc.tsv"],
						required = True)
	parser.add_argument('-t', "--type",
	                    help = '[REQUIRED] type of results',
	                    choices = ["importance", "xval"],
						required = True)
	parser.add_argument('-o', "--output",
	                    help = '[REQUIRED] output file name',
	                    required = True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# collect ML info
#==============================================================
def collect_ML_info (type, extension, path):
	results = {}
	funcs = {}
	labels = {}
	filelist = utilities.find_files(path, extension, None)
	for myfile in filelist:
		if not os.path.isfile(myfile):
			config.logger.info("Warning! This file doesn't exist: " + myfile)
			continue
		myfunc = os.path.basename(myfile)
		mym = re.search("([\S]+)" + extension, myfunc)
		myfunc = mym.group(1)
		myfunc = re.sub("_", ":", myfunc)
		funcs[myfunc] = ""
		open_file = open(myfile, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			if type == "importance":
				myid = info[0]
				myscore = info[-1]
				if not myid in results:
					results[myid] = {}
				results[myid][myfunc] = myscore
			if type == "xval":
				myid = info[0]
				sdicts = []
				d = {}
				for k, v in zip( info[2::2], info[3::2] ):
					d[k] = float( v )
					sdicts.append( d )
				for k in d.keys():
					myv = d[k]
					if not k in results:
						results[k] = {}
					if not myid in results[k]:
						results[k][myid] = {}
					results[k][myid][myfunc] = myv
		open_file.close()
	# foreach function

	return results, funcs


# ==============================================================
# write ML info
# ==============================================================
def write_ML_results (type, results, funcs, outfile):
	title = "feature" + "\t" + "\t".join(sorted(funcs.keys()))
	if type == "importance":
		open_out = open(outfile, "w")
		open_out.write(title + "\n")
		for myid in results.keys():
			mystr = myid
			for myf in sorted(funcs.keys()):
				myv = "NaN"
				if myf in results[myid]:
					myv = results[myid][myf]
			mystr = mystr + "\t" + str(myv)
			open_out.write(mystr + "\n")
		# foreach line
		open_out.close()
	if type == "xval":
		for myc in results.keys():
			outfile1 = re.sub(".tsv", "." + myc + ".tsv", outfile)
			open_out = open(outfile1, "w")
			open_out.write(title + "\n")
			for myid in results[myc].keys():
				mystr = myid
				for myf in sorted(funcs.keys()):
					myv = "NaN"
					if myf in results[myc][myid]:
						myv = results[myc][myid][myf]
					mystr = mystr + "\t" + str(myv)
				open_out.write(mystr + "\n")
			open_out.close()

# write_ML_results


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start collect_ML_results process")

	results, funcs = collect_ML_info (values.type, values.extension, values.path)
	write_ML_results (values.type, results, funcs, values.output)
	
	config.logger.info ("Successfully finish collect_ML_results")

# end: main

if __name__ == '__main__':
	main()
