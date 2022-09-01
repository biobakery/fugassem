#!/usr/bin/env python

"""
Merge prediction results for all taxa

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
Merge prediction results for all taxa
"""

def get_args ():
	parser = argparse.ArgumentParser(description = description)
	parser.add_argument('-i', "--input",
	                    help = '[REQUIRED] input a list file of all outputs',
	                    required = True)
	parser.add_argument('-b', "--basename",
	                    help = '[REQUIRED] basename of outputs',
	                    required = True)
	parser.add_argument('-o', "--output",
	                    help = '[REQUIRED] output file name',
	                    required = True)
	values = parser.parse_args()
	return values
# get_args


#==============================================================
# merge prediction info
#==============================================================
def merge_pred_info (list_file, basename, outfile):
	lists = []
	if not os.path.isfile(list_file):
		config.logger.info("ERROR! This list file doesn't exist: " + list_file)
		return None
	for line in utilities.gzip_bzip2_biom_open_readlines (list_file):
		line = line.strip()
		if not len(line):
			continue
		if not os.path.isfile(line):
			config.logger.info("WARNING! This prediction file doesn't exist: " + line)
			continue
		lists.append(line)

	# merged finalized files
	flag_t = 0
	maps = {}
	open_out = open(outfile, "w")
	for myfile in lists:
		each_f = 0
		mypath = os.path.dirname (myfile)
		mybase = os.path.basename (myfile)
		mytaxon = re.sub(basename + ".", "", mybase)
		mytaxon = re.sub(".finalized_ML.prediction.tsv", "", mytaxon)
		for line in utilities.gzip_bzip2_biom_open_readlines(myfile):
			line = line.strip()
			if not len(line):
				continue
			if each_f == 0:
				each_f = 1
				if flag_t == 0:
					flag_t = 1
					open_out.write("taxon" + "\t" + line + "\n")
				continue
			open_out.write(mytaxon + "\t" + line + "\n")
		mymatch = os.path.join(mypath, "feature_maps.txt")
		if not os.path.isfile(mymatch):
			config.logger.info("WARNING! The feature_maps.txt file doesn't exist: " + mymatch)
			continue
		else:
			for line in utilities.gzip_bzip2_biom_open_readlines(mymatch):
				line = line.strip()
				if not len(line):
					continue
				info = line.split("\t")
				maps[info[0]] = info[1]
	# foreach file
	open_out.close()

	# individual files
	for mytype in maps.keys():
		outfile1 = re.sub("finalized_ML.prediction.tsv", mytype + "_ML.prediction.tsv", outfile)
		flag_t = 0
		open_out = open(outfile1, "w")
		for myfile in lists:
			each_f = 0
			myfile = re.sub("finalized_ML.prediction.tsv", mytype + "_ML.prediction.tsv", myfile)
			if not os.path.isfile(myfile):
				config.logger.info("WARNING! This prediction file doesn't exist: " + myfile)
				continue
			mybase = os.path.basename(myfile)
			mytaxon = re.sub(basename + ".", "", mybase)
			mytaxon = re.sub("." + mytype + "_ML.prediction.tsv", "", mytaxon)
			for line in utilities.gzip_bzip2_biom_open_readlines(myfile):
				line = line.strip()
				if not len(line):
					continue
				if each_f == 0:
					each_f = 1
					if flag_t == 0:
						flag_t = 1
						open_out.write("taxon" + "\t" + line + "\n")
					continue
				open_out.write(mytaxon + "\t" + line + "\n")
		open_out.close()
	# for each type

	# match file
	outfile1 = os.path.join(os.path.dirname(outfile), basename + ".feature_maps.txt")
	utilities.dict_to_file (maps, outfile1)


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start merge_prediction process")

	merge_pred_info (values.input, values.basename, values.output)
	
	config.logger.info ("Successfully finish merge_prediction")

# end: main

if __name__ == '__main__':
	main()
