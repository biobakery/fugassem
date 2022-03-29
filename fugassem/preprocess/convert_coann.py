#!/usr/bin/env python

"""
Convert gene by annotation matrix into gene co-ann matrix

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

try:
	from fugassem import config,utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find FUGAsseM python package." +
		" Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Convert gene by annotation matrix into gene co-ann matrix
"""

def get_args ():
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', "--input",
	                    help = '[REQUIRED] input annotation file',
	                    required = True)
	parser.add_argument('-f', "--format",
	                    help = '[OPTIONAL] format type of the input annotation file [Default: matrix]',
	                    choices = ['matrix', 'vector', 'pair'],
	                    default = "matrix")
	parser.add_argument('-t', "--header",
	                    help = '[OPTIONAL] whether include header or not in the input annotation file [Default: yes]',
						choices= ['yes', 'no'],
						default = 'yes')
	parser.add_argument('-s', "--suffix",
	                    help = '[REQUIRED] suffix of each co-ann',
						required = True)
	parser.add_argument('-o', "--output",
	                    help = 'output refined file',
	                    required = True)
	values = parser.parse_args()

	return values
# get_args


#==============================================================
# collect annotation info
#==============================================================
def collect_annotation_matrix (raw_file):
	#open_file = open(raw_file, "r")
	flag_t = 0
	titles = {}
	anns = {}
	funcs = {}
	families = {}
	for line in utilities.gzip_bzip2_biom_open_readlines(raw_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			for item in info:
				titles[info.index(item)] = item
			continue
		myid = info[0]
		families[myid] = ""
		myindex = 1
		while myindex < len(info):
			myf = titles[myindex]
			myv = info[myindex]
			try:
				myv = float(myv)
			except:
				if myv == "NA":
					myv = float("NaN")
				else:
					config.logger.info("Warning! Invalid value: " + myv)
					myv = float("NaN")
			if not math.isnan(myv): #and myv > 0:
				if not myid in anns:
					anns[myid] = {}
				anns[myid][myf] = myv
				if not myf in funcs:
					funcs[myf] = {}
				funcs[myf][myid] = ""
			myindex = myindex + 1
	# foreach line

	return families, funcs, anns


def collect_annotation_vector (raw_file, header):
	#open_file = open(raw_file, "r")
	anns = {}
	funcs = {}
	families = {}
	titles = {}
	flag_t = 0
	pair_flag = False
	for line in utilities.gzip_bzip2_biom_open_readlines(raw_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			if header == "yes":
				for item in info:
					titles[info.index(item)] = item
				continue
		if len(info) == 3:
			myid1 = info[0]
			myid2 = info[1]
			myv = info[-1]
			families[myid1] = ""
			families[myid2] = ""
			mynew1 = myid1 + "\t" + myid2
			mynew2 = myid2 + "\t" + myid1
			anns[mynew1] = myv
			pair_flag = True
			if not mynew2 in anns:
				anns[mynew2] = myv
			continue
		myid = info[0]
		myf = info[-1]
		families[myid] = ""
		if not myid in anns:
			anns[myid] = {}
		anns[myid][myf] = 1
		if not myf in funcs:
			funcs[myf] = {}
		funcs[myf][myid] = ""
	# foreach line

	return families, funcs, anns, pair_flag


def collect_annotation_pair (raw_file, header):
	coann = {}
	families = {}
	titles = {}
	flag_t = 0
	for line in utilities.gzip_bzip2_biom_open_readlines(raw_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			if header == "yes":
				for item in info:
					titles[info.index(item)] = item
				continue
		myid1 = info[0]
		myid2 = info[1]
		myv = 1
		if len(info) > 2:
			myv = info[-1]
		families[myid1] = ""
		families[myid2] = ""
		mynew1 = myid1 + "\t" + myid2
		mynew2 = myid2 + "\t" + myid1
		coann[mynew1] = myv
		if not mynew2 in coann:
			coann[mynew2] = myv
	# foreach line

	return families, coann


#==============================================================
# get co-annotation info
#==============================================================
def co_annotation (families, funcs, anns, suffix, co_ann, outfile):
	if co_ann:
		coann = anns
	else:
		# get co-annotation info
		coann = {}
		for myid1 in families.keys():
			if myid1 in anns:
				for myf in anns[myid1].keys():
					myv = anns[myid1][myf]
					if myf in funcs:
						for myid2 in funcs[myf].keys():
							if myid1 != myid2:
								mynew = myid1 + "\t" + myid2
								coann[mynew] = myv
								#if not mynew in coann:
								#	coann[mynew] = {}
								#coann[mynew][myf] = myv
						# foreach annotated gene for a given function
				# foreach annotated item for a given gene
			# if gene is annotated
		# foreach gene

	# write co-annotation into file
	open_out = open(outfile, "w")
	title = "ID"
	for myid in sorted(families.keys()):
		title = title + "\t" + myid + "__" + suffix
	open_out.write(title + "\n")
	for myid1 in sorted(families.keys()):
		mystr = myid1
		for myid2 in sorted(families.keys()):
			if myid1 == myid2:
				#myv = len(funcs.keys())
				#myv = 1
				myv = 0
			else:
				mynew = myid1 + "\t" + myid2
				if mynew in coann:
					myv = coann[mynew]
				else:
					myv = 0
			mystr = mystr + "\t" + str(myv)
		open_out.write(mystr + "\n")
	open_out.close()

# co_annotation


#==============================================================
###########  Main processing ############
#==============================================================
def main():	
	### get arguments ###
	values = get_args ()

	config.logger.info ("Start convert_coann process")

	pair_flag = False
	if values.format == "matrix":
		families, funcs, anns = collect_annotation_matrix (values.input)
	if values.format == "vector":
		families, funcs, anns, pair_flag = collect_annotation_vector (values.input, values.header)
	if values.format == "pair":
		families, anns = collect_annotation_pair(values.input, values.header)
		funcs = {}
		pair_flag = True
	co_annotation (families, funcs, anns, values.suffix, pair_flag, values.output)

	config.logger.info ("Successfully finish convert_coann process")

# end: main

if __name__ == '__main__':
	main()
