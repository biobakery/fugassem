#!/usr/bin/env python

"""
Convert gene by DDI annotation file into gene-by-gene co-DDI matrix

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
	from fugassem import config
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
	                    help = '[REQUIRED] input DDI annotation file',
	                    required = True)
	parser.add_argument('-d', "--domain",
	                    help = '[REQUIRED] input domain annotation file',
	                    required = True)
	parser.add_argument('-t', "--header",
	                    help = '[OPTIONAL] whether include header or not in the input annotation file [Default: no]',
						choices= ['yes', 'no'],
						default = 'no')
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
def collect_annotation (raw_file, header):
	open_file = open(raw_file, "r")
	if header == "yes":
		info = open_file.readline().strip().split("\t")
		titles = {}
		for item in info:
			titles[info.index(item)] = item
	anns = {}
	funcs = {}
	families = {}
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
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
	open_file.close()

	return families, funcs, anns


#==============================================================
# get co-DDI info
#==============================================================
def co_ddi (families, domain_ann, ddi_ann, ddi, suffix, outfile):
	# get co-ddi info
	coddi = {}
	for myid1 in families.keys():
		if myid1 in domain_ann:
			myf1 = domain_ann[myid1].keys()
			for myid2 in domain_ann.keys():
				if myid1 == myid2:
					continue
				mynew = myid1 + "\t" + myid2
				myf2 = domain_ann[myid2].keys()
				for i in myf1:
					for j in myf2:
						flag = 0
						myddi = i + ":" + j
						if myddi in ddi:
							flag = 1
						else:
							myddi = j + ":" + i
							if myddi in ddi:
								flag = 1
						if flag == 1:
							if myid1 in ddi_ann:
								if myddi in ddi_ann[myid1]:
									flag += 1
							if myid2 in ddi_ann:
								if myddi in ddi_ann[myid2]:
									flag += 1
							if flag == 3:
								if not mynew in coddi:
									coddi[mynew] = {}
								coddi[mynew][myddi] = 1
					# foreach part2 domain
				# foreach part1 domain
			# foreach gene2
		# if gene is annotated
	# foreach gene1

	# write co-ddi into file
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
				myv = 0
			else:
				mynew = myid1 + "\t" + myid2
				if mynew in coddi:
					#myv = len(coann[mynew].keys())
					myv = 1
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

	config.logger.info ("Start convert_ddi process")

	families1, ddi, ddi_ann = collect_annotation (values.input, values.header)
	families2, domain, domain_ann = collect_annotation (values.domain, values.header)
	co_ddi (families2, domain_ann, ddi_ann, ddi, values.suffix, values.output)

	config.logger.info ("Successfully finish convert_ddi process")

# end: main

if __name__ == '__main__':
	main()
