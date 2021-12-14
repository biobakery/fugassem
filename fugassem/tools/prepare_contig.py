#!/usr/bin/env python

"""
Prepare contig annotation of each protein family

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
		description = "Prepare contig information of features\n",
		prog = "prepare_contig")
	parser.add_argument(
		"-c", "--clustering",
		help = "[REQUIRED] metawibele finalized clustering file\n",
		required = True)
	parser.add_argument(
		"-m", "--mapping",
		help = "[REQUIRED] assembled gene to contig file\n",
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] output file\n",
		required = True)

	return parser.parse_args()


def collect_gene_map (map_file):
	"""
	Collect gene2contig from the mapping file
	Input: the file name of the function names
	Output: maps = {gene:contig, ...}
	"""
	maps = {}
	flag = 0
	titles = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (map_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag == 0:
			flag = 1
			for item in info:
				titles[item] = info.index(item)
			continue
		gid = info[titles["GID"]]
		cid = info[titles["CID"]]
		maps[gid] = cid

	return maps


def collect_contig_info (clstr_file, maps, outfile):
	"""
	Collect cluster2contig info
	Input: the file name of the clustering, gene2contig map
	Output: output cluster2contig info file
	"""

	outfile1 = re.sub(".tsv", ".simple.tsv", outfile)
	open_out = open(outfile, "w")
	open_out1 = open(outfile1, "w")
	for line in utilities.gzip_bzip2_biom_open_readlines (clstr_file):
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			line = re.sub(">", "", line)
			info = line.split(";")
			gid = info[0]
			cluster = info[1]
			if gid in maps:
				cid = maps[gid]
				open_out1.write(cluster + "\t" + cid + "\n")
				open_out.write(cluster + "\t" + cid + "__Contig" + "\n")
	open_out.close()
	open_out1.close()

# function collect_function_info


def main():
	args_value = parse_arguments()

	config.logger.info ("Start prepare_contig process......")

	## get info ##
	maps = collect_gene_map (args_value.mapping)
	collect_contig_info(args_value.clustering, maps, args_value.output)

	config.logger.info("Successfully finished prepare_contig process!")

if __name__ == '__main__':
	main()
