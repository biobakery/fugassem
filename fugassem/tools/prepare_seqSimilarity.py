#!/usr/bin/env python

"""
Prepare sequence similarity information as evidence for prediction

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
import logging
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag


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
		description = "Prepare sequence similarity information as evidence for prediction\n",
		prog = "prepare_seqSimilarity")
	parser.add_argument(
		"-a", "--family",
		help = "[REQUIRED] protein family file\n",
		required = True)
	parser.add_argument(
		"-t", "--type",
		help = "[OPTIONAL] type of sequence similarity info [Default: homology]\n",
		default = "homology")
	parser.add_argument(
		"-d", "--vector",
		help = "[OPTIONAL] type of vector encoding [Default: count]\n",
		choices = ["count", "jaccard", "cosine"],
		default = "count")
	parser.add_argument(
		"-m", "--method",
		help = "[OPTIONAL] function type [Default: GO]\n",
		default = "GO")
	parser.add_argument(
		"-s", "--similarity",
		help = "[REQUIRED] sequence similarity info file\n",
		required = True)
	parser.add_argument(
		"-f", "--function",
		help = "[REQUIRED] original function annotation file\n",
		required = True)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] name of output file\n",
		required = True)

	return parser.parse_args()


def collect_families (family_file):
	"""
	Collect protein families
	Input: the file name of protein family
	Output: families = {cluster1, cluster2, ...}
	"""
	
	config.logger.info ("collect_families process......")
	
	families = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (family_file):
		line = line.strip()
		if not len(line):
			continue
		flag = flag + 1
		if flag == 1:
			continue
		info = line.split("\t")
		families[info[0]] = ""

	return families


# define Jaccard Similarity function
def Jaccard (vec1, vec2):
	intersection = len(list(set(vec1).intersection(vec2)))
	union = (len(vec1) + len(vec2)) - intersection
	return float(intersection) / union

# define Cosine Similarity function
def VectorSize (vec):
	return math.sqrt(sum(math.pow(v,2) for v in vec))

def InnerProduct (vec1, vec2):
	return sum(v1*v2 for v1,v2 in zip(vec1,vec2))

def Cosine (vec1, vec2):
	return InnerProduct(vec1, vec2) / (VectorSize(vec1) * VectorSize(vec2))

# calculate similarity
def collect_similarity_unit (simi_file, dist_type):
	"""
	Collect sequence similarity based on Jaccard index
	Input:
		simi_file: the file name of sequence similarity
		dist_type: type of similarity
	Output: pair = {{c1:c2}, ...}
	"""
	config.logger.info ("collect_similarity_unit process......")
	
	units = {}
	maps = {}
	pairs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (simi_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myu = info[1]
		if not myid in maps:
			maps[myid] = {}
		maps[myid][myu] = ""
		if not myu in units:
			units[myu] = {}
		units[myu][myid] = ""

	if dist_type == "count":
		for myid1 in maps.keys():
			for myu in maps[myid1].keys():
				if myu in units:
					for myid2 in units[myu].keys():
						if myid1 != myid2:
							if not myid1 in pairs:
								pairs[myid1] = {}
							if not myid2 in pairs[myid1]:
								pairs[myid1][myid2] = 1
							else:
								pairs[myid1][myid2] += 1

	if dist_type == "jaccard":
		for myid1 in maps.keys():
			for myu in maps[myid1].keys():
				if myu in units:
					for myid2 in units[myu].keys():
						if myid1 == myid2:
							continue
						if not myid2 in maps:
							continue
						if myid1 in pairs:
							if myid2 in pairs[myid1]:
								continue
						if myid2 in pairs:
							if myid1 in pairs[myid2]:
								continue
						vec1 = list(maps[myid1].keys())
						vec2 = list(maps[myid2].keys())
						try:
							sim = Jaccard (vec1, vec2)
							if not myid1 in pairs:
								pairs[myid1] = {}
							if not myid2 in pairs[myid1]:
								pairs[myid1][myid2] = sim
							if not myid2 in pairs:
								pairs[myid2] = {}
							if not myid1 in pairs[myid2]:
								pairs[myid2][myid1] = sim
						except ValueError:
							config.logger ("Failed in calculating Jaccard similarity for " + myid1 + "\t" + myid2)

	if dist_type == "cosine":
		for myid1 in maps.keys():
			for myu in maps[myid1].keys():
				if myu in units:
					for myid2 in units[myu].keys():
						if myid1 == myid2:
							continue
						if not myid2 in maps:
							continue
						if myid1 in pairs:
							if myid2 in pairs[myid1]:
								continue
						if myid2 in pairs:
							if myid1 in pairs[myid2]:
								continue
						vec1 = []
						vec2 = []
						myunits = set(list(maps[myid1].keys()) + list(maps[myid2].keys()))
						for i in myunits:
							if i in maps[myid1]:
								vec1.append(1)
							else:
								vec1.append(0)
							if i in maps[myid2]:
								vec2.append(1)
							else:
								vec2.append(0)
						try:
							sim = Cosine (vec1, vec2)
							if not myid1 in pairs:
								pairs[myid1] = {}
							if not myid2 in pairs[myid1]:
								pairs[myid1][myid2] = sim
							if not myid2 in pairs:
								pairs[myid2] = {}
							if not myid1 in pairs[myid2]:
								pairs[myid2][myid1] = sim
						except ValueError:
							config.logger("Failed in calculating Cosine similarity for " + myid1 + "\t" + myid2)

	maps = {}
	units = {}

	return pairs


def collect_function (func_file, func_type):
	"""
	Collect function maps
	Input:
		the file name of function annotation
	Output: func = {{f1: c1, c2}, ...}
	"""
	
	config.logger.info ("collect_function process......")
	
	go_types = ["GO", "BP", "MF", "CC"]
	if func_type in go_types:
		try:
			godag = get_godag(config.go_obo)
		except:
			godag = None
			config.logger.info("Error! Failed to fetch GO file.")
	funcs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (func_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myf = info[1]
		if func_type in go_types and godag is not None:
			gosubdag = GoSubDag([myf], godag)
			try:
				go_id, go_term = max(gosubdag.go2obj.items(), key=lambda t: t[1].depth)
				go_ids_chosen = list(go_term.get_all_parents())
				go_ids_chosen.append(myf)
			except:
				continue
		else:
			go_ids_chosen = [myf] 
		if not myid in funcs:
			funcs[myid] = {}
		for i in go_ids_chosen:
			funcs[myid][i] = ""

	return funcs


def build_function_matrix (families, pairs, funcs, dist_type, outfile):
	"""
	Build function association matrix based on sequence similarity
	Input: families, pairs, funcs, dist_type
	Output: output function association info file
	"""
	
	config.logger.info ("build_function_matrix process......")

	anns = {}
	myfuncs = {}
	# requiring at least two annotated members in the co-unit if annotated member is labeled as non-zero value
	for myid in families.keys():
		#if myid in funcs:
		#	if not myid in anns:
		#		anns[myid] = {}
		#	for myf in funcs[myid].keys():
		#		myfuncs[myf] = ""
		#		anns[myid][myf] = 1
		if dist_type == "count":
			if myid in pairs:
				for myid2 in pairs[myid].keys():
					if myid == myid2:
						continue
					myv = float(pairs[myid][myid2])
					if myid2 in funcs:
						for myf in funcs[myid2].keys():
							if not myid in anns:
								anns[myid] = {}
							if not myf in anns[myid]:
								myfuncs[myf] = ""
								anns[myid][myf] = myv
							if float(anns[myid][myf]) < float(myv):
								anns[myid][myf] = myv
				# foreach pair
			# if exits pair
		else:
			if myid in pairs:
				sim = {}
				# collect similarity
				for myid2 in pairs[myid].keys():
					if myid == myid2:
						continue
					myv = float(pairs[myid][myid2])
					if myid2 in funcs:
						for myf in funcs[myid2].keys():
							if not myf in sim:
								sim[myf] = []
							sim[myf].append(myv)
				# average similarity per function
				if len(sim.keys()) > 0:
					for myf in sim.keys():
						myfuncs[myf] = ""
						myv = sum(sim[myf]) * 1.0 / len(sim[myf])
						if not myid in anns:
							anns[myid] = {}
						if not myf in anns[myid]:
							anns[myid][myf] = myv
	# foreach protein family

	# write info into file
	mytitle = "ID"
	for myf in sorted(myfuncs.keys()):
		mytitle = mytitle + "\t" + myf + "__seqSimilarity"
	open_out = open(outfile, "w")
	open_out.write(mytitle + "\n")
	for myid in sorted(families.keys()):
		mystr = myid
		for myf in sorted(myfuncs.keys()):
			if myid in anns:
				if myf in anns[myid]:
					mystr = mystr + "\t" + str(anns[myid][myf])
				else:
					mystr = mystr + "\t0"
			else:
				mystr = mystr + "\t0"
		open_out.write(mystr + "\n")
	open_out.close()

# function build_function_matrix


def main():
	args_value = parse_arguments()

	config.logger.info ("Start prepare_seqSimilarity process......")

	## get info ##
	families = collect_families (args_value.family)
	funcs = collect_function (args_value.function, args_value.method)

	#units = ["pfam", "uniref50", "contig", "DDI", "homology"]
	paris = {}
	pairs = collect_similarity_unit (args_value.similarity, args_value.vector)

	build_function_matrix(families, pairs, funcs, args_value.vector, args_value.output)

	config.logger.info("Successfully finished prepare_seqSimilarity process!")

if __name__ == '__main__':
	main()
