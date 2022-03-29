#!/usr/bin/env python

"""
Split the total MTX abundance table and annotations into stratified taxon-based tables

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
from collections import namedtuple

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
		description = "Split the total MTX abundance table and annotations into stratified taxon-based tables\n",
		prog = "split_taxa")
	parser.add_argument(
		"-i", "--input",
		help = "[REQUIRED] input MTX abundance file\n",
		required = True)
	parser.add_argument(
		"-a", "--annotation",
		help = "[REQUIRED] input function annotation file\n",
		required = True)
	parser.add_argument(
		"-t", "--taxon",
		help = "[REQUIRED] taxonomic level\n",
		choices = ["MSP", "Species", "Genus", "Family", "Order", "Class", "Phylum"],
		required = True)
	parser.add_argument(
		"-p", "--prev",
		help = "[OPTIONAL] the minimum prevalence per feature [ Default: 0.01 ]\n",
		default = 0.01)
	parser.add_argument(
		"-c", "--coverage",
		help = "[REQUIRED] the minimum fraction of annotated genes per taxon [Default: 0.1]\n",
		default = 0.1)
	parser.add_argument(
		"-n", "--number",
		help = "[REQUIRED] the minimum number of total genes per taxon [Default: 500]\n",
		default = 500)
	parser.add_argument(
		"-d", "--abund",
		help = "[OPTIONAL] the minimum detected abundance for each gene [ Default: 0 ] \n",
		default = 0)
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] the output file\n",
		required = True)

	return parser.parse_args()


def collect_ann_func (ann_file):
	"""
	Collect function annotations
	Input:
		ann_file - the file of function annotations
	Output: funcs = {Cluster_XYZ: func1, ...}
	"""

	config.logger.info('collect_ann_func')

	funcs = {}
	for line in utilities.gzip_bzip2_biom_open_readlines (ann_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myid = info[0]
		myf = info[1]
		if not myid in funcs:
			funcs[myid] = {}
		funcs[myid][myf] = ""

	return funcs
# function: collect_ann_func


def check_ann_cov (funcs, abunds):
	"""
	Check the coverage of annotated genes per taxon
	Input:
		funcs - a dictionary of annotations
		abunds - a dictionary of abundance in a given taxon
	Output: the percentage of annotation coverage
	"""

	config.logger.info('check_ann_cov')

	mynum = 0
	for myid in abunds.keys():
		if myid in funcs:
			mynum = mynum + 1
	try:
		cov = float(mynum) / len(abunds.keys())
	except:
		cov = 0

	return cov
# function check_ann_cov

def extract_ann_cov (funcs, genes):
	"""
	Extract subset of function annotation for a give gene list
	Input:
		funcs - a dictionary of annotations
		genes - a list of genes
	Output: a dictionary of annotations for given genes
	"""

	config.logger.info('extract_ann_cov')

	func_sub = {}
	for myid in funcs.keys():
		if myid in genes:
			if not myid in func_sub:
				func_sub[myid] = {}
			for myf in funcs[myid].keys():
				func_sub[myid][myf] = ""

	return func_sub
# function extract_ann_cov


def filter_low_prev_gene (abunds, min_abund, min_prev):
	"""
	Filter genes with prevalence under the cutoff
	Input:
		abunds - a dictionary of abundance info
		min_abund - cutoff for detected abundance per gene per sample
		min_prev - cutoff for minimum prevalence per gene
	Output: abunds_new = {Cluster_XYZ: abunds, ...}
	"""

	config.logger.info ('filter_low_prev_gene')

	abunds_new = {}
	for mykey in abunds:
		myline = abunds[mykey]
		info = myline.split("\t")
		myabund = info[1:len(info)]
		mynum = 0
		for i in myabund:
			if i != "NA" and i != "NaN" and i != "nan":
				myv = float(i)
				if myv > float(min_abund):
					mynum += 1
		myprev = mynum * 1.0 / len(myabund)
		if myprev < float (min_prev):
			continue
		abunds_new[mykey] = myline

	return abunds_new


def write_abunds (nrm_info, abunds, header, outfile):
	"""
	Write refined abundance table
	Input: refined abundance dictionary
	Output: output into file
	"""

	config.logger.info ('write_abunds')

	out_file = open(outfile, 'w')
	out_file.write(header + "\n")
	for myid in sorted(abunds.keys()):
		info = abunds[myid].split("\t")
		if myid in nrm_info:
			info[0] = myid + utilities.c_strat_delim + nrm_info[myid] 
		myline = "\t".join(info)
		out_file.write (myline + "\n")


def write_funcs (funcs, outfile):
	"""
	Output function annotations
	Input: a dictionary of annotations
	Output: output into file
	"""

	config.logger.info ('write_funcs')

	out_file = open(outfile, 'w')
	for myid in sorted(funcs.keys()):
		for myf in sorted(funcs[myid]):
			out_file.write(myid + "\t" + myf + "\n")
	out_file.close()


def write_genelist (genes, outfile):
	"""
	Output gene list
	Input: a dictionary of gene list
	Output: output into file
	"""

	config.logger.info ('write_genelist')

	out_file = open(outfile, 'w')
	for myid in genes:
		out_file.write(myid + "\n")
	out_file.close()


def split_taxa_info (funcs, header, abunds, min_abund, min_prev, min_cov, min_num, taxon_level, out_path):
	"""
	Split abundance and annotation information per taxon
	Input:
		funcs - a dictionary of annotations
		header - header of abundance table
		abunds - a dictionary of abundance table
		min_abund - minimum detected abundance for each gene in each sample
		min_prev - cutoff for minimum prevalence
		min_cov - cutoff for minimum coverage of annotated genes in each taxon
		min_num - cutoff for minimum number of total genes in each taxon
		taxon_level - taxonomic level
		out_path - path of output directory
	Output: batches of abundance and annotation files per taxon
	"""

	config.logger.info ('filter_non_detected_sample')

	# taxa maps
	taxa = {"MSP": "msp__", "Species": "s__", "Genus": "g__", "Family": "f__", "Order": "o__", "Class": "c__", "Phylum": "p__"}
	taxon_level = taxa[taxon_level]

	## collect info
	taxa_info = {}
	nrm_info = {}
	for mykey in abunds.keys():
		myline = abunds[mykey]
		info = myline.split("\t")
		myinfo = info[0].split(utilities.c_strat_delim)
		myid = myinfo[0]
		tmp = myinfo[-1].split(utilities.c_taxon_delim)
		mytaxon = "NA"
		for i in tmp:
			if re.search("^" + taxon_level, i):
				mytaxon = re.sub("^" + taxon_level, "", i)
		if mytaxon == "NA":
			continue
		mynrm_level = tmp[-1]
		mynrm_level = re.sub("[\S]*__", "", mynrm_level)
		mytaxon = utilities.refine_taxon (mytaxon)
		mynrm_level = utilities.refine_taxon (mynrm_level)
		nrm_info[myid] = mynrm_level
		if not mytaxon in taxa_info:
			taxa_info[mytaxon] = {}
		taxa_info[mytaxon][myid] = myline

	## check each taxon and write info into files
	taxa_list = {}
	for mytaxon in taxa_info.keys():
		myabunds = taxa_info[mytaxon]
		pass_flag = True
		if min_prev:
			myabunds = filter_low_prev_gene (myabunds, min_cov, min_prev)
		if min_num:
			try:
				min_num = float(min_num)
				mynum = len(myabunds.keys())
				if mynum >= min_num:
					if pass_flag:
						pass_flag = True
				else:
					pass_flag = False
					continue
			except:
				config.logger.info ("Error! Please provide valid values for minimum number of genes.")
				continue
		else:
			pass_flag = True
		if min_cov:
			try:
				min_cov = float(min_cov)
				mycov = check_ann_cov (funcs, myabunds)
				if mycov >= min_cov:
					if pass_flag:
						pass_flag = True
				else:
					pass_flag = False
					continue
			except:
				config.logger.info ("Error! Please provide valid values for annotation coverage.")
				continue
		else:
			pass_flag = True
		if pass_flag:
			mypath = os.path.join(out_path, mytaxon)
			if not os.path.isdir(mypath):
				os.system("mkdir -p " + mypath)
			myout1 = os.path.join(mypath, mytaxon + config.c_abund_extension)
			myout2 = os.path.join(mypath, mytaxon + config.c_ann_extension)
			myout3 = os.path.join(mypath, mytaxon + config.c_gene_extension)
			genes = sorted(myabunds.keys())
			myfuncs = extract_ann_cov (funcs, genes)
			write_abunds (nrm_info, myabunds, header, myout1)
			write_funcs (myfuncs, myout2)
			write_genelist (genes, myout3)
			taxa_list[mytaxon] = ""
					
	# write file list
	outfile = os.path.join(out_path, config.taxa_abund_list)
	open_out = open(outfile, "w")
	for myid in sorted(taxa_list.keys()):
		open_out.write(myid + "\n")
	open_out.close()

# function: split_taxa_info


def main():
	args_value = parse_arguments()

	config.logger.info ("Start split_taxa process ......")

	## Get info ##
	abunds, header = utilities.read_data_from_file (args_value.input)
	funcs = collect_ann_func (args_value.annotation)
	split_taxa_info(funcs, header, abunds, args_value.abund, args_value.prev, args_value.coverage, args_value.number, args_value.taxon, args_value.output)

	config.logger.info ("Successfully finish split_taxa process!")

if __name__ == '__main__':
	main()
