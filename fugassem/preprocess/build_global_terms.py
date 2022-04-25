#!/usr/bin/env python

"""
Build global terms used for prediction for a community

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
from collections import Counter
import glob

try:
	from fugassem import utilities
	from fugassem import config
	from fugassem.common import geneontology
	from fugassem.common.utils import path2name, warn, qw
	from fugassem.common.dictation import polymap, col2dict
except:
	sys.exit ("fugassem is not installed!")


def parse_arguments():
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "Build global terms used for prediction for a community\n",
		prog = "build_global_terms")
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
		default = "Species")
	parser.add_argument(
		"-p", "--prev",
		help = "[OPTIONAL] the minimum prevalence per feature [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-c", "--coverage",
		help = "[REQUIRED] the minimum fraction of annotated genes per taxon [Default: None ]\n",
		default = None)
	parser.add_argument(
		"-n", "--number",
		help = "[REQUIRED] the minimum number of total genes per taxon [ Default: None ]\n",
		default = None)
	parser.add_argument(
		"-d", "--abund",
		help = "[OPTIONAL] the minimum detected abundance for each gene [ Default: 0 ] \n",
		default = 0)
	parser.add_argument(
		"-b", "--obo",
	    help = "[REQUIRED] GO-provided obo file\n",
	    required = True)
	parser.add_argument(
		"-m", "--namespace",
	    default = None,
	    choices = ["BP", "MF", "CC"],
	    nargs = "+",
	    metavar = "<BP/MF/CC>",
	    help = "trimming option: only terms in this namespace")
	parser.add_argument(
		"-q", "--informative",
	    default = None,
	    metavar = "<number OR fraction of genes>",
	    help = "trimming option: only terms that are informative at a given level")
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
		myf = info[1].split("__")[0]
		myf = myf.split(":")
		if len(myf) > 1:
			myf = ":".join(myf[0:2])
		else:
			myf = myf[0]
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


def write_funcs_combined (funcs, terms, outfile):
	"""
	Output function annotations
	Input: a dictionary of annotations
	Output: output into file
	"""

	config.logger.info ('write_funcs_combined')

	myterms = {}
	for i in terms.keys():
		for j in terms[i].keys():
			myterms[j] = ""
	out_file = open(outfile, 'w')
	for myid in sorted(funcs.keys()):
		for myf in sorted(funcs[myid]):
			if myf in myterms:
				out_file.write(myid + "\t" + myf + "\n")
	out_file.close()



def write_funcs_final (funcs, outfile):
	"""
	Output function annotations
	Input: a dictionary of annotations
	Output: output into file
	"""

	config.logger.info ('write_funcs_final')

	outfile2 = re.sub(".tsv", ".simple.tsv", outfile)
	out_file1 = open(outfile, 'w')
	out_file2 = open(outfile2, 'w')
	for myid in sorted(funcs.keys()):
		for myf in sorted(funcs[myid]):
			out_file1.write(myid + "\t" + myf + "\n")
			myf = myf.split(":")
			if len(myf) > 1:
				myf = ":".join(myf[0:2])
			else:
				myf = myf[0]
			myf = myf.split("__")[0]
			out_file2.write(myid + "\t" + myf + "\n")
	out_file1.close()
	out_file2.close()


def collect_taxa_ann (funcs, header, abunds, min_abund, min_prev, min_cov, min_num, taxon_level, out_path, namespace):
	"""
	Collect function annotations per taxon
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

	config.logger.info('collect_taxa_ann')

	# taxa maps
	taxa = {"MSP": "msp__", "Species": "s__", "Genus": "g__", "Family": "f__", "Order": "o__", "Class": "c__",
	        "Phylum": "p__"}
	taxon_level = taxa[taxon_level]

	## collect info
	taxa_info = {}
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
		mytaxon = utilities.refine_taxon(mytaxon)
		mynrm_level = utilities.refine_taxon(mynrm_level)
		if not mytaxon in taxa_info:
			taxa_info[mytaxon] = {}
		taxa_info[mytaxon][myid] = myline

	## check each taxon and collect annotation
	map_list = []
	for mytaxon in taxa_info.keys():
		myabunds = taxa_info[mytaxon]
		pass_flag = True
		if min_prev:
			myabunds = filter_low_prev_gene(myabunds, min_cov, min_prev)
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
				config.logger.info("Error! Please provide valid values for minimum number of genes.")
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
				config.logger.info("Error! Please provide valid values for annotation coverage.")
				continue
		else:
			pass_flag = True
		if pass_flag:
			mymap = os.path.join(out_path, mytaxon + "." + namespace + ".genemap.tsv")
			map_list.append(mymap)
			genes = sorted (myabunds.keys())
			myfuncs = extract_ann_cov (funcs, genes)
			write_funcs (myfuncs, mymap)

	return map_list

# function: collect_taxa_ann


def select_bug_term (obo, map_file, informative, namespace, terms, outfile, desc):
	# attach genes
	if map_file:
		warn("Process file:", map_file)
		mapping = geneontology.polymap(map_file)
		obo.attach_genes (mapping)
		warn("# of attached genes:", len(obo.attached_genes))

		# informative cut
		if informative:
			threshold = float(informative)
			if threshold < 1:
				warn("Intepretting informative cutoff as fraction of annotated genes")
				threshold *= len(obo.attached_genes)
			threshold = int(threshold)
			obo.set_informative(threshold)
			for term in obo.iter_terms():
				if not term.is_informative:
					term.is_acceptable = False
		# namespace cut
		if namespace:
			for term in obo.iter_terms():
				if term.namespace_short not in namespace:
					term.is_acceptable = False

		# output the new polymap
		fh = open(outfile, "w") if outfile is not None else sys.stdout
		ignore_progeny = False
		for term in obo.iter_terms():
			if term.is_acceptable:
				outline = [str(term)]
				outline += list(term.get_progeny_genes() if not ignore_progeny else term.genes)
				fh.write("\t".join(outline) + "\n")
				myindex = 1
				if desc:
					myt = outline[0]
					myt = re.sub("\s+\[BP\]", "[BP]", myt)
					myt = re.sub("\s+\[MF\]", "[MF]", myt)
					myt = re.sub("\s+\[CC\]", "[CC]", myt)
					if re.search("\[BP\]", myt):
						myt = myt + "__BP" 
					if re.search("\[MF\]", myt):
						myt = myt + "__MF" 
					if re.search("\[CC\]", myt):
						myt = myt + "__CC"
				else:
					myt = outline[0].split(":")
					myt = ":".join(myt[0:2])
				while myindex < len(outline):
					myg = outline[myindex]
					if not myg in terms:
						terms[myg] = {}
					terms[myg][myt] = ""
					myindex = myindex + 1
		fh.close()


def main():
	args_value = parse_arguments()
	basename = "GO"
	if args_value.namespace:
		if len(args_value.namespace) == 1:
			basename = args_value.namespace[0]

	config.logger.info ("Start build_global_terms process ......")
	
	## Get info ##
	config.logger.info("Start collect annotation info ......")
	abunds, header = utilities.read_data_from_file (args_value.input)
	funcs = collect_ann_func (args_value.annotation)
	out_path = os.path.dirname (args_value.output)
	map_path = os.path.join(out_path, "tmp")
	if not os.path.isdir(map_path):
		os.system("mkdir -p " + map_path)
	map_list = collect_taxa_ann (funcs, header, abunds, args_value.abund, args_value.prev, args_value.coverage, args_value.number, args_value.taxon, map_path, basename)
	config.logger.info("Finished collect annotation info ......")

	# Collect bug-specific terms
	# load obo / report rel type
	myterms = {}
	for map_file in map_list:
		config.logger.info("Process " + map_file)
		obo = geneontology.Ontology (args_value.obo)
		myout = re.sub(".tsv", ".term_list.tsv", map_file)
		select_bug_term (obo, map_file, args_value.informative, basename, myterms, myout, None)

	# Combined bug-specific terms
	"""
	<<TEST
	out_path = os.path.dirname (args_value.output)
	map_path = os.path.join(out_path, "tmp")
	map_list = glob.glob(map_path + "/*.genemap.tsv")
	myterms = {}
	for map_file in map_list:
		myout = re.sub(".tsv", ".term_list.tsv", map_file)
		if not os.path.isfile(myout):
			config.logger.info("Error! File doesn't exist: " + myout)
		open_file = open(myout, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			myid = info[0].split(":")
			if len(myid) > 1:
				myid = ":".join(myid[0:2])
			else:
				myid = myid[0]
			i = 1
			while i < len(info):
				myv = info[i]
				if not myv in myterms:
					myterms[myv] = {}
				myterms[myv][myid] = ""
				i = i + 1
		open_file.close()
	TEST
	"""
	
	all_map_file = os.path.join (map_path, "combined.all_bugs." + basename + ".genemap.tsv")
	if os.path.isfile(all_map_file):
		os.system("rm -f " + all_map_file)
	for map_file in map_list:
		os.system("cat " + map_file + " >> " + all_map_file)
	all_file = os.path.join (map_path, "combined.all_bugs_propagated." + basename + ".term_list.tsv")
	obo = geneontology.Ontology (args_value.obo)
	func_all = {}
	select_bug_term (obo, all_map_file, None, args_value.namespace, func_all, all_file, None)
	raw_file = os.path.join (map_path, "combined.bug_specific_informative." + basename + ".genemap.tsv")
	write_funcs_combined (func_all, myterms, raw_file)

	"""
	<<TEST
	term_list = re.sub(".tsv", ".term_list.tsv", args_value.output)
	out_path = os.path.dirname (args_value.output)
	map_path = os.path.join(out_path, "tmp")
	raw_file = os.path.join (map_path, "combined.bug_specific_informative." + basename + ".genemap.tsv")
	TEST
	"""

	final_terms = {}
	term_list = re.sub(".tsv", ".term_list.tsv", args_value.output)
	obo = geneontology.Ontology (args_value.obo)
	select_bug_term (obo, raw_file, args_value.informative, args_value.namespace, final_terms, term_list, True)
	write_funcs_final (final_terms, args_value.output)

	config.logger.info ("Successfully finish build_global_terms process!")

if __name__ == '__main__':
	main()
