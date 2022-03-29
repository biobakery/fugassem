#!/usr/bin/env python

"""
FUGAsseM: stratified_abundance module
Stratify protein families into taxon

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


try:
	from fugassem import config
	from fugassem import utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find FUGAsseM python package." +
	         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Stratify protein families into taxon
"""

def get_args ():	
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', 
						help = 'input abundance file', 
						required = True)
	parser.add_argument('-m', 
						help = 'input taxonomy file', 
						required = True)
	parser.add_argument('-r', 
						help = 'specify taxonomic level', 
						choices = ["Species", "MSP"],
						default = "Species")
	parser.add_argument('-t', 
						help = 'threshold for taxon presence (e.g. >X percent of genes non-zero values)', 
						default = None)
	parser.add_argument('-o', 
						help = 'output stratified abundance file', 
						required = True)
	values = parser.parse_args()
	
	return values
# get_args


#==============================================================
# collect cluster info
#==============================================================
def collect_protein_cluster_info (clust_file):	
	cluster = {}
	open_file = open(clust_file, "r")
	myclust = ""
	myrep = ""
	for line in open_file.readlines():
		line = line.strip()
		if not len(line):
			continue
		if re.search("^>", line):
			mym = re.search("cluster=([\d]+)", line)
			myclust = "Cluster_" + mym.group(1)
			#myclust = mym.group(1)
			mym = re.search("^>([^\;]+)", line)
			myrep = mym.group(1)
			continue
		mym = re.search("^([\S]+)", line)
		myid = mym.group(1)
		cluster[myid] = myclust
	# foreach line
	open_file.close()
	return cluster
# function collect_gene_cluster_info


#==============================================================
# collect taxa info
#==============================================================
def collect_taxa_info (taxa_file, taxa_level):
	taxa = {}
	taxa_num = {}
	titles = {}
	taxa_levels = ["Terminal", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"]
	flag_t = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (taxa_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			for i in info:
				titles[i] = info.index(i)
			continue
		myid = info[0]
		mytaxa = "NA"
		mylevel = "NA"
		mylineage = info[titles["taxa_lineage"]]
		mylineage = re.sub("\.", "", mylineage)
		mylineage = re.sub("\/", "_", mylineage)
		mylineage = re.sub("\:", "_", mylineage)
		mylineage = re.sub("\[", "", mylineage)
		mylineage = re.sub("\]", "", mylineage)
		mylineage = re.sub("\|", utilities.c_taxon_delim, mylineage)
		if taxa_level == "MSP":
			if "msp_name" in titles:
				if info[titles["msp_name"]] == utilities.c_msp_unknown or info[titles["msp_name"]] == "NA" or info[titles["msp_name"]] == "Unclassified":
					continue
				mytaxa = mylineage + ".msp__" + info[titles["msp_name"]]
		else:
			if taxa_level not in taxa_levels:
				sys.exit("Please provide valid taxonomy level!")
			if taxa_level == "Terminal":
				if re.search("t__", mylineage):
					mytaxa = re.sub("\.msp__[\s\S]+", "", mylineage)
			if taxa_level == "Species":
				if re.search("s__", mylineage):
					mytaxa = re.sub("\.t__[\s\S]+", "", mylineage)
			if taxa_level == "Genus":
				if re.search("g__", mylineage):
					mytaxa = re.sub("\.s__[\s\S]+", "", mylineage)
			if taxa_level == "Family":
				if re.search("f__", mylineage):
					mytaxa = re.sub("\.g__[\s\S]+", "", mylineage)
			if taxa_level == "Order":
				if re.search("o__", mylineage):
					mytaxa = re.sub("\.f__[\s\S]+", "", mylineage)
			if taxa_level == "Class":
				if re.search("c__", mylineage):
					mytaxa = re.sub("\.o__[\s\S]+", "", mylineage)
			if taxa_level == "Phylum":
				if re.search("p__", mylineage):
					mytaxa = re.sub("\.c__[\s\S]+", "", mylineage)
			if taxa_level == "Kingdom":
				if re.search("k__", mylineage):
					mytaxa = re.sub("\.p__[\s\S]+", "", mylineage)
		if mytaxa == "NA" or mytaxa == "Unclassified" or mytaxa == utilities.c_msp_unknown:
			continue
		taxa[myid] = mytaxa
		if not mytaxa in taxa_num:
			taxa_num[mytaxa] = {}
		taxa_num[mytaxa][myid] = ""
	# foreach line

	return taxa, taxa_num
# collect_taxa_info


#==============================================================
# sum up to family abundance
#==============================================================
def sum_abundance (taxa, taxa_num, flt_presence, infile, outfile):
	titles = {}
	samples = []
	taxa_presence = {}
	outs = {}
	flag_t = 0
	for line in utilities.gzip_bzip2_biom_open_readlines (infile):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			myindex = 1
			while myindex < len(info):
				item = info[myindex]
				titles[myindex] = item
				samples.append(item)
				myindex = myindex + 1
			continue
		myid = info[0]
		if myid in taxa:
			mytaxa = taxa[myid]
		else:
			continue
		myclust = myid + utilities.c_strat_delim + mytaxa
		if not mytaxa in taxa_presence:
			taxa_presence[mytaxa] = {}
		taxa_presence[mytaxa][myid] = ""
		if not myclust in outs:
			outs[myclust] = {}
		myindex = 1
		while myindex < len(info):
			mycount = info[myindex]
			mys = titles[myindex]
			if not mys in outs[myclust]:
				outs[myclust][mys] = 0
			if re.search("\.", mycount):
				outs[myclust][mys] = outs[myclust][mys] + float(mycount)
			else:
				outs[myclust][mys] = outs[myclust][mys] + int(mycount)
			myindex = myindex + 1
		# foreach sample
	# foreach line

	# threshold for species presence (e.g. >X% of genes non-zero values)
	refined_taxa = {}
	if flt_presence and flt_presence != "None":
		try:
			flt_presence = float(flt_presence)
			for mytaxa in taxa_presence:
				mynum = len(taxa_presence[mytaxa].keys())
				if mytaxa in taxa_num:
					mytotal = len(taxa_num[mytaxa].keys())
					myper = float(mynum) / float(mytotal)
					if myper >= float(flt_presence):
						refined_taxa[mytaxa] = ""
		except:
			config.logger.info ("Error! Please provide valid value of present fraction within taxon")

	# output abundance
	#samples = list(set(samples)) 
	open_out = open(outfile, "w")
	mytitle = "ID\t" + "\t".join(samples)
	open_out.write(mytitle + "\n")
	for myclust in sorted(outs.keys()):
		mystr = myclust
		tmp = myclust.split(utilities.c_strat_delim)
		if len(refined_taxa.keys()) > 0:
			if not tmp[-1] in refined_taxa:
				continue
		for mys in samples:
			mycount = 0
			if mys in outs[myclust]:
				mycount = outs[myclust][mys]
			mystr = mystr + "\t" + str(mycount)
		# foreach samples
		open_out.write(mystr + "\n")
	# foreach protein family
	open_out.close()
# assign_counts


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args()


	config.logger.info("Start stratified_abundance process")
	

	### collect cluster info ###
	config.logger.info("Get cluster info ......starting")
	taxa, taxa_num = collect_taxa_info (values.m, values.r)
	config.logger.info("Get cluster info ......done")

	### sum abundance to protein families ###
	config.logger.info("Assign counts to protein families ......starting")
	sum_abundance (taxa, taxa_num, values.t, values.i, values.o)
	config.logger.info("Assign counts to protein families ......done")

	config.logger.info("Finish stratified_abundance")

# end: main

if __name__ == '__main__':
	main()
