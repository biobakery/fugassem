#!/usr/bin/env python
##########################################################################
# Function: Extract subset of total set 
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 08/23/2021
##########################################################################
import sys
import os
import os.path
import re
import argparse
import logging

try:
	from fugassem import config,utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find FUGAsseM python package." +
	" Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Extract subset of total set
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', 
						help = '[REQUIRED] input total file', 
						required = True)
	parser.add_argument('-l', 
						help = '[REQUIRED] input subset list of rowdict file', 
						required = True)
	parser.add_argument('-t', 
						help = '[OPTIONAL] whether include header or not', 
						choices = ["yes", "no"],
						default = "yes")
	parser.add_argument('--pair', 
						help = 'If speficied, input file including gene pair evidence and should check both of them',
						action = "store_true")
	parser.add_argument('-c', 
						help = '[OPTIONAL] input subset list of coldict file', 
						default = None)    
	parser.add_argument('-o', 
						help = '[REQUIRED] output subset file', 
						required = True)    
	values = parser.parse_args()
	return values
# get_args



#==============================================================
# collect cluster ID
#==============================================================
def collect_subset_info (total_file, sub_row, sub_col, header_flag, pair_flag):
	rows = {}
	if not os.path.isfile(sub_row):
		config.logger.info("Error! Subset list file for rowdict doesn't exist.")
	else:
		open_file = open(sub_row, "r")
		for line in open_file:
			line = line.strip()
			if not len(line):
				continue
			info = line.split()
			rows[info[-1]] = ""
		# foreach line
		open_file.close()
	cols = {}
	if sub_col:
		if not os.path.isfile(sub_col):
			config.logger.info("Error! Subset list file for coldict doesn't exist.")
		else:
			open_file = open(sub_col, "r")
			for line in open_file:
				line = line.strip()
				if not len(line):
					continue
				info = line.split()
				cols[info[-1]] = ""
		open_file.close()

	cluster = {}
	titles = {}
	flag_t = 0
	for line in utilities.gzip_bzip2_biom_open_readlines(total_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag_t == 0:
			flag_t = 1
			if header_flag == "yes":
				title = line
				for i in info:
					titles[info.index(i)] = i
				continue
			else:
				title = "NA"
		myname = info[0].split(utilities.c_strat_delim)
		myname = myname[0]
		if myname in rows:
			if pair_flag:
				if len(info) > 2:
					myid2 = info[1].split(utilities.c_strat_delim)[0]
					if not myid2 in rows:
						continue
			if len(titles.keys()) > 0:
				if sub_col:
					title = titles[0]
					myline = info[0]
					myindex = 1
					while myindex < len(info):
						myc = titles[myindex]
						if myc.split("__")[0] in cols:
							title = title + "\t" + myc
							myline = myline + "\t" + info[myindex]
						myindex = myindex + 1
					line = myline
			if not myname in cluster:
				cluster[myname] = []
			cluster[myname].append(line)
	# foreach line
	open_file.close()

	return rows, cluster, title
# collect_subset_info


#==============================================================
# output subset info
#==============================================================
def output_subset (ids, cluster, title, outfile):
	open_out = open(outfile, "w")
	if title != "NA":
		open_out.write(title + "\n")
	for myid in sorted(ids.keys()):
		if myid in cluster:
			for item in cluster[myid]:
				open_out.write(item + "\n")
	open_out.close()
# output_subset 


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	config.logger.info("Start extract_feature_subset.py")

	config.logger.info("Get cluster info ......starting")
	ids, cluster, title = collect_subset_info (values.i, values.l, values.c, values.t, values.pair)
	config.logger.info("Get cluster info ......done")
	
	config.logger.info("Output info ......starting")
	output_subset (ids, cluster, title, values.o)
	config.logger.info("Output info ......done")

	config.logger.info("### Finish extract_feature_subset.py")

# end: main

if __name__ == '__main__':
	main()
