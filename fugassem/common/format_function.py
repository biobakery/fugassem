#!/usr/bin/env python
##########################################################################
# Function: Format input function list
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 03/23/2022
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
Format input function list
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-i', 
						help = '[REQUIRED] input file', 
						required = True)
	parser.add_argument('-o', 
						help = '[REQUIRED] output file', 
						required = True)    
	values = parser.parse_args()
	return values
# get_args



def format_func_info (infile, outfile):
	if not os.path.isfile(infile):
		config.logger.info("Error! Input function file doesn't exist.")
	else:
		open_out = open(outfile, "w")
		for line in utilities.gzip_bzip2_biom_open_readlines(infile):
			line = line.strip()
			if not len(line):
				continue
			info = line.split("\t")
			info[-1] = info[-1].split()[0]
			if re.search("^GO\:", info[-1]):
				tmp = info[-1].split(":")
				if len(tmp) >= 2:
					info[-1] = ":".join(tmp[0:2])
			mystr = "\t".join(info)
			open_out.write(mystr + "\n")
		open_out.close()


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	config.logger.info("Start format_function process")

	format_func_info (values.i, values.o)	

	config.logger.info("### Finish format_function process")

# end: main

if __name__ == '__main__':
	main()
