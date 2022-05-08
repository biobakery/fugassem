#!/usr/bin/env python

"""
Process FUGAsseM tasks for each taxon

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

import os
import sys
import re
import argparse
import random
import logging
import math
from itertools import repeat

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of FUGAsseM tasks
try:
	from fugassem.tasks import fugassem_preprocessing
	from fugassem.tasks import fugassem_prediction
	from fugassem import config, utilities
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the FUGAsseM python package." +
		         " Please check your install.")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Process the main tasks of FUGAsseM for each taxon
"""


def get_args():
	'''
	Process the main tasks of FUGAsseM for each taxon.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	workflow = Workflow(version = config.version, description = "Process the main tasks of FUGAsseM for each taxon", remove_options=["output"])
	# add the custom arguments to the workflow
	workflow.add_argument("minimum-prevalence",
	                      desc = "minimum prevalence of each gene/protein family in normalized-abund MTX [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-abundance",
	                      desc = "minimum abundance of each gene/protein family in normalized-abund MTX [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-detected",
	                      desc = "minimum detected value for each covariate taxon [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("filtering-zero",
	                      desc = "method for pre-filtering zeros in normalized-abund MTX [ Default: None ]",
	                      choices = ["lenient", "semi-strict", "strict", "None"],
	                      default = "None")
	workflow.add_argument("covariate-taxon",
	                      desc = "covariate-taxon abundance table used for pre-filtering zeros [ Default: None ]",
	                      default = None)
	workflow.add_argument("correlation-method",
	                      desc = "correlation method used for co-expression analysis [ Default: Pearson_SE ]",
	                      choices = ["Pearson", "Spearman", "Kendall", "Pearson_SE"],
	                      default = "Pearson_SE")
	workflow.add_argument("go-level",
	                      desc = "GO informative level used for trimming terms that are informative at a given level [ Default: 50 ]:\n"
	                             "<number OR fraction of genes>: spcify numeric level for triming\n"
	                             "<none>: skip trimming\n" 
	                             "<all>: keep all terms",
	                      default = "none")
	workflow.add_argument("go-mode",
	                      desc = "type of GO set used for prediction: [bug-specific] bug-specific informative terms, [universal] universal informative terms, [union] merged bug informative terms for overall prediction",
	                      choices = ["bug-specific", "universal", "union"],
	                      default = "bug-specific")
	workflow.add_argument("func-type",
	                      desc = "GO catgeroy used for prediction [ Default: GO ]",
	                      choices = ["GO", "BP", "CC", "MF"],
	                      default = "GO")
	workflow.add_argument("ml-type",
	                      desc = "machine learning method for function prediction [ Default: RF ]:\n"
	                             "[RF]: Random Forest, [NB] Naive Bayes, [DT] Decision Tree",
	                      choices = ["RF", "NB", "DT"],
	                      default = "RF")
	workflow.add_argument("vector-list",
	                      desc = "[string] a comma separated list of vector-based evidence files (e.g. gene-over-function type of file) [ Default: None ]",
	                      default = None)
	workflow.add_argument("matrix-list",
	                      desc = "[string] a comma separated list of matrix-based evidence files (e.g. gene-by-gene network type of file) [ Default: None ]",
	                      default = None)
	workflow.add_argument("pair-flag",
	                      desc = "[Boolean] whether matrix-based evidence is paired or not [ Default: yes ]",
	                      default = "yes")
	workflow.add_argument("taxon",
	                      desc = '[REQUIRED] specify taxon name',
	                      required = True)
	workflow.add_argument("basename",
	                      desc="specify the basename for output files [ Default: fugassem ]",
	                      default = "fugassem")
	workflow.add_argument("gene",
	                      desc = '[REQUIRED] input list file of genes for a given taxon',
	                      required = True)
	workflow.add_argument("function",
	                      desc = '[REQUIRED] input list file of functions for a given taxon',
	                      required = True)
	workflow.add_argument("bypass-preprocessing",
	                      desc = "do not run the module for preprocessing evidences and functions for prediction")
	workflow.add_argument("bypass-prediction",
	                      desc = "do not run the module for predicting functions")
	workflow.add_argument("bypass-mtx",
	                      desc = "do not integrate MTX for finalized predicting functions")
	workflow.add_argument("threads",
	                      desc = "number of threads/cores for each task to use",
	                      default = 1)
	workflow.add_argument('memory',
	                      desc = 'The amount of memory to use for each fugassem job. Provided in MB',
	                      default = '10240')
	workflow.add_argument('time',
	                      desc = 'The amount of time to use for each fugassem job. Provided in minute',
	                      default = '600')
	workflow.add_argument("output",
	                      desc = '[REQUIRED] output directary name',
	                      required = True)

	return workflow

# get_args


# ==============================================================
###########  Main processing ############
# ==============================================================
def process (workflow, vector_list, matrix_list):
	### get arguments ###
	args = workflow.parse_args()

	# input files
	abund_file = args.input
	gene_file = args.gene
	func_file = args.function
	go_obo = config.go_obo

	# output files
	#mybasename = args.basename + "." + args.taxon
	#myoutput_dir = os.path.join(args.output, args.taxon)
	mybasename = args.basename
	myoutput_dir = args.output
	preprocess_dir = os.path.join(myoutput_dir, "preprocessing")
	predict_dir = os.path.join(myoutput_dir, "prediction")

	# run preprocessing module
	final_func_file = os.path.join(preprocess_dir, mybasename + ".final_funcs.tsv")
	final_func_smp_file = os.path.join(preprocess_dir, mybasename + ".final_funcs.simple.tsv")
	final_funclist_file = os.path.join(preprocess_dir, mybasename + ".final_funclist.txt")
	feature_list_file = os.path.join(preprocess_dir, mybasename + ".feature_list.txt")
	if not args.bypass_preprocessing or args.bypass_preprocessing == "False":
		config.logger.info("Start to run preprocessing module for " + args.taxon + "......")
		flag = 0
		try:
			if os.stat(abund_file).st_size == 0:
				config.logger.info("ERROR! Valid abundance file is not available: " + abund_file)
				flag = flag + 1
		except:
			#tmp_file = re.sub(".tsv$", ".raw.tsv", abund_file)
			tmp_file = abund_file + ".raw.tsv"
			if os.stat(tmp_file).st_size == 0:
				config.logger.info("ERROR! Valid abundance file is not available: " + tmp_file)
				flag = flag + 1
		try:
			if os.stat(gene_file).st_size == 0:
				config.logger.info("ERROR! Valid gene-list file is not available: " + gene_file)
				flag = flag + 1
		except:
			#tmp_file = re.sub(".txt$", ".raw.txt", gene_file)
			tmp_file = gene_file + ".raw.txt"
			if os.stat(tmp_file).st_size == 0:
				config.logger.info("ERROR! Valid gene-list file is not available: " + tmp_file)
				flag = flag + 1
		if os.stat(func_file).st_size == 0:
			config.logger.info("ERROR! Valid func-list file is not available: " + func_file)
			flag = flag + 1
		if flag == 0:
			feature_list = {}
			feature_list, final_func_smp_file, final_func_file, final_funclist_file = fugassem_preprocessing.preprocessing_task(
				abund_file, gene_file, func_file,
				args.go_level, args.func_type, go_obo,
				args.minimum_prevalence,
				args.minimum_abundance,
				args.minimum_detected,
				args.filtering_zero,
				args.covariate_taxon,
				args.correlation_method,
				vector_list, matrix_list, args.pair_flag,
				preprocess_dir, mybasename, feature_list, feature_list_file,
				workflow, args.threads, args.time, args.memory)
		else:
			config.logger.info("WARNING! Bypass module: preprocessing module is skipped......" + args.taxon)
	else:
		config.logger.info("WARNING! Bypass module: preprocessing module is skipped......")
		if not os.path.isfile(feature_list_file):
			config.logger.info ("ERROR! List file of evidence-feature files is not available. Please run preprocessing module......")
		else:
			feature_list = utilities.file_to_dict(feature_list_file)

	# run prediction module
	if not args.bypass_prediction or args.bypass_prediction == "False":
		if args.bypass_mtx == "False" or not args.bypass_mtx:
			args.bypass_mtx = False
		else:
			args.bypass_mtx = True
		config.logger.info("Start to run prediction module for " + args.taxon + "......")
		prediction_list = {}
		final_pred_file = os.path.join(predict_dir, "finalized", mybasename + ".finalized_ML.prediction.tsv")
		prediction_list, final_pred_file = fugassem_prediction.prediction_task(final_func_file, final_funclist_file,
		                                                              feature_list,
		                                                              args.ml_type, args.func_type,
		                                                              predict_dir, mybasename, prediction_list,
		                                                              workflow, args.threads, args.time, args.memory, args.bypass_mtx)

	else:
		config.logger.info("WARNING! Bypass module: prediction module is skipped......")

	
	## start the workflow
	workflow.go()
	

def main():
	workflow = get_args()
	myraw = os.getcwd()
	values = workflow.parse_args()
	myoutput_dir = os.path.abspath(values.output)
	if not os.path.isdir(myoutput_dir):
		os.system("mkdir -p " + myoutput_dir)
	vector_list = values.vector_list
	matrix_list = values.matrix_list
	if vector_list and vector_list != "None":
		vectors = vector_list.split(",")
		vector_list = []
		for vector_file in vectors:
			vector_file = os.path.abspath(vector_file)
			vector_list.append(vector_file)
		vector_list = ",".join(vector_list)
	if matrix_list and matrix_list != "None":
		matrixs = matrix_list.split(",")
		matrix_list = []
		for matrix_file in matrixs:
			matrix_file = os.path.abspath(matrix_file)
			matrix_list.append(matrix_file)
		matrix_list = ",".join(matrix_list)
	os.chdir(myoutput_dir)
	process(workflow, vector_list, matrix_list)
	os.chdir(myraw)


if __name__ == "__main__":
	main()
