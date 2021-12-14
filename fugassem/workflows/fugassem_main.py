#!/usr/bin/env python

"""
FUGAsseM: Function predictor of Uncharacterized Gene products by Assessing high-dimensional community data in microbiome
1) split MTX into individual taxa (e.g. per species, per MSP)
2) preprocess original annotations and evidences
3) predict functions using machine learning

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
import os, fnmatch
import logging

# import the workflow class from anadama2
from anadama2 import Workflow

# import the library of FUGAsseM tasks
try:
	from fugassem.tasks import preprocessing
	from fugassem.tasks import prediction
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the FUGAsseM python package." +
		         " Please check your install.")

# import the utilities functions and config settings from FUGAsseM
try:
	from fugassem import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the FUGAsseM python package." +
		         " Please check your install.")


VERSION = config.version

def parse_cli_arguments ():
	'''
	Parses any command-line arguments passed into the workflow.
	the version number will appear when running this script with the "--version" option
	the description will appear when running this script with the "--help" option
	create a workflow instance, providing the version number and description
	'''

	working_dir = os.getcwd()
	workflow = Workflow(version = VERSION, description = "FUGAsseM for function prediction", remove_options=["output"])

	# add the custom arguments to the workflow
	workflow.add_argument("taxon-level",
	                      desc = "taxonomic level used for stratification [ Default: Species ]",
	                      choices = ["MSP", "Species", "Genus", "Family", "Order", "Class", "Phylum"],
	                      default = "Species")
	workflow.add_argument("minimum-prevalence",
	                      desc = "minimum prevalence of each gene/protein family in normalized-abund MTX [ Default: 0.01 ]",
	                      default = 0.01)
	workflow.add_argument("minimum-abundance",
	                      desc = "minimum abundance of each gene/protein family in normalized-abund MTX [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-detected",
	                      desc = "minimum detected value for each covariate taxon [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-coverage",
	                      desc = "minimum fraction of annotated genes per taxon [ Default: 0.1 ]",
	                      default = 0.1)
	workflow.add_argument("filtering-zero",
	                      desc = "method for pre-filtering zeros in normalized-abund MTX [ Default: lenient ]",
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
	                      default = 50)
	workflow.add_argument("func-type",
	                      desc = "GO catgeroy used for prediction [ Default: BP ]",
	                      choices = ["BP", "CC", "MF"],
	                      default = "BP")
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
	workflow.add_argument("basename",
	                      desc="specify the basename for output files [ Default: fugassem ]",
	                      default = "fugassem")
	workflow.add_argument("input-annotation",
	                      desc = "input the function annotation file for all gene/protein families",
	                      required = True)
	workflow.add_argument("bypass-preparing-taxa",
	                      desc = "do not run the module for spliting stratified MTX into individual taxa",
	                      action = "store_true")
	workflow.add_argument("bypass-preprocessing",
	                      desc = "do not run the module for preprocessing evidences and functions for prediction",
	                      action = "store_true")
	workflow.add_argument("bypass-prediction",
	                      desc = "do not run the module for predicting functions",
	                      action = "store_true")
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
	                      desc = "provide an output folder which the workflow database and log is written. "
	                             "By default, that be written to the anadama2 folder of users' working directory",
	                      default = working_dir)

	return workflow


def get_method_config (config_file):
	'''
	:param config_file:
	:return: config info for each analysis
	'''
	config.logger.info ("###### Start get_method_config module ######")

	config_items = config.read_user_edit_config_file(config_file)
	family_conf = {}
	domain_motif_conf = {}
	abundance_conf = {}
	integration_conf = {}
	values = ["yes", "no", "Yes", "No"]

	if "global_homology" in config_items:
		for name in config_items["global_homology"].keys():
			myvalue = config_items["global_homology"][name]
			if not myvalue in values:
				config.logger.info ('ERROR! The config value can not be recognized. Please check your config file!')
				continue
			family_conf[name] = myvalue
		# for each method

	if "domain_motif" in config_items:
		for name in config_items["domain_motif"].keys():
			myvalue = config_items["domain_motif"][name]
			if not myvalue in values:
				config.logger.info ('ERROR! The config value can not be recognized. Please check your config file!')
				continue
			domain_motif_conf[name] = myvalue
		# for each method

	if "abundance" in config_items:
		for name in config_items["abundance"].keys():
			myvalue = config_items["abundance"][name]
			abundance_conf[name] = myvalue
		# for each method

	if "integration" in config_items:
		for name in config_items["integration"].keys():
			myvalue = config_items["integration"][name]
			if not myvalue in values:
				config.logger.info ('ERROR! The config value can not be recognized. Please check your config file!')
				continue
			integration_conf[name] = myvalue
		# for each method

	config.logger.info("###### Finish get_method_config module ######")

	return family_conf, domain_motif_conf, abundance_conf, integration_conf



def main(workflow):
	'''
	Function prediction for novel gene products from microbial communities
	'''

	## get arguments
	args = workflow.parse_args()
	tmps_dir = os.path.join (os.getcwd(), "temp")
	tmps_dir = os.path.abspath(tmps_dir)

	## prep settings
	if args.threads:
		try:
			args.threads = int(args.threads)
		except:
			config.logger.info ("Please provide valid number of required threads! Otherwise, will set it as 1")
			args.threads = 1
	else:
		args.threads = 1
	if args.memory:
		try:
			args.memory = int(args.memory)
		except:
			config.logger.info("Please provide valid number of required memory! Otherwise, will set it as 10240")
			args.memory = 10240
	else:
		args.memory = 10240
	if args.time:
		try:
			args.time = int(args.time)
		except:
			config.logger.info("Please provide valid number of required time! Otherwise, will set it as 600")
			args.time = 600
	else:
		args.time = 600
	if args.minimum_prevalence:
		try:
			args.minimum_prevalence = float(args.minimum_prevalence)
		except:
			config.logger.info ("Please provide valid number of minimum prevalence by --minimum-prevalence! Otherwise, will set it as 0.01")
			args.minimum_prevalence = 0.01
	else:
		args.minimum_prevalence = 0.01
	if args.minimum_abundance:
		try:
			args.minimum_abundance = float(args.minimum_abundance)
		except:
			config.logger.info ("Please provide valid number of minimum abundance by --minimum-abundance! Otherwise, will set it as 0")
			args.minimum_abundance = 0
	else:
		args.minimum_abundance = 0
	if args.minimum_detected:
		try:
			args.minimum_detected = float(args.minimum_detected)
		except:
			config.logger.info ("Please provide valid number of minimum abundance by --minimum-detected! Otherwise, will set it as 0")
			args.minimum_detected = 0
	else:
		args.minimum_detected = 0
	if args.minimum_coverage:
		try:
			args.minimum_coverage = float(args.minimum_coverage)
		except:
			config.logger.info ("Please provide valid number of minimum abundance by --minimum-coverage! Otherwise, will set it as 0")
			args.minimum_coverage = 0
	else:
		args.minimum_coverage = 0

	## get all input files
	go_obo = config.go_obo
	basename = config.basename
	mtx_file = os.path.abspath(args.input)
	ann_file = None
	if args.input_annotation:
		ann_file = os.path.abspath(args.input_annotation)
	if args.basename:
		basename = args.basename
	if not os.path.isfile(go_obo):
		sys.exit("Please provide go-basic.obo file in the installed folder!")
	if not os.path.isfile(mtx_file):
		sys.exit("Please provide taxon-stratified MTX file!")
	if not os.path.isfile(ann_file):
		sys.exit("Please provide raw function annotation file!")

	## prepare outputs
	output_dir = os.getcwd()
	if args.output:
		output_dir = os.path.abspath(args.output)
	if not os.path.isdir(output_dir):
		os.system("mkdir -p " + output_dir)
	taxa_file = os.path.join(output_dir, config.taxa_abund_list)
	taxa = {}

	## split MTX into individual taxon
	if not args.bypass_preparing_taxa:
		config.logger.info("Start to run split-taxa module......")
		taxa_file, taxa = preprocessing.preprocess_taxa (mtx_file, ann_file, args.taxon_level,
	                                args.minimum_prevalence, args.minimum_abundance, args.minimum_coverage, output_dir,
	                                args.threads, args.time, args.memory)
	else:
		config.logger.info("WARNING! Bypass module: preparing-taxa module is skipped......")
		if not os.path.isfile(taxa_file):
			sys.exit("ERROR! Taxon-specific MTX files are not available. Please run preprocess_taxa module......")
		taxa = utilities.file_to_dict (taxa_file)

	## add tasks to the workflow
	for mytaxa in sorted(taxa.keys()):
		config.logger.info("####---------- Process " + mytaxa + " -------------####")
		mybasename = basename + "." + mytaxa
		myoutput_dir = os.path.join(output_dir, mytaxa)
		preprocess_dir = os.path.join(myoutput_dir, "preprocessing")
		predict_dir = os.path.join(myoutput_dir, "prediction")
		abund_file = os.path.join(myoutput_dir, mytaxa + config.c_abund_extension)
		gene_file = os.path.join(myoutput_dir, mytaxa + config.c_gene_extension)
		func_file = os.path.join(myoutput_dir, mytaxa + config.c_ann_extension)

		# run preprocessing module
		final_func_file = os.path.join(preprocess_dir, mybasename + ".final_funcs.tsv")
		final_func_smp_file = os.path.join(preprocess_dir, mybasename + ".final_simplified_funcs.tsv")
		final_funclist_file = os.path.join(preprocess_dir, mybasename + ".final_funclist.txt")
		feature_list_file = os.path.join(preprocess_dir, basename + ".feature_list.txt")
		if not args.bypass_preprocessing:
			config.logger.info ("Start to run preprocessing module for " + mytaxa + "......")
			feature_list = {}
			feature_list, final_func_smp_file, final_func_file, final_funclist_file = preprocessing.preprocessing_task (abund_file, gene_file, func_file,
		                                                                                              args.go_level, args.func_type, go_obo,
		                                                                                              args.minimum_prevalence,
		                                                                                              args.minimum_abundance,
		                                                                                              args.minimum_detected,
		                                                                                              args.filtering_zero,
		                                                                                              args.covariate_taxon,
		                                                                                              args.correlation_method,
		                                                                                              args.vector_list, args.matrix_list,
		                                                                                              preprocess_dir, mybasename, feature_list, feature_list_file,
		                                                                                              workflow, args.threads, args.time, args.memory)
		else:
			config.logger.info("WARNING! Bypass module: preprocessing module is skipped......")
			if not os.path.isfile(feature_list_file):
				sys.exit("ERROR! List file of evidence-feature files is not available. Please run preprocessing module......")
			feature_list = utilities.file_to_dict(feature_list_file)

		# run prediction module
		if not args.bypass_prediction:
			config.logger.info("Start to run prediction module for " + mytaxa + "......")
			prediction_list = {}
			final_pred_file = os.path.join(predict_dir, "finalized", mybasename + ".finalized_ML.prediction.tsv")
			prediction_list, final_pred_file = prediction.prediction_task (final_func_file, final_funclist_file, feature_list,
		                                                    args.ml_type, args.func_type,
		                                                    predict_dir, mybasename, prediction_list,
		                                                    workflow, args.threads, args.time, args.memory)

		else:
			config.logger.info("WARNING! Bypass module: prediction module is skipped......")

	## start the workflow
	workflow.go()


if __name__ == "__main__":
	main(parse_cli_arguments())



