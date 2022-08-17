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
from anadama2.tracked import TrackedExecutable, TrackedDirectory

try:
	from fugassem.tasks import fugassem_preprocessing
	from fugassem.tasks import fugassem_process 
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
	                      choices = ["MSP", "Terminal", "Species", "Genus", "Family", "Order", "Class", "Phylum"],
	                      default = "Species")
	workflow.add_argument("minimum-prevalence",
	                      desc = "minimum prevalence of each gene/protein family in normalized-abund MTX [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-abundance",
	                      desc = "minimum abundance of each gene/protein family in normalized-abund MTX [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-detected",
	                      desc = "minimum detected value for each covariate taxon [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-coverage",
	                      desc = "minimum fraction of annotated genes per taxon [ Default: 0 ]",
	                      default = 0)
	workflow.add_argument("minimum-number",
	                      desc = "minimum number of total genes per taxon [ Default: 0 ]",
	                      default = 0)
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
	                      desc = "GO informative level used for trimming terms that are informative at a given level [ Default: none ]:\n"
	                             "<number OR fraction of genes>: specify numeric level for trimming\n"
	                             "<none>: skip trimming\n" 
	                             "<all>: keep all terms",
	                      default = "none")
	workflow.add_argument("go-mode",
	                      desc = "mode of GO set used for prediction: [bug-specific] bug-specific informative terms, [universal] universal informative terms, [union] merged bug informative terms for overall prediction",
	                      choices = ["bug-specific", "universal", "union"],
	                      default = "bug-specific")
	workflow.add_argument("func-type",
	                      desc = "GO category used for prediction [ Default: GO ]",
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
	workflow.add_argument("matrix-pair",
	                      desc = "treat matrix-based evidence as pair-based network", 
	                      action = "store_true")
	workflow.add_argument("basename",
	                      desc="specify the basename for output files [ Default: fugassem ]",
	                      default = "fugassem")
	workflow.add_argument("input-annotation",
	                      desc = "input the function annotation file for all gene/protein families",
	                      required = True)
	workflow.add_argument("bypass-preparing-taxa",
	                      desc = "do not run the module for splitting stratified MTX into individual taxa",
	                      action = "store_true")
	workflow.add_argument("bypass-preprocessing",
	                      desc = "do not run the module for preprocessing evidences and functions for prediction",
	                      action = "store_true")
	workflow.add_argument("bypass-prediction",
	                      desc = "do not run the module for predicting functions",
	                      action = "store_true")
	workflow.add_argument("bypass-mtx",
	                      desc = "do not integrate MTX for finalized predicting functions",
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


def fugassem_main (workflow):
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
			config.logger.info ("Please provide valid number of minimum coverage by --minimum-coverage! Otherwise, will set it as 0.01")
			args.minimum_coverage = 0.01
	else:
		args.minimum_coverage = 0.01
	if args.minimum_number:
		try:
			args.minimum_number = int(args.minimum_number)
		except:
			config.logger.info ("Please provide valid number of minimum number by --minimum-number! Otherwise, will set it as 10")
			args.minimum_number = 10
	else:
		args.minimum_number = 10
	
	if args.matrix_pair:
		pair_flag = "yes"
	else:
		pair_flag = "no"

	if args.go_mode == "bug-specific":
		universal_flag = "no"
		go_level_flag = args.go_level
	else:
		universal_flag = args.go_mode
		go_level_flag = "none"
		if args.go_level == "none":
			config.logger.info ("Warning! No informative-level was specified for building universal or union GO set")

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
		config.logger.info ("Error: go-basic.obo doesn't exist: " + go_obo)
		sys.exit("Please provide go-basic.obo file in the installed folder!")
	if not os.path.isfile(mtx_file):
		config.logger.info ("Error: taxon-stratified MTX file doesn't exist: " + mtx_file)
		sys.exit("Please provide taxon-stratified MTX file!")
	if not os.path.isfile(ann_file):
		config.logger.info ("Error: original function annotation file doesn't exist: " + ann_file) 
		sys.exit("Please provide raw function annotation file!")

	## prepare outputs
	output_dir = os.getcwd()
	if args.covariate_taxon:
		args.covariate_taxon = os.path.abspath(args.covariate_taxon)
	if args.output:
		output_dir = os.path.abspath(args.output)
	output_dir_raw = output_dir
	output_dir = os.path.join(output_dir_raw, "main")
	output_dir_merged = os.path.join(output_dir_raw, "merged")
	if not os.path.isdir(output_dir):
		os.system("mkdir -p " + output_dir)
	if not os.path.isdir(output_dir_merged):
		os.system("mkdir -p " + output_dir_merged)
	taxa_file = os.path.join(output_dir, config.taxa_abund_list)
	taxa = {}

	## split MTX into individual taxon
	if not args.bypass_preparing_taxa:
		config.logger.info("Start to run split-taxa module......")
		if universal_flag != "no":
			raw_func_file = ann_file
			if universal_flag == "universal":
				new_func_file = os.path.join(output_dir, os.path.basename(ann_file) + ".universal.tsv")
				new_func_sim_file = os.path.join(output_dir, os.path.basename(ann_file) + ".universal.simple.tsv")
				ann_file = new_func_file
				final_func_file, final_func_smp_file = fugassem_preprocessing.retrieve_function (raw_func_file, "no",
			                                                                                 args.go_level, args.func_type, go_obo,
	                                                                                         output_dir, new_func_file, new_func_sim_file,
	                                                                                         args.threads, args.time, args.memory)
			if universal_flag == "union":
				new_func_file = os.path.join(output_dir, os.path.basename(ann_file) + ".union.tsv")
				new_func_sim_file = os.path.join(output_dir, os.path.basename(ann_file) + ".union.simple.tsv")
				ann_file = new_func_file
				final_func_file, final_func_smp_file = fugassem_preprocessing.retrieve_union_function (mtx_file, raw_func_file, "no",
			                                                                                args.go_level, args.func_type, go_obo,
			                                                                                args.taxon_level, args.minimum_prevalence, args.minimum_abundance, args.minimum_coverage, args.minimum_number,
			                                                                                output_dir, new_func_file, new_func_sim_file,
	                                                                                        args.threads, args.time, args.memory)
		taxa_file, taxa = fugassem_preprocessing.preprocess_taxa(mtx_file, ann_file, args.taxon_level,
			                                                    args.minimum_prevalence, args.minimum_abundance,
			                                                    args.minimum_coverage, args.minimum_number,
			                                                    output_dir,
			                                                    args.threads, args.time, args.memory)
	else:
		config.logger.info("WARNING! Bypass module: preparing-taxa module is skipped......")
		if not os.path.isfile(taxa_file):
			sys.exit("ERROR! Taxon-specific MTX files are not available. Please run preprocess_taxa module......")
		taxa = utilities.file_to_dict (taxa_file)
	
	if args.bypass_preprocessing and args.bypass_prediction:
		config.logger.info("WARNING! Bypass module: main fugassem modules is skipped......")
	else:
		merged_final_files = []
		## add tasks to the workflow
		for mytaxa in sorted(taxa.keys()):
			config.logger.info("####---------- Process " + mytaxa + " -------------####")
			mybasename = basename + "." + mytaxa
			myoutput_dir = os.path.join(output_dir, mytaxa)
			preprocess_dir = os.path.join(myoutput_dir, "preprocessing")
			predict_dir = os.path.join(myoutput_dir, "prediction")
			abund_file = os.path.join(myoutput_dir, "data", mytaxa + config.c_abund_extension)
			gene_file = os.path.join(myoutput_dir, "data", mytaxa + config.c_gene_extension)
			func_file = os.path.join(myoutput_dir, "data", mytaxa + config.c_ann_extension)
			mylog = os.path.join(myoutput_dir, mytaxa + ".fugassem.log")
			final_func_file = os.path.join(preprocess_dir, mybasename + ".final_funcs.tsv")
			final_func_smp_file = os.path.join(preprocess_dir, mybasename + ".final_funcs.simple.tsv")
			final_funclist_file = os.path.join(preprocess_dir, mybasename + ".final_funclist.txt")
			feature_list_file = os.path.join(preprocess_dir, mybasename + ".feature_list.txt")
			#final_pred_file = os.path.join(predict_dir, "finalized", mybasename + ".finalized_ML.prediction.tsv")
			final_pred_file = os.path.join(myoutput_dir, mybasename + ".finalized_ML.prediction.tsv")
			merged_final_files.append(final_pred_file)

			args_list = [mytaxa,
		               mybasename,
		               args.minimum_prevalence,
		               args.minimum_abundance,
		               args.minimum_detected,
		               args.filtering_zero,
		               args.covariate_taxon,
		               args.correlation_method,
		               go_level_flag,
		               args.func_type,
		               args.ml_type,
		               args.vector_list,
		               args.matrix_list,
		               args.bypass_preprocessing,
		               args.bypass_prediction,
		               args.bypass_mtx,
		               args.threads,
		               args.memory,
		               args.time,
		               myoutput_dir,
		               mylog,
					   pair_flag,
			           args.go_mode]
			target_list = [final_func_file, final_func_smp_file, final_funclist_file, feature_list_file, final_pred_file]
			workflow.add_task_gridable(
				"fugassem_process --input [depends[0]] --gene [depends[1]] --function [depends[2]] "
				"--taxon [args[0]] --basename [args[1]] "
				"--minimum-prevalence [args[2]] --minimum-abundance [args[3]] --minimum-detected [args[4]] --filtering-zero [args[5]] --covariate-taxon [args[6]] "
				"--correlation-method [args[7]] --go-level [args[8]] --go-mode [args[22]] --func-type [args[9]] --ml-type [args[10]] "
				"--vector-list [args[11]] --matrix-list [args[12]] --pair-flag [args[21]] "
				"--bypass-preprocessing [args[13]] --bypass-prediction [args[14]] --bypass-mtx [args[15]] "
				"--threads [args[16]] --memory [args[17]] --time [args[18]] "
				"--output [args[19]] > [args[20]] 2>&1",
				depends = [abund_file, gene_file, func_file, TrackedExecutable("fugassem_process")],
				targets = target_list,
				args = args_list,
				cores = args.threads,
				time = args.time,
				mem = args.memory,
				name = "fugassem_process")
		# foreach taxon

		# combined files
		out_list_file =  os.path.join(output_dir_merged, basename +  ".finalized_ML.prediction.list.txt")
		outfile = os.path.join(output_dir_merged, basename +  ".finalized_ML.prediction.tsv")
		mymatch = os.path.join(output_dir_merged, basename + ".feature_maps.txt")
		mylog = os.path.join(output_dir_merged, "merged_finalized_prediction.log")
		utilities.array_to_file (merged_final_files, out_list_file)
		workflow.add_task(
			"fugassem_merged_prediction --input [depends[0]] --basename [args[0]] --output [targets[0]] > [args[1]] 2>&1",
			depends = utilities.add_to_list (out_list_file, TrackedExecutable("fugassem_merged_prediction")) + merged_final_files,
			targets = [outfile, mymatch],
			args = [basename, mylog],
			cores = args.threads,
			time = args.time,
			mem = args.memory,
			name = "fugassem_merged_prediction")
	# if run process module

	## start the workflow
	workflow.go()


def main():
	fugassem_main (parse_cli_arguments())


if __name__ == "__main__":
	main()
