"""
FUGAsseM Workflow: prediction module
A collection of tasks for function prediction

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
from os import path
import subprocess
import itertools
import re

from anadama2.tracked import TrackedExecutable, TrackedDirectory

# import the utilities functions and config settings from MetaWIBELE
try:
	from fugassem import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the FUGAsseM python package." +
		         " Please check your install.")


def precess_ML (feature_file, func_file, funclist_file, ml_type, feature_suffix, vector,
                output_folder, pred_prefix, prediction_file,
                workflow, threads, time_equation, mem_equation,
				func_type):
	"""
	Train and predict for each function using machine learner

	Args:
		feature_file: feature file for ML.
		func_file: raw annotated function file.
		funclist_file: function list file for classifying.
		ml_type: specify type of ML method, [RF]: Random Forest, [NB] Naive Bayes, [DT] Decision Tree; [ Default: RF ]
		feature_suffix: suffix name of the vector feature. [ Deafult: None]
		vector: [store_true] add vector from metadata to features
		output_folder (string): The path of the output folder.
		pred_prefix: prefix of output.
		prediction_file: prediction file.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): required number of hours defined in the workflow.
		mem_equation (int): required number of GB defined in the workflow.
		func_type: function type

	Requires:
		feature file for ML

	Returns:
		string: the file of prediction file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import precess_ML

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add prepare_function tasks
		prediction_file = precess_ML (
						feature_file,
						func_file,
						funclist_file,
						ml_type,
						feature_suffix,
						vector,
						output_dir,
						pred_prefix,
						prediction_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation,
						func_type)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start precess_ML module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	prediction_folder = os.path.join (main_folder, pred_prefix + "_results")
	prediction_list = os.path.join (prediction_folder, "prediction_list.txt")
	prediction_file = os.path.join(main_folder, os.path.basename(prediction_file))
	xval_file = re.sub(".tsv$", ".xval.tsv", prediction_file)
	xval_yes_file = re.sub(".tsv$", ".xval.yes.tsv", prediction_file)
	xval_no_file = re.sub(".tsv$", ".xval.no.tsv", prediction_file)

	prediction_log = os.path.join (output_folder, "machine_learning.log")
	if vector:
		workflow.add_task (
			"fugassem_machine_learning -f [depends[0]] -l [depends[1]] -t [args[0]] -s [args[1]] --vector -c [args[2]] -o [args[3]] > [args[4]] 2>&1",
			depends = [feature_file, funclist_file, func_file, TrackedExecutable("fugassem_machine_learning")],
			targets = [prediction_list],
			args = [ml_type, feature_suffix, threads, prediction_folder, prediction_log],
			cores = threads,
			time = time_equation,
			mem = mem_equation,
			name = "fugassem_machine_learning")
	else:
		workflow.add_task (
			"fugassem_machine_learning -f [depends[0]] -l [depends[1]] -t [args[0]] -s [args[1]] -c [args[2]] -o [args[3]] > [args[4]] 2>&1",
			depends = [feature_file, funclist_file, func_file, TrackedExecutable("fugassem_machine_learning")],
			targets = [prediction_list],
			args = [ml_type, feature_suffix, threads, prediction_folder, prediction_log],
			cores = threads,
			time = time_equation,
			mem = mem_equation,
			name = "fugassem_machine_learning")

	collect_log = os.path.join(output_folder, "collect_ml_results.log")
	myextension = ".xval.roc.tsv"
	mytype = "xval"
	workflow.add_task(
		"fugassem_collect_ml_results -p [args[0]] -e [args[1]] -t [args[2]] -o [args[3]] > [args[4]] 2>&1",
		depends = [prediction_list, funclist_file, TrackedExecutable("fugassem_collect_ml_results")],
		targets = [xval_yes_file, xval_no_file],
		args = [prediction_folder, myextension, mytype, xval_file, collect_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_collect_ml_results")

	predict_log = os.path.join(output_folder, "predict_function.log")
	workflow.add_task(
		"fugassem_predict_function -a [depends[0]] -l [depends[1]] -b [depends[2]] -m ML -f no -s 0 -e 1 -c [args[0]] -o [targets[0]] > [args[1]] 2>&1",
		depends = [func_file, funclist_file, xval_yes_file, TrackedExecutable("fugassem_predict_function")],
		targets = [prediction_file],
		args = [func_type, predict_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_predict_function")

	return prediction_file


def integrate_base_learning (feature_file, learnlist_file, learning_list,
                output_folder, combined_feature_file,
                workflow, threads, time_equation, mem_equation):
	"""
	Integrate prediction evidences from individual classifier

	Args:
		feature_file: feature file for ML.
		learnlist_file: prediction list file for individual ML learner.
		learning_list: a list of prediction file names.
		output_folder (string): The path of the output folder.
		combined_feature_file: integrated prediction evidence features for next layer of ML.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		feature file for ML

	Returns:
		string: the file of prediction file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import integrate_base_learning

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add prepare_function tasks
		combined_feature_file = integrate_base_learning (
						feature_file,
						learnlist_file,
						learning_list,
						output_dir,
						combined_feature_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start integrate_base_learning module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	combined_feature_file = os.path.join(main_folder, os.path.basename(combined_feature_file))
	integration_log = os.path.join (output_folder, "prepare_nextlayer_ml.log")

	myinputs = []
	myinputs.append(feature_file)
	myinputs.append(learnlist_file)
	myinputs.extend(learning_list)
	workflow.add_task(
		"fugassem_prepare_nextlayer_ml -r [depends[0]] -l [depends[1]] -o [targets[0]] > [args[0]] 2>&1",
		depends = utilities.add_to_list(myinputs, TrackedExecutable("fugassem_prepare_nextlayer_ml")),
		targets = [combined_feature_file],
		args = [integration_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_prepare_nextlayer_ml")

	return combined_feature_file


def prediction_task (func_file, funclist_file, feature_list, ml_type, func_type,
                     output_folder, basename, prediction_list,
                     workflow, threads, time_equation, mem_equation, bypass_mtx):
	"""
	Prediction function using machine learning approaches

	Args:
		func_file: original annotated function file.
		funclist_file: function list file for prediction.
				feature_list: a directory of feature files.
		ml_type: specify type of ML method, [RF]: Random Forest, [NB] Naive Bayes, [DT] Decision Tree; [ Default: RF ]
		func_type: function category, e.g. GO | BP | CC | MF
		output_folder (string): The path of the output folder.
		basename: basename for output files.
		prediction_list: a directory of individual prediction file.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.
		prediction_task: whether to skip MTX for final prediction

	Requires:
		function file for one taxon
		function list file for one taxon

	Returns:
		string: the file of finalized prediction file
		dict: a dictionary of predicted files

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import prediction_task

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		prediction_list, final_prediction_file = prediction_task (
						func_file,
						funclist_file,
						feature_list,
						ml_type,
						func_type,
						output_folder,
						basename,
						prediction_list,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation,
						bypass_mtx,
						func_type)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start prediction_task module #####")

	## prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)

	## process individual ML
	pred_folder_list = {}
	for mytype in feature_list.keys():
		feature_file = feature_list[mytype]
		pred_folder = os.path.join(main_folder, mytype)
		pred_prefix = mytype + "_ML"
		final_pred_file = os.path.join(pred_folder, basename + "." + pred_prefix + ".prediction.tsv")
		if re.search("^vector", mytype):
			precess_ML(feature_file, func_file, funclist_file, ml_type, "vector", "yes",
		           pred_folder, pred_prefix, final_pred_file,
		           workflow, threads, time_equation, mem_equation, func_type)
		else:
			precess_ML(feature_file, func_file, funclist_file, ml_type, None, None,
			           pred_folder, pred_prefix, final_pred_file,
			           workflow, threads, time_equation, mem_equation, func_type)
		prediction_list[mytype] = final_pred_file
		pred_folder_list[os.path.join(pred_folder, pred_prefix + "_results")] = ""

	## process the 2nd layer of ML
	pred_list_file = os.path.join(main_folder, basename + ".predfile_list.tsv")
	if bypass_mtx:
		myindex = os.path.join(main_folder, "coexp", "coexp_ML_results")
		try:
			del pred_folder_list[myindex]
		except:
			pass
	utilities.dict_to_file(pred_folder_list, pred_list_file)
	preds = []
	for i in sorted(prediction_list.keys()):
		preds.append(prediction_list[i])
	feature_file = feature_list["coexp"]
	pred_folder = os.path.join(main_folder, "finalized")
	combined_feature_file = os.path.join(pred_folder, basename + ".integrated.features.trans.tsv")
	final_pred_file = os.path.join(pred_folder, basename + ".finalized_ML.prediction.tsv")
	integrate_base_learning (feature_file, pred_list_file, preds,
	                         pred_folder, combined_feature_file,
	                         workflow, threads, time_equation, mem_equation)
	precess_ML(combined_feature_file, func_file, funclist_file, ml_type, "all", "yes",
		        pred_folder, "integrated_ML", final_pred_file,
		        workflow, threads, time_equation, mem_equation, func_type)

	return prediction_list, final_pred_file
