"""
FUGAsseM: preprocessing module
A collection of tasks for preprocessing inputs

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
import subprocess
import itertools
import re

from anadama2.tracked import TrackedExecutable, TrackedDirectory

# import the utilities functions and config settings from FUGAsseM
try:
	from fugassem import utilities, config
except ImportError:
	sys.exit("CRITICAL ERROR: Unable to find the FUGAsseM python package." +
		         " Please check your install.")


def preprocess_taxa (abunds, anns, taxa_level, min_prev, min_abund, min_cov, min_num, output_folder, threads, time_equation, mem_equation):
	"""
	This set of tasks will split MTX into stratified taxa

	Args:
		abunds: MTX file.
		anns: original function annotation file.
		taxa_level: taxonomic level for stratification.
		min_prev: the minimum prevalence per gene [ Default: None ].
		min_abund: the minimum detected abundance for each gene [ Default: 0 ].
		min_cov: the minimum fraction of annotated genes per taxon [ Default: 0.1 ].
		min_num: the minimum number of total genes per taxon [ Default: 500 ].
		output_folder (string): The path of the output folder.
		threads (int): The number of threads/cores for clustering to use.

	Requires:
		MTX file
		annotation file

	Returns:
		dict: stratified taxa list.
		string: the file of stratified-taxon abunds files.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import preprocess_taxa

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_taxa tasks
		taxa_file, taxa = preprocess_taxa (abunds,
						anns,
						taxa_level,
						min_prev,
						min_abund,
						min_cov,
						min_num,
						output_dir,
                        args.threads,
                        time_equation,
                        mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start preprocess_taxa module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)

	taxa_file = os.path.join(main_folder, config.taxa_abund_list)
	output_log = taxa_file + ".log"

	# run tasks to preprocess taxa
	utilities.run_task (
		"fugassem_split_taxa -i [depends[0]] -a [depends[1]] -t [args[0]] -p [args[1]] -c [args[2]] -n [args[3]] -d [args[4]] -o [args[5]] > [args[6]] 2>&1",
		depends = [abunds, anns, TrackedExecutable("fugassem_split_taxa")],
		targets = [taxa_file],
		args = [taxa_level, min_prev, min_cov, min_num, min_abund, main_folder, output_log],
		cores = threads,
		name = "fugassem_split_taxa")

	# collect taxa list
	taxa = utilities.file_to_dict (taxa_file)

	return taxa_file, taxa


def refine_abundance (abund_file, gene_file, min_prev, min_abund, min_detected, zero_flt, taxa_abund_file,
                      output_folder, final_abund_file, final_gene_file, final_family_file,
                      workflow, threads, time_equation, mem_equation):
	"""
	Refine abundance of each taxon

	Args:
		abund_file: MTX abundance file for one taxon.
		gene_file: gene list file
		min_prev: the minimum prevalence per gene [ Default: 0.01 ].
		min_abund: the minimum detected abundance for each gene [ Default: 0 ].
		min_detected: the minimum detected value for each covariate feature [ Default: 0 ]
		zero_flt: pre-filtering approach [ Default: None ], choices: ["lenient", "semi-strict", "strict", "None"]
		taxa_abund_file: species abundance table used for technical-zero filtering [ Default: None ]
		output_folder (string): The path of the output folder.
		final_abund_file: finalized MTX abundance file for one taxon.
		final_gene_file: finalized gene list file.
		final_family_file: finaized family list file.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw MTX abundance file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined abunds file.
		string: gene list file

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import refine_abundance

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add refine_abundance tasks
		abund_file, gene_file = refine_abundance (
						abund_file,
						gene_file,
						min_prev = 0.01,
						min_abund = 0,
						min_detected = 0,
						zero_flt = None,
						taxa_abund_file = None,
						output_dir,
						final_abund_file,
						final_gene_file,
						final_family_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start refine_abundance module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	#gene_file = re.sub(c_abund_extension + "$", c_gene_extension, abund_file)
	gene_file_raw = gene_file
	abund_file_raw = abund_file
	raw_gene = re.sub(".txt$", ".raw.txt", gene_file)
	raw_file = re.sub(".tsv$", ".raw.tsv", abund_file)
	gene_file = os.path.join(main_folder, os.path.basename(gene_file))
	abund_file = os.path.join(main_folder, os.path.basename(abund_file))
	refined_file = os.path.join(main_folder, os.path.basename(re.sub(".tsv$", ".refined.tsv", abund_file)))
	smoothed_file =  os.path.join(main_folder, os.path.basename(re.sub(".tsv$", ".smoothed.tsv", abund_file)))
	transformed_file = os.path.join(main_folder, os.path.basename(re.sub(".tsv$", ".transformed.tsv", abund_file)))
	final_gene_file = os.path.join(main_folder, os.path.basename(final_gene_file))
	final_abund_file = os.path.join(main_folder, os.path.basename(final_abund_file))
	smoothing_log = final_abund_file + ".smoothing.log"
	refining_log = final_abund_file + ".refining.log"

	# add tasks to workflow
	"""
	workflow.add_task ("mv -f [depends[0]] [targets[0]]",
			depends = [gene_file_raw],
	        targets = [raw_gene],
	        cores = 1,
	        name = "mv__backup_genelist")

	workflow.add_task ("mv -f [depends[0]] [targets[0]]",
			depends = [abund_file_raw],
	        targets = [raw_file],
	        cores = 1,
	        name = "mv__backup_abundance")
	"""

	workflow.add_task (
		"fugassem_abundance_smoothing -i [depends[0]] -t log -o [targets[0]] > [args[0]] 2>&1",
		depends = [abund_file_raw, TrackedExecutable("fugassem_abundance_smoothing")],
		targets = [abund_file, refined_file, transformed_file],
		args = [smoothing_log],
		cores = threads,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_abundance_smoothing")

	workflow.add_task ("mv -f [depends[0]] [targets[0]]",
			depends = [abund_file],
	        targets = [smoothed_file],
	        cores = 1,
	        name = "mv__backup_smoothing")

	workflow.add_task ("rm -f [depends[0]]",
			depends = [refined_file],
	        cores = 1,
	        name = "rm__smoothing_intermediate")

	workflow.add_task (
		"fugassem_refine_abunds -i [depends[0]] -r [depends[1]] -c [args[0]] -m [args[1]] "
		"-p [args[2]] -a [args[3]] -d [args[4]] -o [targets[0]] > [args[5]] 2>&1",
		depends = [transformed_file, abund_file_raw, TrackedExecutable("fugassem_refine_abunds")],
		targets = [final_abund_file],
		args = [taxa_abund_file, zero_flt, min_prev, min_abund, min_detected, refining_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_refine_abunds")

	workflow.add_task ('sed -e 1d [depends[0]] | awk -F \"\\t\" \'{print $1}\' > [targets[0]]',
			depends = [final_abund_file],
			targets = [final_gene_file],
	        cores = 1,
	        name = "collect__genelist")

	workflow.add_task ('less [depends[0]] | awk -F \"\\t\" \'{print $1}\' > [targets[0]]',
			depends = [final_abund_file],
			targets = [final_family_file],
	        cores = 1,
	        name = "collect__familylist")

	return final_abund_file, final_gene_file, final_family_file


def calculate_covariation (abund_file, corr_method, output_folder, final_corr_file, final_abund_file, evidence_list, evidence_type,
                           workflow, threads, time_equation, mem_equation):
	"""
	Calculate co-variation of abundance/expression per taxon

	Args:
		abund_file: refined abundance for one taxon.
		corr_method: method for corrleation, e.g. Pearson_SE | Pearson | Spearman | Kendall
		output_folder (string): The path of the output folder.
		final_corr_file: finalized correlation without NaN file for one taxon.
		final_abund_file: finalized abundance without NaN file for one taxon.
		evidence_list: a dictionary recording files of evidences
		evidence_type: type of evidence
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		refined abundance file for one taxon

	Returns:
		string: the file of correlation file.
		string: the file of imputed-nan abundance file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import calculate_covariation

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		final_corr_file, final_abund_file = calculate_covariation (
						abund_file,
						corr_method = "Pearson_SE",
						output_dir,
						final_corr_file,
						final_abund_file,
						evidence_list,
						evidence_type,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start calculate_covariation module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	final_corr_file = os.path.join(main_folder, os.path.basename(final_corr_file))
	tmp_corr_file = re.sub(".tsv$", ".raw_corr.tsv", final_corr_file)
	final_abund_file = os.path.join(main_folder, os.path.basename(final_abund_file))
	corr_log = tmp_corr_file + ".corr.log"
	corr_impute_log = final_corr_file + ".imputed_nan.log"
	abund_impute_log = final_abund_file + ".imputed_nan.log"

	workflow.add_task(
		"fugassem_calculate_correlation -i [depends[0]] -m [args[0]] -c [args[1]] -o [targets[0]] > [args[2]] 2>&1",
		depends = [abund_file, TrackedExecutable("fugassem_calculate_correlation")],
		targets = [tmp_corr_file],
		args = [corr_method, threads, corr_log],
		cores = threads,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_calculate_correlation")

	workflow.add_task(
		"fugassem_impute_nan -i [depends[0]] -m zero -o [targets[0]] > [args[0]] 2>&1",
		depends = [abund_file, TrackedExecutable("fugassem_impute_nan")],
		targets = [final_abund_file],
		args = [abund_impute_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_impute_nan")

	workflow.add_task(
		"fugassem_impute_nan -i [depends[0]] -m mean -o [targets[0]] > [args[0]] 2>&1",
		depends = [tmp_corr_file, TrackedExecutable("fugassem_impute_nan")],
		targets = [final_corr_file],
		args = [corr_impute_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_impute_nan")

	if not evidence_type in evidence_list:
		evidence_list[evidence_type] = final_corr_file

	return final_corr_file, final_abund_file


def preprocess_evidence (evi_file, gene_file, header, output_folder, final_evidence_file, evidence_list, evidence_type,
                      workflow, threads, time_equation, mem_equation):
	"""
	Preprocess evidence for a list of genes

	Args:
		evi_file: raw function file for one taxon.
		gene_file: gene list file
		header: whether raw function file includes header
		output_folder (string): The path of the output folder.
		final_evidence_file: finalized function file for one taxon.
		evidence_list: a dictionary recording files of evidences
		evidence_type: type of evidence
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw function file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined function file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import preprocess_evidence

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		final_evidence_file = preprocess_evidence (
						evi_file,
						gene_file,
						header,
						output_dir,
						final_evidence_file,
						evidence_list,
						evidence_type,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start preprocess_evidence module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	final_evidence_file = os.path.join(main_folder, os.path.basename(final_evidence_file))
	evidence_log = final_evidence_file + ".prep_evidence.log"

	workflow.add_task (
		"fugassem_extract_feature_subset -i [depends[0]] -l [depends[1]] -t [args[0]] -o [targets[0]] > [args[1]] 2>&1",
		depends = [evi_file, gene_file, TrackedExecutable("fugassem_extract_feature_subset")],
		targets = [final_evidence_file],
		args = [header, evidence_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_extract_feature_subset")

	if not evidence_type in evidence_list:
		evidence_list[evidence_type] = final_evidence_file

	return final_evidence_file


def preprocess_function (func_file, gene_file, header, output_folder, final_func_file,
                      workflow, threads, time_equation, mem_equation):
	"""
	Preprocess function for a list of genes

	Args:
		func_file: raw function file for one taxon.
		gene_file: gene list file
		header: whether raw function file includes header
		output_folder (string): The path of the output folder.
		final_func_file: finalized function file for one taxon.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw function file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined function file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import preprocess_function

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		func_file = preprocess_function (
						func_file,
						gene_file,
						header,
						output_dir,
						final_func_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start preprocess_function module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	final_func_file = os.path.join(main_folder, os.path.basename(final_func_file))
	func_log = final_func_file + ".prep_func.log"

	workflow.add_task (
		"fugassem_extract_feature_subset -i [depends[0]] -l [depends[1]] -t [args[0]] -o [targets[0]] > [args[1]] 2>&1",
		depends = [func_file, gene_file, TrackedExecutable("fugassem_extract_feature_subset")],
		targets = [final_func_file],
		args = [header, func_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_extract_feature_subset")

	return final_func_file


def refine_function (func_file, gene_file, header, go_level, func_type, go_obo, output_folder, final_func_file, final_func_smp_file, final_funclist_file,
                     workflow, threads, time_equation, mem_equation):
	"""
	Refine function for a list of genes

	Args:
		func_file: raw function file for one taxon.
		gene_file: gene list file
		header: whether raw function file includes header
		go_level: trimming option, only terms that are informative at a given level
				[number OR fraction of genes]: Default: 50
				[none]: skip trimming
				[all]: keep all terms
		func_type: function category, e.g. BP | CC | MF
		go_obo: go-basic obo file
		output_folder (string): The path of the output folder.
		final_func_file: finalized function file for one taxon.
		final_func_file: finalized simplified function file for one taxon.
		final_funclist_file: list file of functions, e.g. BP_function_list.txt
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw function file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined function file.
		string: the file of refined simplified function file.
		string: list file of functions

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import refine_function

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		final_func_file, final_func_smp_file, final_funclist_file = refine_function (
						func_file,
						gene_file,
						header = "no,
						go_level = 50,
						func_type = "BP",
						go_obo,
						output_dir,
						final_func_file,
						final_func_smp_file,
						final_funclist_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start refine_function module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	func_info_file = re.sub(".tsv$", "." + func_type + "_list.tsv", func_file)
	final_func_file = os.path.join(main_folder, os.path.basename(final_func_file))
	final_func_smp_file = os.path.join(main_folder, os.path.basename(final_func_smp_file))
	final_funclist_file = os.path.join(main_folder, os.path.basename(final_funclist_file))
	func_log1 = final_func_file + ".geneontology.log"
	func_log2 = final_func_file + ".refined_func.log"

	if not go_obo or go_level == "none":
		workflow.add_task("ln -s [depends[0]] [targets[0]]",
		                  depends = [func_file],
		                  targets = [final_func_file],
		                  cores = 1,
		                  name = "lm__final_function_file")
		workflow.add_task("ln -s [depends[0]] [targets[0]]",
		                  depends = [func_file],
		                  targets = [final_func_smp_file],
		                  cores = 1,
		                  name = "lm__final_simplified_function_file")
	else:
		if go_level == "all":
			workflow.add_task (
				"fugassem_geneontology [depends[0]] --mapping [depends[1]] --namespace [args[0]] --outfile [targets[0]] > [args[1]] 2>&1",
				depends = [go_obo, func_file, TrackedExecutable("fugassem_geneontology")],
				targets = [func_info_file],
				args = [func_type, func_log1],
				cores = threads,
				time = time_equation,
				mem = mem_equation,
				name = "fugassem_geneontology")
		else:
			workflow.add_task (
				"fugassem_geneontology [depends[0]] --mapping [depends[1]] --namespace [args[0]] --informative [args[1]] --outfile [targets[0]] > [args[2]] 2>&1",
				depends = [go_obo, func_file, TrackedExecutable("fugassem_geneontology")],
				targets = [func_info_file],
				args = [func_type, go_level, func_log1],
				cores = threads,
				time = time_equation,
				mem = mem_equation,
				name = "fugassem_geneontology")

		workflow.add_task(
			"fugassem_refine_anns -a [depends[0]] -l [depends[1]] -t [args[0]] -o [targets[0]] > [args[1]] 2>&1",
			depends = [func_file, func_info_file, TrackedExecutable("fugassem_refine_anns")],
			targets =[final_func_file, final_func_smp_file],
			args = [func_type, func_log2],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "fugassem_refine_anns")

	workflow.add_task(
		'sed -e \'/familyID\\t/d\' [depends[0]] | awk -F \"\\t\" \'{print $2}\' | sort | uniq > [targets[0]]',
		depends = [final_func_smp_file],
		targets = [final_funclist_file],
		cores = 1,
		name = "collect__function_list")

	return final_func_file, final_func_smp_file, final_funclist_file


def preprocess_feature (raw_file, new_file, header, feature_type, output_folder, final_feature_file,
                      workflow, threads, time_equation, mem_equation):
	"""
	Preprocess features for machine learning

	Args:
		func_file: raw feature file.
		new_file: new feature file for combining.
		header: whether the new-feature file has header or not [Default: None].
		feature_type: data type of new feature, choices = ["matrix", "homology", "vector", "function"]
		output_folder (string): The path of the output folder.
		final_feature_file: finalized feature file.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw function file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined function file.

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import preprocess_feature

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		final_feaure_file = preprocess_feature (
						raw_file,
						new_file,
						header,
						feature_type,
						output_dir,
						final_feature_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start preprocess_feature module #####")

	# prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)
	final_feature_file = os.path.join(main_folder, os.path.basename(final_feature_file))
	feature_log = final_feature_file + ".prep_feature.log"

	workflow.add_task (
		"fugassem_collect_features -r [depends[0]] -n [depends[1]] -f [args[0]] -t [args[1]] -o [targets[0]] > [args[2]] 2>&1",
		depends = [raw_file, new_file, TrackedExecutable("fugassem_collect_features")],
		targets = [final_feature_file],
		args = [feature_type, header, feature_log],
		cores = 1,
		time = time_equation,
		mem = mem_equation,
		name = "fugassem_collect_features")

	return final_feature_file


def preprocessing_task (abund_file, gene_file, func_file, go_level, func_type, go_obo,
                min_prev, min_abund, min_detected, zero_flt, taxa_abund_file, corr_method,
                vector_list, matrix_list,
                output_folder, basename, feature_list, feature_list_file,
                workflow, threads, time_equation, mem_equation):
	"""
	Preprocess annotated function and evidences for a list of genes per taxon

	Args:
		abund_file: raw abundance file for one taxon.
		gene_file: gene list file for one taxon.
		func_file: raw function file for one taxon.
		go_level: trimming option, only terms that are informative at a given level
				[number OR fraction of genes]: Default: 50
				[none]: skip trimming
				[all]: keep all terms
		func_type: function category, e.g. BP | CC | MF
		go_obo: go-basic obo file
		min_prev: the minimum prevalence per gene [ Default: 0.01 ].
		min_abund: the minimum detected abundance for each gene [ Default: 0 ].
		min_detected: the minimum detected value for each covariate feature [ Default: 0 ]
		zero_flt: pre-filtering approach [ Default: None ], choices: ["lenient", "semi-strict", "strict", "None"]
		taxa_abund_file: species abundance table used for technical-zero filtering [ Default: None ]
		corr_method: method for corrleation, e.g. Pearson_SE | Pearson | Spearman | Kendall
		vector_list: [string] comma separated list of vector-based evidence files
		marix_list: [string] comma separated list of matrix-based evidence files
		output_folder (string): The path of the output folder.
		feature_list: a directory of feature files.
		feature_list_file: a list file of features files.
		basename: basename for output files.
		workflow (anadama2.workflow): An instance of the workflow class.
		threads (int): The number of threads/cores for clustering to use.
		time_equation (int): requred number of hours defined in the workflow.
		mem_equation (int): requred number of GB defined in the workflow.

	Requires:
		raw abundance file for one taxon
		raw function file for one taxon
		raw gene list file for one taxon

	Returns:
		string: the file of refined function file.
		string: the file of refined simplified function file.
		string: list file of functions
		dict: a dictionary of finalized feature file for machine learning

	Example:
		from anadama2 import Workflow
		from fugassem.tasks.preprocessing import preprocessing_task

		# create an anadama2 workflow instance
		workflow=Workflow()

		# add preprocess_function tasks
		feature_list, final_func_file, final_func_smp_file, final_funclist_file = preprocessing_task (
						abund_file,
						gene_file,
						func_file,
						go_level = 50,
						func_type = "BP",
						go_obo,
						min_prev = 0.01,
						min_abund = 0,
						min_detected = 0,
						zero_flt = None,
						taxa_abund_file = None,
						corr_method = "Pearson_SE",
						vector_list,
						matrix_list,
						output_dir,
						basename,
						feature_list,
						feature_list_file,
						workflow,
                        args.threads,
                        args.time_equation,
                        args.mem_equation)
		# run the workflow
		workflow.go()
	"""

	config.logger.info("###### Start preprocessing_task module #####")

	## prep I/O files
	main_folder = output_folder
	if not os.path.isdir(main_folder):
		os.system("mkdir -p " + main_folder)

	# abundance
	final_gene_file = os.path.join(main_folder, basename + ".final_genelist.txt")
	refined_abund_file = os.path.join(main_folder, basename + ".refined_abunds.tsv")
	final_abund_file = os.path.join(main_folder, basename + ".final_abunds.tsv")
	final_corr_file = os.path.join(main_folder, basename + ".corr_abunds.tsv")
	final_family_file = os.path.join(main_folder, basename + ".proteinfamilies.tsv")

	# annotated function
	final_func_file = os.path.join(main_folder, basename + ".final_funcs.tsv")
	final_func_smp_file = os.path.join(main_folder, basename + ".final_funcs.simple.tsv")
	final_funclist_file = os.path.join(main_folder, basename + ".final_funclist.txt")

	evidence_list = {}
	## refine abundance
	refine_abundance (abund_file, gene_file, min_prev, min_abund, min_detected, zero_flt,
	                taxa_abund_file, main_folder, refined_abund_file, final_gene_file, final_family_file,
	                workflow, threads, time_equation, mem_equation)

	## calculate co-variation
	calculate_covariation (refined_abund_file, corr_method, main_folder, final_corr_file, final_abund_file,
	                      evidence_list, "exp",
	                      workflow, threads, time_equation, mem_equation)
	#evidence_type = "exp"
	#if not evidence_type in feature_list:
	#	feature_list[evidence_type] = final_corr_file

	## refine annotated function
	refine_function (func_file, final_gene_file, "no", go_level, func_type, go_obo, main_folder,
	                final_func_file, final_func_smp_file, final_funclist_file,
	                workflow, threads, time_equation, mem_equation)

	## collect expression evidence
	evidence_type = "coexp"
	final_feature_file = os.path.join(main_folder, basename + "." + evidence_type + ".features.tsv")
	final_feature_file_trans = re.sub(".tsv$", ".trans.tsv", final_feature_file)
	preprocess_feature(final_corr_file, final_func_smp_file, None, "function", main_folder, final_feature_file,
	                   workflow, threads, time_equation, mem_equation)
	workflow.add_task(
			"fugassem_transpose < [depends[0]] > [targets[0]]",
			depends = [final_feature_file, TrackedExecutable("fugassem_transpose")],
			targets = [final_feature_file_trans],
			cores = 1,
			time = time_equation,
			mem = mem_equation,
			name = "fugassem_transpose")
	if not evidence_type in feature_list:
		feature_list[evidence_type] = final_feature_file_trans

	## vector-based evidence files (gene-over-function type of file)
	mynum = 0
	if vector_list and vector_list != "None":
		vectors = vector_list.split(",")
		for vector_file in vectors:
			vector_file = os.path.abspath(vector_file)
			mynum = mynum + 1
			evidence_type = "vector" + str(mynum)
			config.logger.info("Preprocess vector evidence: " + evidence_type + "\t" + os.path.basename(vector_file))
			final_evidence_file = os.path.join(main_folder, basename + "." + evidence_type + ".tsv")
			preprocess_evidence (vector_file, final_gene_file, "yes", main_folder, final_evidence_file,
		                 evidence_list, evidence_type,
		                 workflow, threads, time_equation, mem_equation)

			final_feature_file = os.path.join(main_folder, basename + "." + evidence_type + ".features.tsv")
			tmp_feature_file = final_feature_file + ".tmp"
			final_feature_file_trans = re.sub(".tsv$", ".trans.tsv", final_feature_file)
			preprocess_feature (final_family_file, final_evidence_file, "yes", "vector", main_folder, tmp_feature_file,
			                   workflow, threads, time_equation, mem_equation)
			preprocess_feature (tmp_feature_file, final_func_smp_file, None, "function", main_folder, final_feature_file,
			                   workflow, threads, time_equation, mem_equation)
			workflow.add_task(
				"fugassem_transpose < [depends[0]] > [targets[0]]",
				depends = [final_feature_file, TrackedExecutable("fugassem_transpose")],
				targets = [final_feature_file_trans],
				cores = 1,
				time = time_equation,
				mem = mem_equation,
				name = "fugassem_transpose")
			if not evidence_type in feature_list:
				feature_list[evidence_type] = final_feature_file_trans
		# foreach vector file
	# if exists vector list

	## matrix-based evidence files (gene-by-gene network type of file)
	mynum = 0
	if matrix_list and matrix_list != "None":
		matrixs = matrix_list.split(",")
		for matrix_file in matrixs:
			matrix_file = os.path.abspath(matrix_file)
			mynum = mynum + 1
			evidence_type = "matrix" + str(mynum)
			config.logger.info("Preprocess matrix evidence: " + evidence_type + "\t" + os.path.basename(matrix_file))
			final_evidence_file = os.path.join(main_folder, basename + "." + evidence_type + ".tsv")
			preprocess_evidence (matrix_file, final_gene_file, "no", main_folder, final_evidence_file,
		                        evidence_list, evidence_type,
		                        workflow, threads, time_equation, mem_equation)

			final_evidence_matrix = os.path.join(main_folder, basename + "." + evidence_type + ".matrix.tsv")
			convert_log = final_evidence_matrix + ".log"
			workflow.add_task(
				"fugassem_convert_coann -i [depends[0]] -s [args[0]] -f [args[1]] -t [args[2]] -o [targets[0]] > [args[3]] 2>&1",
				depends = [final_evidence_file, TrackedExecutable("fugassem_convert_coann")],
				targets = [final_evidence_matrix],
				args = [evidence_type, "vector", "no", convert_log],
				cores = 1,
				time = time_equation,
				mem = mem_equation,
				name = "fugassem_convert_coann")

			final_feature_file = os.path.join(main_folder, basename + "." + evidence_type + ".features.tsv")
			tmp_feature_file = final_feature_file + ".tmp"
			final_feature_file_trans = re.sub(".tsv$", ".trans.tsv", final_feature_file)
			preprocess_feature(final_family_file, final_evidence_matrix, "yes", "matrix", main_folder, tmp_feature_file,
		                workflow, threads, time_equation, mem_equation)
			preprocess_feature(tmp_feature_file, final_func_smp_file, None, "function", main_folder, final_feature_file,
		                workflow, threads, time_equation, mem_equation)
			workflow.add_task(
				"fugassem_transpose < [depends[0]] > [targets[0]]",
				depends = [final_feature_file, TrackedExecutable("fugassem_transpose")],
				targets = [final_feature_file_trans],
				cores = 1,
				time = time_equation,
				mem = mem_equation,
				name = "fugassem_transpose")
			if not evidence_type in feature_list:
				feature_list[evidence_type] = final_feature_file_trans
		# foreach matrix file
	# if exists matrix list

	# write evidence-features files into a list file
	utilities.dict_to_file(feature_list, feature_list_file)
	
	return feature_list, final_func_smp_file, final_func_file, final_funclist_file
