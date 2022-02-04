"""
FUGAsseM: config module
Configuration settings

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
import argparse
import re

# try to import the python2 ConfigParser
# if unable to import, then try to import the python3 configparser
try:
    import ConfigParser as configparser
except ImportError:
    import configparser

import logging


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """ 
Config fugassem
"""

def get_args (): 
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument('-c', "--config",
                        help = 'input global config file',
                        required = False,
						default = "none")
	
	values = parser.parse_args()
	return values
# get_args


def log_settings():
	"""
	Write to the log file the config settings for the run
	"""

	lines = []
	lines.append("DATABASE SETTINGS")
	lines.append("nucleotide database folder = " + nucleotide_database)
	lines.append("protein database folder = " + protein_database)
	if pathways_database_part1:
		lines.append("pathways database file 1 = " + pathways_database_part1)
		lines.append("pathways database file 2 = " + pathways_database_part2)
	else:
		lines.append("pathways database file = " + pathways_database_part2)
	lines.append("utility mapping database folder = " + utility_mapping_database)
	lines.append("")

	lines.append("RUN MODES")
	lines.append("resume = " + str(resume))
	lines.append("verbose = " + str(verbose))
	lines.append("bypass prescreen = " + str(bypass_prescreen))
	lines.append("bypass nucleotide index = " + str(bypass_nucleotide_index))
	lines.append("bypass nucleotide search = " + str(bypass_nucleotide_search))
	lines.append("bypass translated search = " + str(bypass_translated_search))
	lines.append("translated search = " + translated_alignment_selected)
	lines.append("pick frames = " + pick_frames_toggle)
	lines.append("threads = " + str(threads))
	lines.append("")

	lines.append("SEARCH MODE")
	lines.append("search mode = " + search_mode)
	lines.append("identity threshold = " + str(identity_threshold))
	lines.append("")

	lines.append("ALIGNMENT SETTINGS")
	lines.append("evalue threshold = " + str(evalue_threshold))
	lines.append("prescreen threshold = " + str(prescreen_threshold))
	lines.append("translated subject coverage threshold = " + str(translated_subject_coverage_threshold))
	lines.append("translated query coverage threshold = " + str(translated_query_coverage_threshold))
	lines.append("")

	lines.append("PATHWAYS SETTINGS")
	lines.append("minpath = " + minpath_toggle)
	lines.append("xipe = " + xipe_toggle)
	lines.append("gap fill = " + gap_fill_toggle)
	lines.append("")

	lines.append("INPUT AND OUTPUT FORMATS")
	lines.append("input file format = " + input_format)
	lines.append("output file format = " + output_format)
	lines.append("output max decimals = " + str(output_max_decimals))
	lines.append("remove stratified output = " + str(remove_stratified_output))
	lines.append("remove column description output = " + str(remove_column_description_output))
	lines.append("log level = " + log_level)
	lines.append("")

	logger.info("\nRun config settings: \n\n" + "\n".join(lines))


def update_user_edit_config_file_single_item(section, name, value):
	"""
	Update the settings to the user editable config file for one item
	"""

	new_config_items = {section: {name: value}}

	update_user_edit_config_file(new_config_items)

	print("MetaWIBELE configuration file updated: " + section + " : " + name + " = " + str(value))


def update_user_edit_config_file(new_config_items):
	"""
	Update the settings to the user editable config file
	"""

	config = configparser.RawConfigParser()

	# start with the current config settings
	config_items = read_user_edit_config_file(full_path_user_edit_config_file)

	# update with the new config items
	for section in new_config_items:
		for name, value in new_config_items[section].items():
			if section in config_items:
				if name in config_items[section]:
					config_items[section][name] = value
				else:
					sys.exit("ERROR: Unable to add new name ( " + name +
					         " ) to existing section ( " + section + " ) to " +
					         " config file: " + full_path_user_edit_config_file)
			else:
				sys.exit("ERROR: Unable to add new section ( " + section +
				         " ) to config file: " + full_path_user_edit_config_file)

	for section in config_items:
		config.add_section(section)
		for name, value in config_items[section].items():
			value = str(value)
			if "file" in section or "folder" in section:
				# convert to absolute path if needed
				if not os.path.isabs(value):
					value = os.path.abspath(value)
			config.set(section, name, value)

	try:
		file_handle = open(full_path_user_edit_config_file, "wt")
		config.write(file_handle)
		file_handle.close()
	except EnvironmentError:
		sys.exit("Unable to write to the MetaWIBELE config file.")


def read_user_edit_config_file(full_path_user_edit_config_file):
	"""
	Read the settings from the config file
	"""

	config = configparser.ConfigParser()

	try:
		config.read(full_path_user_edit_config_file)
	except EnvironmentError:
		sys.exit("Unable to read from the config file: " + full_path_user_edit_config_file)

	# read through all of the sections
	config_items = {}
	for section in config.sections():
		config_list = config.items(section)
		config_items[section] = {}
		for name, value in config_list:
			if "file" in section or "folder" in section:
				# if not absolute path, then return absolute path relative to this folder
				if not os.path.isabs(value):
					value = os.path.abspath(os.path.join(os.path.dirname(full_path_user_edit_config_file), value))
			config_items[section][name] = value

	return config_items


def get_item(config_items, section, name, type=None):
	"""
	Get the item from the dictionary of section/names from the user edit config file
	"""

	# try to obtain the value from the config dictionary
	try:
		value = config_items[section][name]
	except KeyError:
		sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
		         " . \nItem not found. \nItem should be in section (" + section + ") with name (" + name + ").")

	# if present, try to change the value type
	if type:
		try:
			if type == "string":
				value = str(value)
			elif type == "int":
				value = int(value)
			elif type == "float":
				value = float(value)
			elif type == "bool":
				if value in ["False", "false", "F", "f"]:
					value = False
				elif value in ["True", "true", "T", "t"]:
					value = True
				else:
					raise ValueError
		except ValueError:
			sys.exit("CRITICAL ERROR: Unable to load value from " + full_path_user_edit_config_file +
			         " . \nItem found in section (" + section + ") with name (" + name + "). " +
			         "\nItem is not of type (" + type + ").")

	return value


## default option for MetaWIBELE ##
version = '0.2.1'
log_level_choices = ["DEBUG","INFO","WARNING","ERROR","CRITICAL"]
log_level = log_level_choices[1]
verbose = 'DEBUG'

# name global logging instance
logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s: %(message)s',
					level = getattr(logging, log_level), filemode='w', datefmt='%m/%d/%Y %I:%M:%S %p')


## I/O files and extensions
# feature and annotation
c_abund_extension = ".abunds.tsv"
c_ann_extension = ".anns.tsv"
c_gene_extension = ".genelist.txt"
taxa_abund_list = "taxa_abunds_files.txt"


## User config file ##
fugassem_install_directory = os.path.dirname(os.path.abspath(__file__))
config_directory = os.path.join(fugassem_install_directory, "configs")
database_directory = os.path.join(fugassem_install_directory, "data")
basename = "fugassem"


## Databases ##
# GO database
#go_database_dir = get_item (config_items, "database", "go_db", "string")
go_database_dir = os.path.join(database_directory, "GO")
go_obo = None
go_obo_all = None
if os.path.exists (go_database_dir):
	if not go_database_dir.lower() == "none" and not go_database_dir == "":
		files = [os.path.abspath(x) for x in os.listdir(go_database_dir)]
		for i in files:
			myname = os.path.basename(i)
			if myname == "go-basic.obo":
				go_obo  = os.path.join(go_database_dir, myname)
			if myname == "go.obo.gz":
				go_obo_all  = os.path.join(go_database_dir, myname)


# abundance
abundance_detection_level = 0 
