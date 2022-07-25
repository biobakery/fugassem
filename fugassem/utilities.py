#!/usr/bin/env python
import sys
import os
import os.path
import re
import argparse
import subprocess
import csv
import gzip
import bz2
import time
import math
import csv

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------
# indicator of a comment line or the header
# the last line in the file with this indicator is the header
GENE_TABLE_COMMENT_LINE = "#"

# the extension used for biom files
BIOM_FILE_EXTENSION = ".biom"

c_many_bytes = 1e8
c_zip_multiplier = 10

# taxa and family
c_str_unknown = "NO_NAME"
c_ungrouped = "UNGROUPED"
c_unmapped = "UNMAPPED"
c_msp_unknown = "msp_unknown"
c_unintegrated = "UNINTEGRATED"
c_topsort = {
	c_unmapped: 0,
	c_ungrouped: 1,
	c_unintegrated: 2,
	"UniRef50_unknown": 3,
	"UniRef90_unknown": 4,
}

PROTEIN_FAMILY_ID = "familyID"
PROTEIN_ID = "seqID"
c_metedata_delim = "."	 # nested metadata, e.g. CD.dysbiosis
c_strat_delim = "|" 	 # strantified item, e.g. Cluster_1000010|Bacteroides dorei
c_taxon_delim = "."	 	 # taxonomic lineage, e.g. g__Faecalibacterium.s__Faecalibacterium_prausnitzii.t__Faecalibacterium_prausnitzii_A2-165
c_multiname_delim = ";"	 # multiple ietms, e.g. PF00482;PF01841
c_name_delim = ": "
c_msp_unknown = "msp_unknown"

# setting for RF
c_folds = 5
c_trees = 100
c_rseed = 1201

# setting for DT
c_depth = 10

# ---------------------------------------------------------------
# utilities used for I/O data
# ---------------------------------------------------------------
def read_data_from_file(infile, header_flag="yes"):
	data = {}
	open_file = open(infile, "r")
	if header_flag == "yes":
		header = open_file.readline().strip()
	else:
		header = None
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		data[info[0]] = line
	# foreach line
	open_file.close()

	return data, header


def file_to_dict(infile):
	"""
	read data from file into dictionary variable
	"""

	data = {}
	open_file = open(infile, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		data[info[0]] = info[-1]
	return data


def dict_to_file(dict_data, outfile):
	"""
	write data in dictionary into file
	"""

	open_out = open(outfile, "w")
	for mydata in sorted(dict_data.keys()):
		myv = dict_data[mydata]
		if myv == "":
			open_out.write(mydata + "\n")
		else:
			open_out.write(mydata + "\t" + myv + "\n")
	open_out.close()


def is_file_exist(file_path):
	""" Check if file is not empty by confirming if its size is more than 0 bytes"""
	# Check if file exist and it is not empty
	status = 0
	if os.path.exists(file_path):
		if time.time() - os.stat(file_path).st_mtime > 60:
			status = 1
	return status



def find_files(folder, extension=None, exit_if_not_found=None):
	""" Return the files in the given folder with the extension if provided

	Args:
		folder (string): A path to a folder
		extension (string): The file extension to search for (optional)
		exit_if_not_found (bool): Indicator to check if files exist (optional)

	Requires:
		None

	Returns:
		list: A list of files in the folder

	Example:
		files = find_files("examples","fastq")
	"""
	# get all of the files in the folder
	files = [os.path.join(folder, file) for file in os.listdir(folder)]
	files = list(filter(lambda file: os.path.isfile(file), files))

	# filter to only files with extension
	if extension:
		files = list(filter(lambda file: file.endswith(extension), files))

	if exit_if_not_found:
		if not files:
			message = "ERROR: No files were found in the folder " + folder
			if extension:
				message += " with extension " + extension
			sys.exit(message + " .\n")

	return files


def name_files(names, folder, subfolder=None, tag=None, extension=None, create_folder=None):
	""" Return a list of file names based on the names and folders provided

	Args:
		names (list or string): A list of basenames or files.
		folder (string): The path to the folder.
		subfolder (string): The subfolder to use with the files (optional).
		tag (string): The tag to add to the file basenames (optional).
		extension (string): The extension to use for the files (optional).
		create_folder (bool): Create the folder and subfolder if they do not exist (optional).

	Requires:
		None

	Returns:
		list: A list of file names (or string if input is string).

	Example:
		files = name_files(["file1","file2"], "output")
	"""

	# if names is a list, convert to string
	was_string = False
	if isinstance(names, basestring):
		was_string = True
		names = [names]

	# get the basenames from the files
	names = [os.path.basename(name) for name in names]

	# use the full path to the folder
	folder = os.path.abspath(folder)

	# get the name of the full folder plus subfolder if provided
	if subfolder:
		folder = os.path.join(folder, subfolder)

	# add the extension if provided, and replace existing
	if extension:
		names = [os.path.splitext(name)[0] + "." + extension for name in names]

	# add the tag to the names, if provided
	if tag:
		names = [os.path.splitext(name)[0] + "_" + tag + os.path.splitext(name)[1] for name in names]

	files = [os.path.join(folder, name) for name in names]

	if create_folder:
		create_folders(os.path.dirname(files[0]))

	# if the input was originally a string, convert from list
	if was_string:
		files = files[0]

	return files


def name_task(sample, software):
	""" Name the task based on the sample name and software """

	return software + "____" + os.path.basename(sample)


def add_to_list(items, new_item):
	""" Add the value to the list/tuple. If the item is not a list, create a new
		list from the item and the value

	Args:
		items (list, string or tuple): Single or multiple items
		new_item (string): The new value

	Returns:
		(list): A list of all values
	"""

	if isinstance(items, tuple):
		items = [i for i in items]

	if not isinstance(items, list):
		items = [items]

	return items + [new_item]


def create_folders(folder):
	""" Create folder if it does not exist

	Args:
		folder (string): The full path to the folder.

	Requires:
		None

	Returns:
		None

	Example:
		create_folders("new_folder")
	"""

	try:
		if not os.path.exists(folder):
			os.makedirs(folder)
	except EnvironmentError:
		print("Warning: Unable to create folder: " + folder)


def get_files(folder, extension):
	""" Return paths to all files in a folder with a given extension """

	for file in os.listdir(folder):
		file = os.path.join(folder, file)
		if os.path.isfile(file) and file.endswith(extension):
			yield file


def find_exe_in_path(exe):
	"""
	Check that an executable exists in $PATH
	"""

	paths = os.environ["PATH"].split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path, exe)
		if os.path.exists(fullexe):
			if os.access(fullexe, os.X_OK):
				return True
	return False


def find_exe_in_path(exe):
	"""
	Check that an executable exists in $PATH
	"""

	paths = os.environ["PATH"].split(os.pathsep)
	for path in paths:
		fullexe = os.path.join(path, exe)
		if os.path.exists(fullexe):
			if os.access(fullexe, os.X_OK):
				return True
	return False


def biom_to_tsv(biom_file, new_tsv_file, taxonomy=None):
	"""
	Convert from a biom to tsv file
	"""

	cmd = ["biom", "convert", "-i", biom_file, "-o", new_tsv_file, "--to-tsv"]

	# check if taxonomy is set (can be set to zero)
	if taxonomy != None:
		cmd += ["--header-key", "taxonomy"]

	try:
		if os.path.isfile(new_tsv_file):
			os.unlink(new_tsv_file)
		p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
	except (EnvironmentError, subprocess.CalledProcessError):
		command = " ".join(cmd)
		sys.exit("Unable to convert biom file to tsv" + "\n" + command)


def tsv_to_biom(tsv_file, biom_file):
	"""
	Convert from a biom to tsv file
	"""

	cmd = ["biom", "convert", "-i", tsv_file, "-o", biom_file, "--table-type", "Gene table", "--to-hdf5"]

	try:
		if os.path.isfile(biom_file):
			os.unlink(biom_file)
		p_out = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
	except (EnvironmentError, subprocess.CalledProcessError):
		command = " ".join(cmd)
		sys.exit("Unable to convert tsv file to biom" + "\n" + command)


def process_gene_table_with_header(gene_table, allow_for_missing_header=None):
	"""
	Process through the header portion of the gene table file
	"""

	# try to open the file
	try:
		lines = gzip_bzip2_biom_open_readlines(gene_table)
	except EnvironmentError:
		sys.exit("Unable to read file: " + gene_table)

	# find the headers
	header = ""
	first_data_line = ""
	for line in lines:
		if line[0] == GENE_TABLE_COMMENT_LINE:
			header = line
		else:
			first_data_line = line
			break

	if not header and not allow_for_missing_header:
		sys.exit("File does not have a required header: " + gene_table +
		         " . Please add a header which includes the indicator: " +
		         GENE_TABLE_COMMENT_LINE)

	# provide the header, if one was found
	if header:
		yield header

	# provide the first data line
	yield first_data_line

	# now provide the remaining lines
	for line in lines:
		yield line


def sample_names (files, extension, pair_identifier = None):
	""" Return the basenames of the files, without any extensions, as the sample names

	Args:
		files (list): A list of files (with or without the full paths)
		extension (string): The extension for all files.
		pair_identifier (string): The string in the file basename to identify
			the first pair in the set (optional).

	Requires:
		None

	Returns:
		list: A list of sample names (file basenames)

	Example:
		names = sample_names(["1.R1.fq", "1.R2.fq"],".fq")

	"""

	# add period to extension if not included
	#if not extension.startswith("."):
	#	extension = "." + extension

	# if files is a string, convert to a list
	convert = False
	if isinstance(files, str):
		files = [files]
		convert = True

	samples = [os.path.basename(file).replace(extension, "") for file in files]

	# remove the pair_idenifier from the sample name, if provided
	if pair_identifier:
		# only remove the last instance of the pair identifier
		samples = [pair_identifier.join(sample.split(pair_identifier)[:-1]) if pair_identifier in sample else sample for
		           sample in samples]

	if convert:
		samples = samples[0]

	return samples


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# utilities used by the rename, renorm, run tasks, etc. scripts
# ---------------------------------------------------------------
def run_task (command, **keywords):
	""" Run the task command, formatting command with keywords. The command stdout
		and stderr are written to the workflow log.

	Args:
		command (string): A string to execute on the command line. It can be
			formatted the same as a task command.

	Returns:
		(int): Return code from command.
	"""

	from anadama2.helpers import format_command
	from anadama2.helpers import sh

	# format the command to include the items for this task
	command = format_command(command, **keywords)

	# run the command
	return_code = sh(command)()

	return return_code


# ---------------------------------------------------------------
# ---------------------------------------------------------------
# helper classes
# ---------------------------------------------------------------

class Table():
	"""
	Very basic table class; would be more efficient using numpy 2D array
	"""

	def __init__(self, path):
		self.anchor = None
		self.colheads = []
		self.rowheads = []
		self.data = []
		self.is_stratified = False
		if path is None:
			rows = csv.reader(sys.stdin, dialect='excel-tab')
			path = "STDIN"
			# print( "Loading table from: <STDIN>", file=sys.stderr )
			print("Loading table from: <STDIN>")
		else:
			rows = [line.split("\t") for line in process_gene_table_with_header(path, True)]
			# print( "Loading table from:", path, file=sys.stderr )
			print("Loading table from:" + path)
			size_warn(path)
		for row in rows:
			try:
				if self.anchor is None:
					self.anchor = row[0]
					self.colheads = row[1:]
				else:
					self.rowheads.append(row[0])
					self.data.append(row[1:])
			except IndexError:
				# ignore empty lines in input file
				pass
		for rowhead in self.rowheads:
			if c_strat_delim in rowhead:
				# print( "  Treating", path, "as stratified output, e.g.",
				#       rowhead.split( c_strat_delim ), file=sys.stderr )
				self.is_stratified = True
				break

	def write(self, path=None, unfloat=False):
		""" If the output file has the biom extension, then write a biom output file """

		# get the rows of output to write
		rows = self.write_rows(unfloat)

		# check if the output file specified is a biom file based on the extension
		try:
			biom_file = path.endswith(BIOM_FILE_EXTENSION)
		except (AttributeError, ValueError):
			biom_file = False

		if biom_file:
			write_biom(path, rows)
		else:
			write_tsv(path, rows)

	def write_rows(self, unfloat=False):
		""" Yield the values to write to the output file """
		yield [self.anchor] + self.colheads
		for i in range(len(self.rowheads)):
			values = self.data[i][:]
			if unfloat:
				values = list(map(lambda x: "%.6g" % (x), values))
			yield [self.rowheads[i]] + values


def write_tsv(path, rows):
	""" Write the output in tsv (possibly compressed) format to a file or stdout """
	fh = try_zip_open(path, write=True) if path is not None else sys.stdout
	writer = csv.writer(fh, delimiter="\t", lineterminator="\n")

	for row in rows:
		writer.writerow(row)


def write_biom(path, rows):
	""" Write the file in biom format """

	try:
		import biom
	except ImportError:
		sys.exit("Could not find the biom software." +
		         " This software is required since the input file is a biom file.")

	try:
		import numpy
	except ImportError:
		sys.exit("Could not find the numpy software." +
		         " This software is required since the input file is a biom file.")

	try:
		import h5py
	except ImportError:
		sys.exit("Could not find the h5py software." +
		         " This software is required since the input file is a biom file.")

	# reformat the rows into a biom table
	samples = next(rows)[1:]
	ids = []
	data = []
	for row in rows:
		ids.append(row[0])
		data.append(row[1:])

	table = biom.Table(numpy.array(data), ids, samples)

	# write a h5py biom table
	with h5py.File(path, 'w') as file_handle:
		table.to_hdf5(file_handle, "humann2 utility script")


class Ticker():
	def __init__(self, iterable, step=100, pad="  "):
		self.count = 0
		self.total = len(iterable)
		self.step = 100
		self.pad = pad

	def tick(self):
		self.count += 1
		if self.count % self.step == 0:
			self.report()

	def report(self):
		frac = self.count / float(self.total)
		# print( self.pad+"{:.1f}%".format( 100 * frac ), file=sys.stderr, end="\r" )
		print(self.pad + "{:.1f}%".format(100 * frac))


# ---------------------------------------------------------------
# helper functions
# ---------------------------------------------------------------

def size_warn(path):
	m = 1 if ".gz" not in path else c_zip_multiplier
	if m * os.path.getsize(path) > c_many_bytes:
		# print( "  This is a large file, one moment please...", file=sys.stderr )
		print("  This is a large file, one moment please...")


def try_zip_open(path, write=None):
	"""
	open an uncompressed or gzipped file; fail gracefully
	"""
	fh = None

	# set the open mode
	if write:
		open_mode = "w"
	elif path.endswith(".bz2"):
		open_mode = "r"
	else:
		open_mode = "rt"

	try:
		if path.endswith(".gz"):
			fh = gzip.open(path, open_mode)
		elif path.endswith(".bz2"):
			fh = bz2.BZ2File(path, open_mode)
		else:
			fh = open(path, open_mode)
	except EnvironmentError:
		sys.exit("Problem opening file: " + path)
	return fh


def read_biom_table(path):
	"""
	return the lines in the biom file
	"""

	try:
		import biom
	except ImportError:
		sys.exit("Could not find the biom software." +
		         " This software is required since the input file is a biom file.")

	try:
		tsv_table = biom.load_table(path).to_tsv().split("\n")
	except (EnvironmentError, TypeError):
		sys.exit("ERROR: Unable to read biom input file.")

	return tsv_table


def gzip_bzip2_biom_open_readlines(path):
	"""
	return the lines in the opened file for tab delimited text, gzip, bzip2 and biom files
	"""

	# if the file is biom, convert to text and return lines
	if path.endswith(BIOM_FILE_EXTENSION):
		for line in read_biom_table(path):
			yield line
	else:
		with try_zip_open(path) as file_handle:
			for line in file_handle:
				if path.endswith(".bz2"):
					# convert the line to text from binary
					yield line.decode('utf-8').rstrip()
				else:
					yield line.rstrip()


def load_polymap(path, start=0, skip=None, allowed_keys=None, allowed_values=None):
	"""
	Load a file like:
	A 1 2
	B 1
	B 3
	C 1 2 4
	To a nested dict structure:
	{A:{1:1, 2:1}, B:{1:1, 3:1}, C:{1:1, 2:2, 4:1}
	Inner values are not important (set to 1)
	"""
	polymap = {}
	# print( "Loading mapping file from:", path, file=sys.stderr )
	print("Loading mapping file from:" + path)
	size_warn(path)
	for line in gzip_bzip2_biom_open_readlines(path):
		row = line.split("\t")
		key = row[start]
		if allowed_keys is None or key in allowed_keys:
			for i, value in enumerate(row):
				if i != start and (skip is None or i not in skip):
					if allowed_values is None or value in allowed_values:
						polymap.setdefault(key, {})[value] = 1
	return polymap


def fsplit(feature):
	items = feature.split(c_strat_delim)
	stratum = None if len(items) == 1 else items[1]
	items = items[0].split(c_name_delim)
	name = None if len(items) == 1 else items[1]
	feature = items[0]
	return feature, name, stratum


def fjoin(feature, name=None, stratum=None):
	if name is not None:
		feature = c_name_delim.join([feature, name])
	if stratum is not None:
		feature = c_strat_delim.join([feature, stratum])
	return feature


def fsort(features):
	# force 1|A to come before 11
	features = sorted(features, key=lambda f: f.split(c_strat_delim))
	# force special features to the top (defined above)
	default = 1 + max(c_topsort.values())
	features = sorted(features, key=lambda f: c_topsort.get(fsplit(f)[0], default))
	return features


# ==============================================================
# utilities used for basis statistics info
# ==============================================================
def mean(data):
	"""Return the sample arithmetic mean of data."""
	n = len(data)
	if n < 1:
		raise ValueError('mean requires at least one data point')
	return sum(data) / float(n)


def _ss(data):
	"""Return sum of square deviations of sequence data."""
	c = mean(data)
	ss = sum((x - c) ** 2 for x in data)
	return ss


def stddev(data, ddof=0):
	"""Calculates the population standard deviation
	by default; specify ddof=1 to compute the sample
	standard deviation."""
	n = len(data)
	if n < 1:
		raise ValueError('variance requires at least two data points')
	ss = _ss(data)
	pvar = ss / (n - ddof)
	return pvar ** 0.5


def remove_duplicate(duplicate):
	"""
	remove duplicates in list
	"""
	final_list = []
	for num in duplicate:
		if num not in final_list:
			final_list.append(num)
	return final_list


def refine_taxon (myt):
	myt = re.sub("\.", "", myt)
	myt = re.sub("\(", "", myt)
	myt = re.sub("\)", "", myt)
	myt = re.sub("\[", "", myt)
	myt = re.sub("\]", "", myt)
	myt = re.sub("=", "", myt)
	myt = re.sub("\:", "_", myt)
	myt = re.sub("\/", "_", myt)
	myt = re.sub("\s+", "_", myt)
	myt = re.sub("___", "_", myt)
	myt = re.sub("__", "_", myt)

	return myt
