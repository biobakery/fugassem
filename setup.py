"""
FUGAsseM setup

To run: python setup.py install
"""

import os
import sys

# check for either of the required versions
# required python versions (3.7+)
required_python_version_major = [3]
required_python_version_minor = [7]
pass_check = False
try:
	for major, minor in zip(required_python_version_major, required_python_version_minor):
		if (sys.version_info[0] == major and sys.version_info[1] >= minor):
			pass_check = True
except (AttributeError, IndexError):
	sys.exit("CRITICAL ERROR: The python version found (version 1) " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")

if not pass_check:
	sys.exit("CRITICAL ERROR: The python version found (version " +
	         str(sys.version_info[0]) + "." + str(sys.version_info[1]) + ") " +
	         "does not match the version required (version " +
	         str(required_python_version_major) + "." +
	         str(required_python_version_minor) + "+)")


try:
	import setuptools
except ImportError:
	sys.exit("Please install setuptools.")

# check setuptools version
required_setuptools_version_major = 1
try:
	setuptools_version = setuptools.__version__
	setuptools_version_major = int(setuptools_version.split(".")[0])
	if setuptools_version_major < required_setuptools_version_major:
		sys.exit("CRITICAL ERROR: The setuptools version found (version " +
		         setuptools_version + ") does not match the version required " +
		         "(version " + str(required_setuptools_version_major) + "+)." +
					" Please upgrade your setuptools version.")
except (ValueError, IndexError, NameError):
	sys.exit("CRITICAL ERROR: Unable to call setuptools version. Please upgrade setuptools.")

from setuptools.command.install import install as _install

import distutils

# try to import urllib.request.urlretrieve for python3
try:
	from urllib.request import urlretrieve
except ImportError:
	from urllib import urlretrieve

from glob import glob
import re
import tarfile
import subprocess
import shutil
import zipfile
import tempfile
import re
import time


VERSION = "0.3.4"
AUTHOR = "FUGAsseM Development Team"
MAINTAINER = "Yancong Zhang"
MAINTAINER_EMAIL = "zhangyc201211@gmail.com"


def byte_to_megabyte(byte):
	"""
	Convert byte value to megabyte
	"""

	return byte / (1024.0 ** 2)


class ReportHook():
	def __init__(self):
		self.start_time = time.time()

	def report(self, blocknum, block_size, total_size):
		"""
		Print download progress message
		"""

		if blocknum == 0:
			self.start_time = time.time()
			if total_size > 0:
				print("Downloading file of size: " + "{:.2f}".format(byte_to_megabyte(total_size)) + " MB\n")
		else:
			total_downloaded = blocknum * block_size
			status = "{:3.2f} MB ".format(byte_to_megabyte(total_downloaded))

			if total_size > 0:
				percent_downloaded = total_downloaded * 100.0 / total_size
				# use carriage return plus sys.stdout to overwrite stdout
				try:
					download_rate = total_downloaded / (time.time() - self.start_time)
					estimated_time = (total_size - total_downloaded) / download_rate
				except ZeroDivisionError:
					download_rate = 0
					estimated_time = 0
				estimated_minutes = int(estimated_time / 60.0)
				estimated_seconds = estimated_time - estimated_minutes * 60.0
				status += "{:3.2f}".format(percent_downloaded) + " %  " + \
				          "{:5.2f}".format(byte_to_megabyte(download_rate)) + " MB/sec " + \
				          "{:2.0f}".format(estimated_minutes) + " min " + \
				          "{:2.0f}".format(estimated_seconds) + " sec "
			status += "        \r"
			sys.stdout.write(status)


def download(url, download_file):
	"""
	Download a file from a url
	"""

	try:
		print("Downloading " + url)
		file, headers = urlretrieve(url, download_file, reporthook=ReportHook().report)
		# print final return to start new line of stdout
		print("\n")
	except EnvironmentError:
		print("WARNING: Unable to download " + url)


def download_unpack_tar(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		tarfile_handle = tarfile.open(download_file)
		tarfile_handle.extractall(path=folder)
		tarfile_handle.close()
	except (EnvironmentError, tarfile.ReadError):
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


def download_unpack_zip(url, download_file_name, folder, software_name):
	"""
	Download the url to the file and decompress into the folder
	"""

	# Check for write permission to the target folder
	if not os.access(folder, os.W_OK):
		print("WARNING: The directory is not writeable: " +
		      folder + " . Please modify the permissions.")

	download_file = os.path.join(folder, download_file_name)

	download(url, download_file)

	error_during_extract = False

	try:
		zipfile_handle = zipfile.ZipFile(download_file)
		zipfile_handle.extractall(path=folder)
		zipfile_handle.close()
	except EnvironmentError:
		print("WARNING: Unable to extract " + software_name + ".")
		error_during_extract = True

	if not error_during_extract:
		try:
			os.unlink(download_file)
		except EnvironmentError:
			print("WARNING: Unable to remove the temp download: " + download_file)


setuptools.setup(
	name="fugassem",
	author=AUTHOR,
	author_email=MAINTAINER_EMAIL,
	version=VERSION,
	license="MIT",
	description="FUGAsseM: Function predictor of Uncharacterized Gene products by Assessing high-dimensional community data in microbiome",
	long_description="FUGAsseM is a computational tool based on a “guild by association” approach to predict functions of novel gene products in the context of microbial communities. "
	                 "It uses machine learning methods to predict functions of microbial proteins by integrating multiple types of community-based data, "
	                 "including co-expression patterns from metatranscriptomics, co-localization patterns from metagenomic assemblies, homologous-based annotation, domain-based annotations, etc.",
	url="https://github.com/biobakery/fugassem",
	keywords=['microbial', 'microbiome', 'bioinformatics', 'microbiology', 'metagenomic', 'metatranscriptomic',
	         'function prediction', 'unknown proteins', 'microbial dark matter',
	          'fugassem', 'anadama2'],
	platforms=['Linux', 'MacOS'],
	classifiers=[
		"Programming Language :: Python",
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Operating System :: MacOS",
		"Operating System :: Unix",
		"Programming Language :: Python :: 3.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	packages=setuptools.find_packages(),
	entry_points={
		'console_scripts': [
			'fugassem = fugassem.fugassem:main',
			'fugassem_main = fugassem.workflows.fugassem_main:main',
			'fugassem_process = fugassem.tasks.fugassem_process:main',
			'fugassem_calculate_correlation = fugassem.preprocess.calculate_correlation:main',
			'fugassem_collect_features = fugassem.preprocess.collect_features:main',
			'fugassem_convert_DDI = fugassem.preprocess.convert_DDI:main',
			'fugassem_convert_coann = fugassem.preprocess.convert_coann:main',
			'fugassem_impute_nan = fugassem.preprocess.impute_nan:main',
			'fugassem_refine_abunds = fugassem.preprocess.refine_abunds:main',
			'fugassem_refine_anns = fugassem.preprocess.refine_anns:main',
			'fugassem_split_taxa = fugassem.preprocess.split_taxa:main',
			'fugassem_build_union_terms = fugassem.preprocess.build_union_terms:main',

			'fugassem_collect_ml_results = fugassem.predict.collect_ml_results:main',
			'fugassem_machine_learning= fugassem.predict.machine_learning:main',
			'fugassem_predict_function = fugassem.predict.predict_function:main',
			'fugassem_prepare_nextlayer_ml = fugassem.predict.prepare_nextlayer_ml:main',
			'fugassem_select_feature = fugassem.predict.select_feature:main',
			'fugassem_merged_prediction = fugassem.predict.merged_prediction:main',

			'fugassem_generate_stratified_mtx_input = fugassem.tools.generate_stratified_mtx_input:main',
			'fugassem_generate_annotation_input= fugassem.tools.generate_annotation_input:main',
			'fugassem_abundance_RPK = fugassem.tools.abundance_RPK:main',
			'fugassem_abundance_normalization = fugassem.tools.abundance_normalization:main',
			'fugassem_abundance_smoothing = fugassem.tools.abundance_smoothing:main',
			'fugassem_collect_seqSimilarity_results = fugassem.tools.collect_seqSimilarity_results:main',
			'fugassem_gene_abundance = fugassem.tools.gene_abundance:main',
			'fugassem_gene_abundance_indexRef = fugassem.tools.gene_abundance_indexRef:main',
			'fugassem_gene_catalog_abundance = fugassem.tools.gene_catalog_abundance:main',
			'fugassem_prepare_annotation = fugassem.tools.prepare_annotation:main',
			'fugassem_prepare_contig = fugassem.tools.prepare_contig:main',
			'fugassem_prepare_seqSimilarity = fugassem.tools.prepare_seqSimilarity:main',
			'fugassem_prepare_coann_pairs = fugassem.tools.prepare_coann_pairs:main',
			'fugassem_stratified_abundance = fugassem.tools.stratified_abundance:main',
			'fugassem_sum_to_family_abundance = fugassem.tools.sum_to_family_abundance:main',
			'fugassem_performance_vis = fugassem.tools.performance_vis:main',

			'fugassem_extract_feature_subset = fugassem.common.extract_feature_subset:main',
			'fugassem_format_function = fugassem.common.format_function:main',
			'fugassem_geneontology = fugassem.common.geneontology:main',
			'fugassem_transpose = fugassem.common.transpose:main',
		]},
	package_data={
		'fugassem': [
			'workflows/*.py',
			'tasks/*.py',
			'tools/generate_*.py',
			'data/GO/*'
		]},
	scripts=glob('fugassem/workflows/*.py'),
	zip_safe=False
)
