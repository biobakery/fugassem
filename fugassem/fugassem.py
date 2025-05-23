#!/usr/bin/env python

"""
FUGAsseM: Function predictor of Uncharacterized Gene products by Assessing high-dimensional community data in microbiome


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
import argparse
import subprocess

# import the workflow class from anadama2
from anadama2 import Workflow

VERSION = "0.3.8"

WORKFLOW_FOLDER = "workflows"
WORKFLOW_EXTENSION = ".py"
WORKFLOW_DEFAULT = "fugassem_workflow"


def find_workflows():
	""" Search for installed workflows """

	workflow_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), WORKFLOW_FOLDER)
	workflows = {}
	for file in os.listdir(workflow_folder):
		# look for files with the expected extension
		if file.endswith(WORKFLOW_EXTENSION):
			# do not need to add full path as these are also installed as executable scripts
			workflows[file.replace(WORKFLOW_EXTENSION, "")] = file

	return workflows


def parse_arguments(args, workflows):
	"""
	Parse the arguments from the user
	"""
	parser = argparse.ArgumentParser(
		description = "FUGAsseM: Function predictor of uncharacterized gene products in communities\n",
		formatter_class = argparse.ArgumentDefaultsHelpFormatter,
		prog = "fugassem")
	parser.add_argument(
		"--version",
		action = "version",
		version = "%(prog)s v"+ VERSION)	
	parser.add_argument(
		"workflow",
		choices = workflows,
		default = WORKFLOW_DEFAULT)
	
	return parser.parse_args(args)


def run_workflow(args, workflow):
	""" Run the workflow with the arguments provided """

	try:
		command = [workflow] + args[2:]
		subprocess.call(command)
	except (subprocess.CalledProcessError, EnvironmentError):
		sys.exit("Error: Unable to run workflow: " + " ".join(command))


def main():
	# find workflows
	workflows = find_workflows()

	# parse the arguments (only the first two as the rest are for the workflow)
	if len(sys.argv) < 2:
		sys.argv = [sys.argv[0], WORKFLOW_DEFAULT]
	if not sys.argv[1] in workflows:
		tmp = sys.argv[1:]
		sys.argv[1] = WORKFLOW_DEFAULT
		sys.argv = sys.argv[:2] + tmp
	args = parse_arguments(sys.argv[1:2], workflows.keys())
	
	# run the workflow (providing all of the arguments)
	run_workflow(sys.argv, workflows[args.workflow])


if __name__ == "__main__":
	main()
