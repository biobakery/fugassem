#!/usr/bin/env python
##########################################################################
# Function: Calculate evaluation measures and visualize results
# Author: Yancong Zhang (zhangyc201211@gmail.com)
# Date: 05/03/2022
##########################################################################
import sys
import os
import os.path
import re
import argparse
import logging
import math
import numpy as np

from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf

try:
	from fugassem import utilities
	from fugassem import config
except:
	sys.exit("fugassem is not installed!")	


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
def get_args ():
	"""
	Parse the arguments from the user
	"""

	parser = argparse.ArgumentParser(
		description = "Calculate evaluation measures and visualize results\n",
		prog = "performance_vis")
	parser.add_argument(
		"-i", "--input",
		help = "[REQUIRED] input prediction file",
		required = True)
	parser.add_argument(
		"-t", "--type",
		help = "[OPTIONAL] type of true label",
		choices = ["raw_ann", "GT_ann", "new_ann"],
		default = "raw_ann")
	parser.add_argument(
		"-p", "--plot",
		help = "[OPTIONAL] plotting types",
		choices = ["auc", "aupr", "all"],
		default = "auc")
	parser.add_argument(
		"-o", "--output",
		help = "[REQUIRED] prefix name of output files",
		required = True)

	return parser.parse_args()


def collect_prediction (pred_file, label_type):
	"""
	Collect prediction info
	Input:
		pred_file - prediction file
	Output: dictionary of prediction info
	"""

	config.logger.info('collect_prediction')

	if not os.path.isfile(pred_file):
		config.logger.info ("ERROR! Please provide valid prediction file: " + pred_file)
		return 0
	stat_probs = {}
	stat_label = {}
	overall_probs = {}
	overall_label = {}
	head = 0
	titles = {}
	flag = 0
	for line in utilities.gzip_bzip2_biom_open_readlines(pred_file):
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		if flag == 0:
			flag = 1
			for i in info:
				titles[i] = info.index(i)
			continue
		myf = info[titles["func"]]
		myc = info[titles["category"]]
		myc = re.sub("UniRef90_GO_", "", myc)
		myp = info[titles["score"]]
		try:
			myp = float(myp)
		except:
			myp = float("NaN")
		if label_type in titles:
			myl = info[titles[label_type]]
			if myl == "yes" or myl == "1":
				myl = 1
			if myl == "no" or myl == "0":
				myl = 0
		else:
			config.logger.info ("Error! Label column doesn't exist: " + label_type)
			return None, None, None, None
		myid = myf + "\t" + myc
		if not myid in stat_probs:
			stat_probs[myid] = []
		stat_probs[myid].append(myp)
		if not myid in stat_label:
			stat_label[myid] = []
		stat_label[myid].append(myl)
		if not myc in overall_probs:
			overall_probs[myc] = []
		overall_probs[myc].append(myp)
		if not myc in overall_label:
			overall_label[myc] = []
		overall_label[myc].append(myl)

	return stat_probs, stat_label, overall_probs, overall_label


def calculate_measure (testy, probsy):
	"""
	Calculate measurements
	Input:
		testy - the true outcomes (0,1)
		probsy -  the predicted probabilities for the 1 class
	Output: measures
	"""

	config.logger.info('calculate_measure')

	myauc = roc_auc_score(testy, probsy)
	fpr, tpr, thresholds1 = roc_curve(testy, probsy)
	precision, recall, thresholds2 = precision_recall_curve(testy, probsy)
	aupr = auc (recall, precision)
	try:
		f1 = f1_score(testy, probsy)
	except:
		f1 = []
		i = 0
		while i < len(precision):
			x = precision[i]
			y = recall[i]
			try:
				if (x+y) != 0:
					z = 2 * (x * y) / (x + y)
				else:
					z = float("NaN")
			except:
				z = float("NaN")
			f1.append(z)
			i = i + 1
		f1 = np.array(f1)

	return thresholds1, tpr, fpr, myauc, thresholds2, precision, recall, f1, aupr


def quantify_prediction (stat_probs, stat_label, overall_probs, overall_label, output):
	"""
	Process function assignment process
	Input:
		stat_probs - prediction values per function
		stat_label - true labels per function
		overall_probs - overall prediction values
		overall_label - overall labels
		output - prefix name of output files
	Output: quantitative evaluations
	"""
	config.logger.info('quantify_prediction')
	
	# each function
	stat_each_auc = {}
	stat_each_aupr = {}
	outfile = output + ".stat.tsv"
	open_out = open(outfile, "w")
	open_out.write("func\tfunc_type\tthresholds\ttpr\tfpr\tauc\tauc_label\tprecision\trecall\tF_measure\taupr\n")
	for myid in sorted(stat_probs):
		if myid in stat_label:
			testy = stat_label[myid]
			probsy = stat_probs[myid]
			lnum = np.unique(np.array(testy))
			if len(lnum) < 2:
				continue
			thresholds1, tpr, fpr, auc, thresholds2, precision, recall, f1, aupr = calculate_measure (testy, probsy)
			stat_each_auc[myid] = [fpr, tpr, auc]
			stat_each_aupr[myid] = [recall, precision, aupr]

			threds = [x for x in thresholds1 if x in thresholds2]
			#print(myid)
			#print(str(thresholds1.size) + "\t" + str(thresholds2.size) + "\t" + str(len(threds)))
			for item in threds:
				myi = np.where (thresholds1 == item)
				mytpr = tpr[myi][0]
				myfpr = fpr[myi][0]
				myauc = auc
				myauc_label = "AUC = " + "{:.3f}".format(myauc)
				myj = np.where (thresholds2 == item)
				myp = precision[myj][0]
				myr = recall[myj][0]
				myf1 = f1[myj][0]
				myaupr = aupr
				mystr = myid + "\t" + str(item) + "\t" + str(mytpr) + "\t" + str(myfpr) + "\t" + str(myauc) + "\t" + str(myauc_label) + "\t" + str(myp) + "\t" + str(myr) + "\t" + str(myf1) + "\t" + str(myaupr)
				open_out.write(mystr + "\n")
	open_out.close()

	# overall
	stat_all_auc = {}
	stat_all_aupr = {}
	outfile = output + ".overall.stat.tsv"
	open_out = open(outfile, "w")
	open_out.write("func_type\tthresholds\ttpr\tfpr\tauc\tauc_label\tprecision\trecall\tF_measure\taupr\n")
	for myid in sorted(overall_probs):
		if myid in overall_label:
			testy = overall_label[myid]
			probsy = overall_probs[myid]
			lnum = np.unique(np.array(testy))
			if len(lnum) < 2:
				continue
			thresholds1, tpr, fpr, auc, thresholds2, precision, recall, f1, aupr = calculate_measure(testy, probsy)
			myid_new = myid
			stat_all_auc[myid_new] = [fpr, tpr, auc]
			stat_all_aupr[myid_new] = [recall, precision, aupr]

			threds = [x for x in thresholds1 if x in thresholds2]
			for item in threds:
				myi = np.where (thresholds1 == item)
				mytpr = tpr[myi][0]
				myfpr = fpr[myi][0]
				myauc = auc
				myauc_label = "AUC = " + "{:.3f}".format(myauc)
				myj = np.where (thresholds2 == item)
				myp = precision[myj][0]
				myr = recall[myj][0]
				myf1 = f1[myj][0]
				myaupr = aupr
				mystr = myid + "\t" + str(item) + "\t" + str(mytpr) + "\t" + str(myfpr) + "\t" + str(myauc) + "\t" + str(myauc_label) + "\t" + str(myp) + "\t" + str(myr) + "\t" + str(myf1) + "\t" + str(myaupr)
				open_out.write(mystr + "\n")
	open_out.close()

	return stat_each_auc, stat_each_aupr, stat_all_auc, stat_all_aupr


def vis_measures (measures, mylabel, xlab, ylab, outfile):
	"""
	Visualize quantified measures
	Input:
		measures - quantified measures
		label - plot title
		outfile - output file name
	Output: plot quntatitve measures
	"""
	config.logger.info('vis_measures')

	plt.rcParams.update({'figure.max_open_warning': 0})
	lw = 2
	figs = list()
	for myid in sorted(measures.keys()):
		myt = myid.split("\t")[-1]
		myid_new = myid.split("\t")[0]
		fpr, tpr, auc = measures[myid]
		auc = "{:.3f}".format(auc)
		fig = plt.figure()
		fig.set_size_inches(5, 5)
		color = "blue"
		if myt == "BP":
			color = "blue"
		if myt == "MF":
			color = "red"
		if myt == "CC":
			color = "orange"
		label = mylabel + str(auc) + ")"
		plt.plot(fpr, tpr, label = label, color = color)
		plt.plot([0, 1], [0, 1], color = "black", linestyle="--")
		plt.xlim([0.0, 1.0])
		plt.ylim([0.0, 1.05])
		plt.xlabel(xlab)
		plt.ylabel(ylab)
		plt.title(myid_new)
		plt.legend(loc="lower right")
		figs.append(fig)

	mypdf = pdf.PdfPages(outfile)
	for fig in figs:  ## will open an empty extra figure :(
		mypdf.savefig(fig)
	mypdf.close()


#==============================================================
###########  Main processing ############
#==============================================================
def main():
	
	### get arguments ###
	values = get_args ()

	config.logger.info("Start performance_vis......")

	stat_probs, stat_label, overall_probs, overall_label = collect_prediction (values.input, values.type)
	stat_each_auc, stat_each_aupr, stat_all_auc, stat_all_aupr = quantify_prediction (stat_probs, stat_label, overall_probs, overall_label, values.output)

	if values.plot == "auc" or values.plot == "all":
		outfile1 = values.output + ".each.auc.pdf"
		#outfile2 = values.output + ".overall.auc.pdf"
		mylabel = "ROC curve (area = "
		xlab = "False Positive Rate"
		ylab  = "True Positive Rate"
		vis_measures(stat_each_auc, mylabel, xlab, ylab, outfile1)
		#vis_measures(stat_all_auc, mylabel, xlab, ylab, outfile2)
	if values.plot == "aupr" or values.plot == "all":
		outfile1 = values.output + ".each.aupr.pdf"
		#outfile2 = values.output + ".overall.aupr.pdf"
		mylabel = "Precision-Recall curve (area = "
		xlab = "Recall"
		ylab  = "Precision"
		vis_measures(stat_each_aupr, mylabel, xlab, ylab, outfile1)
		#vis_measures(stat_all_aupr, mylabel, xlab, ylab, outfile2)

	config.logger.info("Finish performance_vis")

# end: main

if __name__ == '__main__':
	main()
