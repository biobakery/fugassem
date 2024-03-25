#!/usr/bin/env python

"""
Train and predict for each function using machine learner

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
from math import sqrt
import logging
import math
from itertools import repeat
import multiprocessing as mp
#from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
import joblib
import numpy as np
from numpy import array as a
from numpy import mean, std, log10
from numpy.random import shuffle
from sklearn.ensemble import RandomForestClassifier as rfc
from sklearn.tree import DecisionTreeClassifier as dtc
from sklearn.naive_bayes import GaussianNB as nbc
from sklearn.model_selection import StratifiedKFold as skf
from sklearn.impute import SimpleImputer
from sklearn.model_selection import RandomizedSearchCV as rsc

try:
	from fugassem.common import table2
	from fugassem.common import utils 
	from fugassem import config, utilities
	from fugassem.predict import select_feature
except:
	sys.exit ("fugassem is not installed!")


# ---------------------------------------------------------------
# Description and arguments
# ---------------------------------------------------------------
description = """
Train and predict machine learner for function prediction
"""

def get_args():
	parser = argparse.ArgumentParser(
		description = "Process ML for function prediction\n",
		prog = "machine_learning")
	parser.add_argument("--feature", '-f',
	                    help='[REQUIRED] input feature file',
	                    required = True)
	parser.add_argument("--list", '-l',
	                    help='[REQUIRED] input list file of functions',
	                    required = True)
	parser.add_argument("--type", '-t',
	                    help='[REQUIRED] specify type of ML method, [RF]: Random Forest, [NB] Naive Bayes, [DT] Decision Tree',
	                    choices = ["RF", "NB", "DT"],
	                    default = "RF")
	parser.add_argument("--vector",
	                    help='[OPTIONAL] add vector from metadata to features',
	                    action = "store_true")
	parser.add_argument("--suffix", '-s',
	                    help='[OPTIONAL] suffix name of the vector feature',
	                    default = None)
	parser.add_argument("--perm",
	                    help='[OPTIONAL] do permutation analysis',
	                    action = "store_true")
	parser.add_argument("--redundancy", "-r",
						help = "[OPTIONAL] maximum correlation between features for filtering more redundant features in training set [ Default: None ]\n",
						default = None)
	parser.add_argument("--method", "-m",
						help = "[OPTIONAL] correlation methods, [ Default: Pearson ]\n",
						choices = ["Pearson", "Spearman", "Kendall", "Pearson_SE"],
						default = "Pearson")
	parser.add_argument( "--core", "-c",
						help = "[OPTIONAL] number of threads, [ Default: 1 ]\n",
						default = 1)
	parser.add_argument("--importance", "-i",
						help = "[OPTIONAL] minimum prediction importance of features from 1st round prediction for filtering less importance features in training set [ Default: None ]\n",
						default = None)
	parser.add_argument("--hyperpara", "-x",
	                    help='[OPTIONAL] do hyperparameter tuning for Random Forest',
	                    action = "store_true")
	parser.add_argument("--imppath", "-p",
	                    help = "[OPTIONAL] path for prediction importance files. It should be customized when '--importance' parameter is used [ Default: None ]\n",
						default = None)
	parser.add_argument("--outdir", '-o',
	                    help='[REQUIRED] output directary name',
	                    required = True)
	values = parser.parse_args()

	return values


# get_args


# ---------------------------------------------------------------
# prepare data
# ---------------------------------------------------------------
def collect_list(list_file):
	list = {}
	open_file = open(list_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		info[0] = re.sub(":", "_", info[0])
		list[info[0]] = ""
	open_file.close()

	return list


def prepare_data(p_features, perm, func_name, vector_name):
	# feature loading
	features = table2.table(p_features)
	try:
		funcs = features.rowdict(func_name)
	except:
		return None, None, None, None	
	features.delete_row(func_name)
	if vector_name:
		myvec_table = features.grep("ID", vector_name, in_place = False)
	else:
		myvec_table = None
	features.head("Feature", invert=True)
	if myvec_table:
		features.augment(myvec_table)
	features.float()

	#debug
	#config.logger.info (features.rowheads)

	# load function info
	funcs = {k: ("yes" if v == "yes" else "no") for k, v in funcs.items()}
	try:
		funcs = {k: ("yes" if int(v) == 1 else "no") for k, v in funcs.items()}
	except:
		pass

	# convert values / labels to numpy format
	X1 = a([a(col) for colhead, col in features.iter_cols()])
	y1 = a([funcs[x] for x in features.colheads])

	if perm:
		shuffle(y1)

	return features, funcs, X1, y1


# ---------------------------------------------------------------
# cross validate
# ---------------------------------------------------------------
def hyperpara_tuning (rf, X_train, y_train, cores):
	"""
	Use hyperparameter tuning algorithms that help to fine-tune the Machine learning models
	"""
	# number of trees
	n_estimators = [int(x) for x in np.linspace(start=100, stop=2000, num=10)]
	# number of features to consider at every split
	max_features = ['sqrt', 'log2', None]
	# maximum number of levels in tree
	max_depth = [int(x) for x in np.linspace(10, 110, num=11)]
	max_depth.append(None)
	# minimum number of samples required to split a node
	#min_samples_split = [2, 5, 10]
	# minimum number of samples required at each leaf node
	#min_samples_leaf = [1, 2, 4]
	# method of selecting samples for training each tree
	#bootstrap = [True, False]
	# create the random grid
	param_grid = {'n_estimators': n_estimators,
				  'max_features': max_features,
				  'max_depth': max_depth}

	# random search of parameters, using five fold cross validation, search across 100 different combinations, and use specified cores
	min_num = 5
	y_members = {}
	for i in y_train:
		if not i in y_members:
			y_members[i] = 1
		else:
			y_members[i] + 1
	y_members = min([y_members[i] for i in y_members.keys()])
	if y_members < min_num:
		min_num = y_members 
	if min_num < 2:
		return rf

	random_search = rsc (estimator = rf, param_distributions = param_grid, cv = min_num)
	random_search.fit(X_train, y_train)
	best_estimator = random_search.best_estimator_

	# Update the model
	updated_rf = rfc (n_estimators = best_estimator["n_estimators"],
					  max_features = best_estimator["max_features"],
					  max_depth =  best_estimator["max_depth"])

	return updated_rf


def fetch_genes (index, features):
	""" collect gene names """
	genes = {features.colheads[i]: i for i in index}
	
	return genes


def balance_train(train, features, funcs):
	""" balance class labels in the training set """
	genes = {features.colheads[i]: i for i in train}
	is_a = {}
	for s in genes:
		is_a.setdefault(funcs[s], []).append(s)
	bneck = min([len(v) for v in is_a.values()])
	train2 = []
	for v in is_a.values():
		random.shuffle(v)
		train2 += [genes[s] for s in v[0:bneck]]

	return a(train2)


def weight_train(labels):
	""" generates genes weights to handle classes of different count;
	doesn't appear to work as expected """
	counts = {k: 0 for k in labels}
	for k in labels:
		counts[k] += 1
	w = a([1 / float(counts[k]) for k in labels])
	return w


def select_nonredunt_feature (features, train, func_name, redu_level, corr_method, cores):
	""" Select non-redundant features in the training set;
	 remove redundant features features"""

	genes = {features.colheads[i]: i for i in train}
	myfeature = {}
	for i in features.rowheads:
		myv = features.row(i)
		myfeature[i] = "\t".join([str(myv[j]) for j in train])

	# select non-redundant features
	config.logger.info (func_name + ": raw size (before non-redundant feature selection) " + str(features.size()))
	select_features = select_feature.flt_redundant_feature(redu_level, myfeature, genes, corr_method, cores)
	features = features.select("ID", select_features.keys(), in_place=False)
	config.logger.info (func_name + ": refined size (after non-redundant feature selection) " + str(features.size()))

	return features


def select_impt_feature(features, train, func_name, imp_level, imp_file):
	""" Select important features in the training set;
	 remove non-important features based on the previous prediction importance reporting"""
	
	config.logger.info('select_impt_feature')

	genes = {features.colheads[i]: i for i in train}
	myfeature = {}
	for i in features.rowheads:
		myv = features.row(i)
		myfeature[i] = "\t".join([str(myv[j]) for j in train])

	# collect importance
	flt = {}
	if not os.path.isfile(imp_file):
		config.logger("Error! The prediction importance file doesn't exit: " + imp_file)
		return features
	open_file = open(imp_file, "r")
	for line in open_file:
		line = line.strip()
		if not len(line):
			continue
		info = line.split("\t")
		myf = info[0]
		myv = info[-1]
		try:
			myv = float(myv)
			if math.isnan(myv):
				flt[myf] = ""
			if myv <= float(imp_level):
				flt[myf] = ""
		except:
			flt[myf] = ""
	open_file.close()


	# select important features
	config.logger.info (func_name + ": raw size (before important feature selection) " + str(features.size()))
	select_features = {}
	for myf in myfeature.keys():
		if myf in flt:
			continue
		select_features[myf] = myfeature[myf]
	features = features.select("ID", select_features.keys(), in_place=False)
	config.logger.info (func_name + ": refined size (after important feature selection) " + str(features.size()))

	return features


def learning(ml_type, func_name, features, funcs, X1, y1, redu_level, corr_method, cores, imp_level, imp_file, hyper):
	conf = {}
	for e in y1:
		for o in y1:
			conf[(e, o)] = 0

	roc = []
	imp = []
	x_num = 0
	myskf = skf(utilities.c_folds)
	for train, test in myskf.split(X1, y1):
		x_num += 1
		myfeatures = features
		#r = rfc(n_estimators=utilities.c_trees, random_state=utilities.c_rseed)
		#config.logger.info ({k: list(y1[train]).count(k) for k in set(y1[train])})
		config.logger.info ("Fold #" + str(x_num))
		config.logger.info ("Size of training set: " + str(len(train)))
		config.logger.info ("Size of testing set: " + str(len(test)))

		# select features
		s_flag = 0
		if imp_level:
			try:
				imp_level = float(imp_level)
				myfeatures = select_impt_feature(myfeatures, train, func_name, imp_level, imp_file)
				s_flag = 1
			except ValueError:
				sys.exit("Error! Please use valid value for minimum prediction importance by '--importance' parameter")
		if redu_level:
			try:
				redu_level = float(redu_level)
				cores = int(cores)
				myfeatures = select_nonredunt_feature(myfeatures, train, func_name, redu_level, corr_method, cores)
				s_flag = 1
			except ValueError:
				sys.exit("Error! Please use valid value for maximum redundant level by '--redundancy' parameter")
		if s_flag == 1:
			# convert values / labels to numpy format
			X1 = a([a(col) for colhead, col in myfeatures.iter_cols()])
			y1 = a([funcs[x] for x in myfeatures.colheads])
		
		# Final features
		config.logger.info (func_name + ": final size (after feature selection) " + str(myfeatures.size()))

		# select tree size
		mysize = utilities.c_trees * (int(len(myfeatures.rowheads) / 1000) + 1)
		if ml_type == "RF":
			config.logger.info (func_name + ": selected number of trees " + str(mysize))
			r = rfc(n_estimators=mysize, random_state=utilities.c_rseed)
		if ml_type == "NB":
			r = nbc()
		if ml_type == "DT":
			r = dtc(random_state=utilities.c_rseed)

		# with balancing
		train = balance_train(train, myfeatures, funcs)
		if hyper and ml_type == "RF":
			try:
				r = hyperpara_tuning (r, X1[train], y1[train], cores)
			except:
				config.logger.info ("Error to run RandomizedSearchCV to tune parameters")
		try:
			r.fit( X1[train], y1[train] )
		except:
			config.logger.info ("Error for fitting model: " + func_name)
			continue
		config.logger.info ("Label summary in training set:")
		config.logger.info ({k: list(y1[train]).count(k) for k in set(y1[train])})


		# output info
		for e, o in zip(y1[test], r.predict(X1[test])):
			conf[(e, o)] += 1
		mygenes = fetch_genes (test, myfeatures)
		for myg, actual, weights in zip(mygenes, y1[test], r.predict_proba(X1[test])):
			outline = [myg, actual]
			for c, w in zip(r.classes_, weights):
				outline += [c, w]
			roc.append(map(str, outline))

		if ml_type != "NB":
			#x_num += 1
			items = [(i, n) for i, n in zip(r.feature_importances_, myfeatures.rowheads)]
			items.sort(reverse=True)
			for i, n in items:
				imp.append(str(x_num) + "\t" + "\t".join([n, str(i)]))
	
	return conf, roc, imp


# ---------------------------------------------------------------
# Write out learning results
# ---------------------------------------------------------------
def write_cross_validation(conf, roc, imp, outdir, prefix):
	myfile = os.path.join(outdir, prefix + ".xval.conf.tsv")
	with open(myfile, "w") as fh:
		for (e, o) in sorted(conf):
			count = conf[(e, o)]
			fh.write("{}\t{}\t{}".format(e, o, count) + "\n")

	myfile = os.path.join(outdir, prefix + ".xval.roc.tsv")
	with open(myfile, "w") as fh:
		for items in roc:
			fh.write("\t".join(items) + "\n")

	if len(imp) == 0:
		return None
	myfile = os.path.join(outdir, prefix + ".xval.impt.tsv")
	with open(myfile, "w") as fh:
		fh.write("fold\tfeature\timportance\n")
		for items in imp:
			fh.write(items + "\n")


# ---------------------------------------------------------------
# write out feature importance from the full model
# ---------------------------------------------------------------
def write_importance_results(ml_type, features, funcs, X1, y1, outdir, prefix, cores, hyper):
	# run full model
	if ml_type == "RF":
		mysize = utilities.c_trees * (int(len(features.rowheads) / 1000) + 1)
		r = rfc(n_estimators=mysize, random_state=utilities.c_rseed)
	if ml_type == "DT":
		r = dtc(random_state=utilities.c_rseed)

	# with balancing
	index = a(balance_train(range(len(y1)), features, funcs))
	if hyper and ml_type == "RF":
		try:
			r = hyperpara_tuning(r, X1[index], y1[index], cores)
		except:
			config.logger.info ("Error to run RandomizedSearchCV to tune parameters")
	try:
		r.fit(X1[index], y1[index])
	except:
		config.logger.info ("Error for fitting model for overall functions")
		return None

	# save model
	mymodel = os.path.join(outdir, prefix + ".random_forest.joblib")
	joblib.dump(r, mymodel)

	# output
	myfile = os.path.join(outdir, prefix + ".importance.tsv")
	with open(myfile, "w") as fh:
		items = [(i, n) for i, n in zip(r.feature_importances_, features.rowheads)]
		items.sort(reverse=True)
		for i, n in items:
			fh.write("\t".join([n, str(i)]) + "\n")


# ---------------------------------------------------------------
# process each function
# ---------------------------------------------------------------
def process_function (ml_type, myfunc, vector, suffix, imppath, p_features, perm, redu_level, corr_method, cores, imp_level, hyper, outdir):
	myout = re.sub(":", "_", myfunc)
	config.logger.info("----" + myfunc + "----")
	config.logger.info("Prepare data for learning")
	if vector:
		myvector = myfunc
	else:
		myvector = None
	if suffix:
		if suffix == "all":
			myvector =  myfunc + "__"
	if imppath:
		imp_file = os.path.join(imppath, myfunc + ".importance.tsv")
	else:
		imp_file = None
	features, funcs, X1, y1 = prepare_data(p_features, perm, myfunc, myvector)
	if not features:
		config.logger.info("WARNING: skip this function since no data available! " + myfunc)
		return None
	config.logger.info("Learning by cross validation")
	conf, roc, imp = learning(ml_type, myfunc, features, funcs, X1, y1, redu_level, corr_method, cores, imp_level, imp_file, hyper)
	config.logger.info("Write out cross-validation results")
	write_cross_validation(conf, roc, imp, outdir, myout)
	if ml_type != "NB":
		config.logger.info("Write out feature importance from the full model")
		write_importance_results(ml_type, features, funcs, X1, y1, outdir, myout, cores, hyper)


# ==============================================================
###########  Main processing ############
# ==============================================================
def main():
	### get arguments ###
	args = get_args()
	if not args.perm:
		random.seed(utilities.c_rseed)
	if args.hyperpara:
		hyper = True
	else:
		hyper = False

	p_features = args.feature
	ml_type = args.type
	outdir = args.outdir
	perm = args.perm
	vector = args.vector
	suffix = args.suffix
	if args.imppath:
		imppath = os.path.abspath(args.imppath)
	else:
		imppath = None
	redu_level = args.redundancy
	imp_level = args.importance
	corr_method = args.method
	cores = args.core
	func_list = collect_list(args.list)

	config.logger.info("Start ML process")
	config.logger.info("ML type: " + ml_type)
	if not os.path.isdir(outdir):
		os.system("mkdir -p " + outdir)
	try:
		cores = int(cores)
	except:
		cores = 1
	pool = ThreadPool(cores)
	pool.starmap(process_function, zip(repeat(ml_type), func_list.keys(), repeat(vector), repeat(suffix), repeat(imppath),
	                                   repeat(p_features), repeat(perm), repeat(redu_level), repeat(corr_method),
	                                   repeat(cores), repeat(imp_level), repeat(hyper), repeat(outdir)))
	pool.close()
	pool.join()
	
	# collect predictions
	myfile = os.path.join(outdir, "prediction_list.txt")
	open_file = open(myfile, "w")
	files = [os.path.abspath(x) for x in os.listdir(outdir)]
	for x in files:
		if re.search(".xval.roc.tsv", x):
			open_file.write(x + "\n")
	open_file.close()

	config.logger.info("Successfully finish ML process")

# end: main

if __name__ == '__main__':
	main()
