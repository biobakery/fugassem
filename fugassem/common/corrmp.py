import ctypes
from multiprocessing.sharedctypes import RawArray
import numpy as np
import pandas as pd
import multiprocessing as mp
from typing import *
from scipy.special import btdtr
from scipy.stats import linregress
import math
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from sklearn.metrics.pairwise import pairwise_distances
import math

shared_variables_dictionary = {}

ArrayLike = Union[pd.DataFrame, np.ndarray]


def remove_missing_values(x, y):
	'''
	Given x and y all in numpy arrays, remove pairs that contain missing values
	'''
	# nan != nan = TRUE
	nas = np.logical_or(x != x, y != y)
	return (x[~nas], y[~nas])


def remove_zero_values(x, y, method):
	'''
	Given x and y all in numpy arrays, remove pairs that contain zero values
	'''
	zeros = np.array([False for i in range(len(x))])
	if method == "none":
		return (x, y)
	if method == "lenient":
		if np.logical_or(sum(x) == 0, sum(y) == 0):
			return (np.array([]), np.array([]))
	if method == "semi_strict":
		zeros = np.logical_and(x == 0, y == 0)
	if method == "strict":
		zeros = np.logical_or(x == 0, y == 0)
	x = x[~zeros]
	y = y[~zeros]
	if len(x) < 2 or len(y) < 2:
		return (np.array([]), np.array([]))
	return (x, y)


def remove_zero_values_log(x, y, method, log_method="log"):
	'''
	Given x and y all in numpy arrays, remove pairs that contain zero values and log transformed
	'''
	zeros = np.array([False for i in range(len(x))])
	if method == "none" or method == "lenient":
		if np.logical_or(sum(x) == 0, sum(y) == 0):
			return (np.array([]), np.array([]))
	if method == "semi_strict":
		zeros = np.logical_and(x == 0, y == 0)
	if method == "strict":
		zeros = np.logical_or(x == 0, y == 0)
	x = x[~zeros]
	y = y[~zeros]
	if len(x) == 0 or len(y) == 0:
		return (np.array([]), np.array([]))
	if len(x[x != 0]) == 0 or len(y[y != 0]) == 0:
		return (np.array([]), np.array([]))
	s_x = min(x[x != 0]) / 2
	x = np.where(x == 0, s_x, x)
	s_y = min(y[y != 0]) / 2
	y = np.where(y == 0, s_y, y)
	if len(x) < 2 or len(y) < 2:
		return (np.array([]), np.array([]))
	if log_method == "log":
		return (np.log(x), np.log(y))
	elif log_method == "log2":
		return (np.log2(x), np.log2(y))
	elif log_method == "log10":
		return (np.log10(x), np.log10(y))
	else:
		return (np.log(x), np.log(y))


def jaccard(x):
	sim = pd.DataFrame(1 - pairwise_distances(x.T.to_numpy(), metric='jaccard'),
					   index=x.columns, columns=x.columns)
	return sim


class NaNCorrMp:

	@staticmethod
	def _init_worker(X: RawArray, X_finite_mask: RawArray, X_corr: RawArray,
					 X_shape: Tuple[int, int], X_corr_shape: Tuple[int, int],
					 X_p_value: RawArray = None,
					 X_se: RawArray = None) -> None:
		shared_variables_dictionary['X'] = X
		shared_variables_dictionary['X_finite_mask'] = X_finite_mask
		shared_variables_dictionary['X_corr'] = X_corr
		shared_variables_dictionary['X_shape'] = X_shape
		shared_variables_dictionary['X_corr_shape'] = X_corr_shape
		shared_variables_dictionary['X_p_value'] = X_p_value
		shared_variables_dictionary['X_p_value_shape'] = X_corr_shape
		shared_variables_dictionary['X_se'] = X_se
		shared_variables_dictionary['X_se_shape'] = X_corr_shape

	@staticmethod
	def calculate(X: ArrayLike, n_jobs: int = -1, chunks: int = 500, zero: str = "none", transform: str = "log", corr: str = "pearson") -> ArrayLike:
		return NaNCorrMp._calculate(X=X, n_jobs=n_jobs, chunks=chunks, add_p_values=False, zero=zero, transform=transform, corr=corr)

	@staticmethod
	def calculate_with_pvalue(X: ArrayLike, n_jobs: int = -1, chunks: int = 500, zero: str = "none", transform: str = "log", corr: str = "pearson") -> Tuple[
		ArrayLike, ArrayLike]:
		return NaNCorrMp._calculate(X=X, n_jobs=n_jobs, chunks=chunks, add_p_values=True, zero=zero, transform=transform, corr=corr)

	@staticmethod
	def _calculate(X: ArrayLike, n_jobs: int, chunks: int, add_p_values: bool, zero: str, transform: str, corr: str) -> Union[
		ArrayLike, Tuple[ArrayLike, ArrayLike]]:
		X_array = X.to_numpy(dtype=np.float64, copy=True).transpose() if type(X) == pd.DataFrame else X
		X_raw = RawArray(ctypes.c_double, X_array.shape[0] * X_array.shape[1])
		X_np = np.frombuffer(X_raw, dtype=np.float64).reshape(X_array.shape)
		np.copyto(X_np, X_array)

		finite_mask_data = np.isfinite(X_array)
		finite_mask_raw = RawArray(ctypes.c_bool, X_array.shape[0] * X_array.shape[1])
		finite_mask_np = np.frombuffer(finite_mask_raw, dtype=np.bool).reshape(X_array.shape)
		np.copyto(finite_mask_np, finite_mask_data)

		X_corr = np.ndarray(shape=(X_array.shape[0], X_array.shape[0]), dtype=np.float64)
		X_corr_raw = RawArray(ctypes.c_double, X_corr.shape[0] * X_corr.shape[1])
		X_corr_np = np.frombuffer(X_corr_raw, dtype=np.float64).reshape(X_corr.shape)

		if add_p_values:
			X_p_value = np.ndarray(shape=X_corr.shape, dtype=np.float64)
			X_p_value_raw = RawArray(ctypes.c_double, X_p_value.shape[0] * X_p_value.shape[1])
			X_p_value_np = np.frombuffer(X_p_value_raw, dtype=np.float64).reshape(X_corr.shape)
		else:
			X_p_value_np = None
			X_p_value_raw = None
			X_p_value_np = None

		arguments = ((j, i, zero, transform, add_p_values, corr) for i in range(X_array.shape[0]) for j in range(i))
		processes = n_jobs if n_jobs > 0 else None
		# worker_function = NaNCorrMp._set_correlation_with_p_value if add_p_values else NaNCorrMp._set_correlation
		worker_function = NaNCorrMp._set_correlation
		with mp.Pool(processes=processes,
					 initializer=NaNCorrMp._init_worker,
					 initargs=(X_raw, finite_mask_raw, X_corr_raw, X_np.shape, X_corr_np.shape, X_p_value_raw)) \
				as pool:
			list(pool.imap_unordered(worker_function, arguments, chunks))

		for i in range(X_corr_np.shape[0]):
			X_corr_np[i][i] = 1.0

		if add_p_values:
			if type(X) == pd.DataFrame:
				return (
					pd.DataFrame(X_corr_np, columns=X.columns, index=X.columns),
					pd.DataFrame(X_p_value_np, columns=X.columns, index=X.columns)
				)
			else:
				return X_corr_np, X_p_value_np

		if type(X) == pd.DataFrame:
			return pd.DataFrame(X_corr_np, columns=X.columns, index=X.columns)
		else:
			return X_corr_np

	@staticmethod
	def _set_correlation(arguments: Tuple[int, int, str, bool, bool, str]) -> None:
		j, i, zero, transform, add_p_values, corr_m = arguments
		if add_p_values:
			X_np, X_corr_np, finite_mask, X_p_value_np = NaNCorrMp._get_global_variables(get_p_value=True)
		else:
			X_np, X_corr_np, finite_mask = NaNCorrMp._get_global_variables()
		finites = finite_mask[i] & finite_mask[j]
		x = X_np[i][finites]
		y = X_np[j][finites]
		corr, p_value = NaNCorrMp._corr(x, y, zero, transform, corr_m)
		X_corr_np[i][j] = corr
		X_corr_np[j][i] = corr
		if add_p_values:
			X_p_value_np[i][j] = p_value
			X_p_value_np[j][i] = p_value

	@staticmethod
	def _get_global_variables(get_p_value: bool = False) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray],
																  Tuple[
																	  np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
		X_np = np.frombuffer(shared_variables_dictionary['X'], dtype=np.float64).reshape(
			shared_variables_dictionary['X_shape'])
		X_corr_np = np.frombuffer(shared_variables_dictionary['X_corr'], dtype=np.float64).reshape(
			shared_variables_dictionary['X_corr_shape'])
		finite_mask = np.frombuffer(shared_variables_dictionary['X_finite_mask'], dtype=bool).reshape(
			shared_variables_dictionary['X_shape'])
		if not get_p_value:
			return X_np, X_corr_np, finite_mask
		else:
			X_p_value_np = np.frombuffer(shared_variables_dictionary['X_p_value'], dtype=np.float64).reshape(
				shared_variables_dictionary['X_p_value_shape'])
			return X_np, X_corr_np, finite_mask, X_p_value_np

	@staticmethod
	def _corr(x: np.ndarray, y: np.ndarray, zero: str, transform: str, corr_m: str) -> Tuple[float, float, np.ndarray]:
		x, y = remove_missing_values(x, y)
		if transform != "none":
			x, y = remove_zero_values_log(x, y, zero, transform)
		else:
			x, y = remove_zero_values(x, y, zero)
		if len(x) < 5:
			return float("nan"), float("nan")
		# mx, my = x.mean(), y.mean()
		# xm, ym = x - mx, y - my
		# r_num = np.add.reduce(xm * ym)
		# r_den = np.sqrt((xm * xm).sum() * (ym * ym).sum())
		# r = r_num / r_den
		# return max(min(r, 1.0), -1.0), x
		corr = float("nan")
		pval = float("nan")
		if corr_m == "pearson":
			corr, pval = pearsonr(x, y)
		if corr_m == "spearman":
			res = spearmanr(x, y)
			corr = res[0]
			pval = res[1]
		if corr == "" or corr == None:
			corr = float("nan")
		if pval == "" or pval == None:
			pval = float("nan")
		return corr, pval

	@staticmethod
	def _set_correlation_with_p_value(arguments: Tuple[int, int, str, str, str]) -> None:
		j, i, zero, transform, corr_m = arguments
		X_np, X_corr_np, finite_mask, X_p_value_np = NaNCorrMp._get_global_variables(get_p_value=True)
		finites = finite_mask[i] & finite_mask[j]
		x = X_np[i][finites]
		y = X_np[j][finites]
		corr, p_value = NaNCorrMp._corr(x, y, zero, transform, corr_m)
		X_corr_np[i][j] = corr
		X_corr_np[j][i] = corr
		# if len(x) == 0:
		#	p_value = float("nan")
		# else:
		#	p_value = NaNCorrMp._p_value(corr, len(x))
		X_p_value_np[i][j] = p_value
		X_p_value_np[j][i] = p_value

	@staticmethod
	def _p_value(corr: float, observation_length: int) -> float:
		ab = observation_length / 2 - 1
		if ab == 0:
			p_value = 1.0
		else:
			p_value = 2 * btdtr(ab, ab, 0.5 * (1 - abs(np.float64(corr))))
		return p_value


	## new functions for correlation with corr, pvalue, SE
	@staticmethod
	def calculate_with_pvalue_se (X: ArrayLike, n_jobs: int = -1, chunks: int = 500, zero: str = "none", transform: str = "log", corr: str = "pearson") -> Tuple[
		ArrayLike, ArrayLike]:
		return NaNCorrMp._calculate_se(X=X, n_jobs=n_jobs, chunks=chunks, add_p_values=True, zero=zero, transform=transform, corr=corr)

	@staticmethod
	def _calculate_se(X: ArrayLike, n_jobs: int, chunks: int, add_p_values: bool, zero: str, transform: str, corr: str) -> Union[
		Tuple[ArrayLike, ArrayLike], Tuple[ArrayLike, ArrayLike, ArrayLike]]:
		X_array = X.to_numpy(dtype=np.float64, copy=True).transpose() if type(X) == pd.DataFrame else X
		X_raw = RawArray(ctypes.c_double, X_array.shape[0] * X_array.shape[1])
		X_np = np.frombuffer(X_raw, dtype=np.float64).reshape(X_array.shape)
		np.copyto(X_np, X_array)

		finite_mask_data = np.isfinite(X_array)
		finite_mask_raw = RawArray(ctypes.c_bool, X_array.shape[0] * X_array.shape[1])
		finite_mask_np = np.frombuffer(finite_mask_raw, dtype=np.bool).reshape(X_array.shape)
		np.copyto(finite_mask_np, finite_mask_data)

		X_corr = np.ndarray(shape=(X_array.shape[0], X_array.shape[0]), dtype=np.float64)
		X_corr_raw = RawArray(ctypes.c_double, X_corr.shape[0] * X_corr.shape[1])
		X_corr_np = np.frombuffer(X_corr_raw, dtype=np.float64).reshape(X_corr.shape)

		if add_p_values:
			X_p_value = np.ndarray(shape=X_corr.shape, dtype=np.float64)
			X_p_value_raw = RawArray(ctypes.c_double, X_p_value.shape[0] * X_p_value.shape[1])
			X_p_value_np = np.frombuffer(X_p_value_raw, dtype=np.float64).reshape(X_corr.shape)
			X_se = np.ndarray(shape=X_corr.shape, dtype=np.float64)
			X_se_raw = RawArray(ctypes.c_double, X_se.shape[0] * X_se.shape[1])
			X_se_np = np.frombuffer(X_se_raw, dtype=np.float64).reshape(X_corr.shape)
		else:
			X_p_value_np = None
			X_p_value_raw = None
			X_p_value_np = None
			X_se_np = None
			X_se_raw = None
			X_se_np = None

		arguments = ((j, i, zero, transform, add_p_values) for i in range(X_array.shape[0]) for j in range(i))
		processes = n_jobs if n_jobs > 0 else None
		worker_function = NaNCorrMp._set_correlation_with_pvalue_se if add_p_values else NaNCorrMp._set_correlation
		with mp.Pool(processes=processes,
					 initializer=NaNCorrMp._init_worker,
					 initargs=(
							 X_raw, finite_mask_raw, X_corr_raw, X_np.shape, X_corr_np.shape, X_p_value_raw, X_se_raw)) \
				as pool:
			list(pool.imap_unordered(worker_function, arguments, chunks))

		for i in range(X_corr_np.shape[0]):
			X_corr_np[i][i] = 1.0

		if add_p_values:
			if type(X) == pd.DataFrame:
				return (
					pd.DataFrame(X_corr_np, columns=X.columns, index=X.columns),
					pd.DataFrame(X_p_value_np, columns=X.columns, index=X.columns),
					pd.DataFrame(X_se_np, columns=X.columns, index=X.columns)
				)
			else:
				return X_corr_np, X_p_value_np, X_se_np

		if type(X) == pd.DataFrame:
			return (pd.DataFrame(X_corr_np, columns=X.columns, index=X.columns),
					pd.DataFrame(X_se_np, columns=X.columns, index=X.columns)
					)
		else:
			return X_corr_np, X_se_np

	@staticmethod
	def _set_correlation_with_pvalue_se(arguments: Tuple[int, int, str, str, bool]) -> None:
		j, i, zero, transform, add_p_values = arguments
		X_np, X_corr_np, finite_mask, X_p_value_np, X_se_np = NaNCorrMp._get_global_variables_se(get_p_value=True)
		finites = finite_mask[i] & finite_mask[j]
		x = X_np[i][finites]
		y = X_np[j][finites]
		corr, p_value, std_err = NaNCorrMp._corr_linregress(x, y, zero, transform)
		X_corr_np[i][j] = corr
		X_corr_np[j][i] = corr
		X_p_value_np[i][j] = p_value
		X_p_value_np[j][i] = p_value
		X_se_np[i][j] = std_err
		X_se_np[j][i] = std_err

	@staticmethod
	def _get_global_variables_se(get_p_value: bool = False) -> Union[Tuple[np.ndarray, np.ndarray, np.ndarray],
																	 Tuple[
																		 np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
		X_np = np.frombuffer(shared_variables_dictionary['X'], dtype=np.float64).reshape(
			shared_variables_dictionary['X_shape'])
		X_corr_np = np.frombuffer(shared_variables_dictionary['X_corr'], dtype=np.float64).reshape(
			shared_variables_dictionary['X_corr_shape'])
		finite_mask = np.frombuffer(shared_variables_dictionary['X_finite_mask'], dtype=bool).reshape(
			shared_variables_dictionary['X_shape'])
		if not get_p_value:
			return X_np, X_corr_np, finite_mask
		else:
			X_p_value_np = np.frombuffer(shared_variables_dictionary['X_p_value'], dtype=np.float64).reshape(
				shared_variables_dictionary['X_p_value_shape'])
			X_se_np = np.frombuffer(shared_variables_dictionary['X_se'], dtype=np.float64).reshape(
				shared_variables_dictionary['X_se_shape'])
			return X_np, X_corr_np, finite_mask, X_p_value_np, X_se_np

	@staticmethod
	def _corr_linregress(x: np.ndarray, y: np.ndarray, zero: str, transform: str) -> float:
		x, y = remove_missing_values(x, y)
		if transform != "none":
			x, y = remove_zero_values_log(x, y, zero, transform)
		else:
			x, y = remove_zero_values(x, y, zero)
		if len(x) < 5:
			return float("nan"), float("nan"), float("nan")
		a = x.tolist()
		b = y.tolist()
		try:
			mystat = linregress(a, b)
			myslope = mystat[0]
			myintercept = mystat[1]
			myr = max(min(mystat[2], 1.0), -1.0)
			myp = mystat[3]
			myse = mystat[4]
		except:
			myr = float("nan")
			myp = float("nan")
			myse = float("nan")
		if myr == "" or myr == None:
			myr = float("nan")
		if myp == "" or myp == None:
			myp = float("nan")
		if myse == "" or myse == None:
			myse = float("nan")
		
		return myr, myp, myse
