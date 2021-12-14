#! /usr/bin/env python

"""
perform tests on the table2 object
( assumes all entries are numeric )
"""

# imports
import sys
from math import log
from scipy.stats import spearmanr

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_fSmallNumber   = 1e-3
c_fAllowedZeroes = 0.05

# ---------------------------------------------------------------
# utilities
# ---------------------------------------------------------------

def funcTestSum ( aValues ):
    s = sum( aValues )
    if 1 - s > c_fSmallNumber:
        sys.stderr.write ("WARNING: Sum of vector < 1:" + str(s) + "\n")
    if s - 1 > c_fSmallNumber:
        sys.stderr.write ("WARNING: Sum of vector > 1:" + str(s) + "\n")

def funcTestZeroes ( aValues ):
    f = aValues.count( 0 ) / float( len( aValues ) )
    if f > c_fAllowedZeroes:
        sys.stderr.write ("WARNING: Non-trivial vector zero fraction:" + str(f) + "\n")

# ---------------------------------------------------------------
# one sample functions
# ---------------------------------------------------------------

dictOneSampleFunctions = {}

def funcPielou ( aValues ):
    funcTestSum( aValues )
    # don't penalize zero values
    aEntropic = [-k * log( k ) for k in aValues if k > 0]
    return sum( aEntropic ) / log( len( aEntropic ) )
dictOneSampleFunctions["pielou"] = funcPielou

# ---------------------------------------------------------------
# two sample functions
# ---------------------------------------------------------------

dictTwoSampleFunctions = {}

def funcBrayCurtis ( aValues1, aValues2, verbose=True ):
    """
    Two equivalent formulas ( I implemented #1 )
    ( 1 ) 1 - 2 * sum( min( x[i], y[i] ) ) / total
    ( 2 ) sum( abs( x[i] - y[i] ) ) / total
    """
    if verbose:
        funcTestSum( aValues1 )
        funcTestSum( aValues2 )
    # if data are compositional, then bc = 1 - share
    share = sum( [min( v1, v2 ) for v1, v2 in zip( aValues1, aValues2 )] )
    total = sum( [v1 + v2 for v1, v2 in zip( aValues1, aValues2 )] )
    return 1 - ( 2 * share / float( total ) )
dictTwoSampleFunctions["bray_curtis"] = funcBrayCurtis

def funcSpearman ( aValues1, aValues2, verbose=True ):
    if verbose:
        funcTestZeroes( aValues1 )
        funcTestZeroes( aValues2 )
    return spearmanr( aValues1, aValues2 )[0]
dictTwoSampleFunctions["spearman"] = funcSpearman

def funcFilteredSpearman ( aValues1, aValues2 ):
    aDoublezero = [ True if aValues1[k] == 0 and aValues2[k] == 0 else False for k in range( len( aValues1 ) ) ]
    aValues1 = [ aValues1[k] for k in range( len( aDoublezero ) ) if not aDoublezero[k] ]
    aValues2 = [ aValues2[k] for k in range( len( aDoublezero ) ) if not aDoublezero[k] ]
    return spearmanr( aValues1, aValues2 )
dictTwoSampleFunctions["filtered_spearman"] = funcSpearman

# ---------------------------------------------------------------
# execute test, return dictionary of results
# ---------------------------------------------------------------

def one_sample ( tableData, func=None ):
    dictResults = {}
    for sample, col in tableData.iter_cols():
        dictResults[sample] = dictOneSampleFunctions[func]( col )
    return dictResults

def two_sample ( tableData, func=None ):
    dictResults = {}
    for sample1, col1 in tableData.iter_cols():
        for sample2, col2 in tableData.iter_cols():
            if sample1 < sample2:
                dictResults[( sample1, sample2 )] = dictTwoSampleFunctions[func]( col1, col2 )
    return dictResults

# ---------------------------------------------------------------
# testing
# ---------------------------------------------------------------

if __name__ == "__main__":
    x = [1,2,3,4,5]
    y = [1,2,3,4,5]
    print (x, y, funcBrayCurtis( x, y ))
    x = [0.1, 0.2, 0.3, 0.4]
    y = [0.11, 0.19, 0.23, 0.47]
    print (x, y, funcBrayCurtis( x, y ))
    x = [0.1, 0.2, 0.3, 0.4]
    y = [0.6, 0.1, 0.2, 0.1]
    print (x, y, funcBrayCurtis( x, y ))
