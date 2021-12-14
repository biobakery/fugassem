#!/usr/bin/env python

"""
Slice columns from tables into dictionaries
"""

import sys
import csv
from fugassem.common.utils import try_open, warn, die

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

c_default = None

# ---------------------------------------------------------------
# shared methods
# ---------------------------------------------------------------

def get_value( row, value=None, func=None ):
    """ formats a value line for assignment to keys """
    # note: value (index) can be 0, so "if value" not safe
    if func is not None and value is not None:
        # get value, apply func, return
        return func( row[value] )
    elif value is not None:
        # just return value
        return row[value]
    elif func is not None:
        # apply function to whole row
        return func( row )
    else:
        return c_default

def read_csv( path, headers=False ):
    """ shared file reader """
    fh = try_open( path )
    if headers:
        fh.readline()
    return csv.reader( fh, dialect="excel-tab" )

def test_insert( d, key, value, verbose=True ):
    """ inserts a key, value pair but warns if we're overwriting something """
    if verbose and key in d and d[key] != value:
        warn( "Overwriting <{}>:<{}> with <{}>".format( key, d[key], value ) )
    d[key] = value

def sort_pair( a, b ):
    """ Sort a pair of keys for smart tuplekey insertion """
    return ( a, b ) if a < b else ( b, a )

def nesteddict2tupledict( dd, sort=True ):
    """
    Convert a nested dict [a][b] into a tupledict [( a, b )]
    default behavior is to perform standard sort on a, b to make keys unique
    """
    tupledict = {}
    for key1, inner in dd.items( ):
        for key2, value in inner.items( ):
            key = sort_pair( key1, key2 ) if sort else ( key1, key2 )
            test_insert( tupledict, key, value )
    return tupledict

# ---------------------------------------------------------------
# methods for loading dictionaries from files
# ---------------------------------------------------------------

def col2dict( path, key=0, value=None, 
              func=None, headers=False, verbose=True ):
    """
    From a tab-delimitted text file
    Return a nested dictionary where one columns is the KEY
    Additional VALUE column is optional
    """    
    d = {}
    for row in read_csv( path, headers ):
        test_insert( d, row[key], get_value( row, value, func ), verbose=verbose )
    return d

def col2dict2( path, key1=0, key2=1, value=None, 
               func=None, headers=False, 
               mirror=False, tupledict=False, verbose=True ):
    """
    From a tab-delimitted text file
    Return a nested dictionary where two columns are 
    the outer/inner key ( KEY1, KEY2 )
    Additional VALUE column is optional
    """
    dd = {}
    for row in read_csv( path, headers ):
        inner = dd.setdefault( row[key1], {} )
        test_insert( inner, row[key2], get_value( row, value, func ), verbose=verbose )
        if mirror:
            inner = dd.setdefault( row[key2], {} )
            test_insert( inner, row[key1], get_value( row, value, func ), verbose=verbose )
    return dd if not tupledict else nesteddict2tupledict( dd, sort=mirror )

def polymap( path, key=0, skip=None, headers=False, reverse=False, sets=False ):
    """
    input like:
    key1 val1 val2 val3
    key2 val4 val5
    output like:
    [key][value]
    """
    if skip is None:
        skip = []
    skip.append( key )
    dd = {}
    for row in read_csv( path, headers ):
        if key > len( row ) - 1:
            warn( "skipping unindexable short row:", row )
            continue
        inner = dd.setdefault( row[key], {} )
        for i, item in enumerate( row ):
            if i not in skip:
                inner[item] = c_default
    if reverse:
        dd2 = {}
        for key, inner in dd.items():
            for key2 in inner:
                dd2.setdefault( key2, {} )[key] = c_default
        dd = dd2
    if sets:
        dd = {k:set( inner ) for k, inner in dd.items( )}
    return dd

# ---------------------------------------------------------------
# other methods
# ---------------------------------------------------------------

def ddjoin( dd1, dd2 ):
    """ [k1][k2] and [k2][k3]=value yields [k1][k3]=value """
    dd = {}
    for k1, inner1 in dd1.items():
        for k2 in inner1:
            for k3, value in dd2.get( k2, {} ).items( ):
                dd.setdefault( k1, {} )[k3] = value
    return dd
