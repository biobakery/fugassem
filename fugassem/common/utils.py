#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import csv
import gzip
import bz2
from collections import defaultdict, OrderedDict
from textwrap import fill

#csv.field_size_limit( sys.maxsize )

# ---------------------------------------------------------------
# warnings
# ---------------------------------------------------------------

def warn( *args ):
    script = "?"
    if sys.argv[0] != "":
        script = os.path.split( sys.argv[0] )[1].upper()
    args = ["WARNING ({}):".format( script )] + list( args )
    print( " ".join( map( str, args ) ), file=sys.stderr )

def say( *args, **kwargs ):
    destination = kwargs.get( "file", sys.stderr )
    print( " ".join( map( str, args ) ), file=destination )

def die( *args ):
    args = ["LETHAL ERROR:"] + list( args )
    say( *args )
    sys.exit( "EXITING." )

# ---------------------------------------------------------------
# file i/o
# ---------------------------------------------------------------

def try_open( path, mode="r", *args, **kwargs ):
    """ open a (possibly compressed?) file; fail gracefully """
    fh = None
    try:
        # 1) load as gzip file
        if path.endswith( ".gz" ):
            say( "Treating", path, "as gzip file" )
            # python 2/3 switching
            if sys.version_info.major == 3:
                opener = gzip.open
                mode = "rt" if mode == "r" else mode
            else:
                opener = gzip.GzipFile
            fh = opener( path, mode=mode, *args, **kwargs )
        # 2) load as bz2 file
        elif path.endswith( ".bz2" ):
            say( "Treating", path, "as bzip2 file" )
            # python 2/3 switching
            if sys.version_info.major == 3:
                opener = bz2.open
                mode = "rt" if mode == "r" else mode
            else:
                opener = bz2.BZ2File
            fh = opener( path, mode=mode, *args, **kwargs )
        # 3) load as regular file
        else:
            fh = open( path, mode=mode, *args, **kwargs )
    except:
        die( "Problem opening", path )
    return fh

def which( program ):
    """
    Adapted from:
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    ret = None
    def is_exe( fpath ):
        return os.path.isfile( fpath ) and os.access( fpath, os.X_OK )
    fpath, fname = os.path.split( program )
    if fpath and is_exe( program ):
        ret = program
    else:
        for path in os.environ["PATH"].split( os.pathsep ):
            path = path.strip( '"' )
            exe_file = os.path.join( path, program )
            if is_exe( exe_file ):
                ret = exe_file
    return ret

def iter_lines( path, skip=0, verbose=False ):
    """ easy file loading """
    with try_open( path ) as fh:
        counter = 0
        for line in fh:
            counter += 1
            if verbose and counter % ( 1e5 ) == 0:
                say( path, "lines processed (millions): {:.1f}".format( counter / 1e6 ) )
            if counter > skip:
                yield line.rstrip( )

def reader( fh ):
    """ my favorite options for csv reader """
    for row in csv.reader( fh, csv.excel_tab ):
        yield row

def iter_rows( path, skip=0, verbose=True ):
    """ easy table loading """
    lens = set( )
    with try_open( path ) as fh:
        counter = 0
        for row in reader( fh ):
            counter += 1
            if verbose and counter % ( 1e5 ) == 0:
                say( "<{}> lines processed (millions): {:.1f}".format( path, counter / 1e6 ) )
            if counter > skip:
                lens.add( len( row ) )
                yield row
    if len( lens ) != 1 and verbose:
        warn( "unusual row lengths:", lens )

def iter_rowdicts( path, verbose=True ):
    """ easy table loading """
    lens = []
    with try_open( path ) as fh:
        counter = 0
        headers = None
        for row in reader( fh ):
            counter += 1
            if verbose and counter % ( 1e5 ) == 0:
                say( "<{}> lines processed (millions): {:.1f}".format( path, counter / 1e6 ) )
            if headers is None:
                headers = row
            elif len( row ) == len( headers ):
                rowdict = {k:v for k, v in zip( headers, row )}               
                yield rowdict
            elif verbose:
                warn( "header/row mistmatch at line", counter )

# ---------------------------------------------------------------
# text manipulation
# ---------------------------------------------------------------

def smartwrap( text, charlim ):
    """ this function was already available in python... """
    return fill( text, width=charlim )

def path2name( path ):
    """ given '/blah/blah/blah/hello.txt.gz' returns 'hello'"""
    return os.path.split( path )[1].split( "." )[0]

def rebase( path, newext=None, newdir=None ):
    olddir, name = os.path.split( path )
    if newext is not None:
        name = name.split( "." )[0:-1]
        name = ".".join( name + [newext] )
    newdir = newdir if newdir is not None else olddir
    return os.path.join( newdir, name )

def tprint( *args, **kwargs ):
    """ coerce list of items to strings then print with tabs between """
    fh = kwargs.get( "file", sys.stdout )
    print( "\t".join( map( str, args ) ), file=fh )

def write_rowdict( headers, values=None, fh=None ):
    if values is None:
        tprint( *headers, file=fh )
    else:
        row = []
        for k in headers:
            if k not in values:
                die( "Bad rowdict: couldn't find:", k, "in", values.keys( ) )
            row.append( values[k] )
        tprint( *row, file=fh )

def col2list( filename, index=0, limit=None, func=None, headers=False ):
    """ quickly load a column from a file """
    ret = []
    counter = 0
    for row in iter_rows( filename ):
        if headers:
            headers = False
        else:
            ret.append( row[index] )
            counter += 1
        if limit is not None and counter >= limit:
            break
    if func is not None:
        ret = [func( k ) for k in ret]
    return ret

def qw( multiline_string, as_dict=False, delim="\s+", comment="#" ):
    """~perl's qw->list function"""
    ret = [k for k in multiline_string.split( "\n" ) if k != "" and k[0] != comment]
    if as_dict:
        d = OrderedDict( )
        for line in ret:
            items = re.split( delim, line )
            if len( items ) != 2:
                die( "bad qw as_dict line:", line )
            else:
                d[items[0]] = items[1]
        ret = d
    return ret

def shorten( string, n=10, dummy="[...]" ):
    t = 2 * n + len( dummy )
    if len( string ) >= t:
        string = string[0:n] + dummy + string[-n:]
    return string
    
# ---------------------------------------------------------------
# dictionary methods
# ---------------------------------------------------------------

def sorteditems( dictData, reverse=False, limit=None ):
    """ return k, v pairs in v-sorted order """
    iCounter = 0
    for key in sorted( dictData.keys(), key=lambda x: dictData[x], reverse=reverse ):
        iCounter += 1
        if limit is not None and iCounter >= limit:
            break
        else:
            yield ( key, dictData[key] )

def sortedby( iterable, sorter, reverse=False ):
    if len( iterable ) != len( sorter ):
        die( "sortedby on iterables of non-equal lengths" )
    pairs = [[s, i] for s, i in zip( sorter, iterable )]
    pairs.sort( reverse=reverse )
    return [pair[1] for pair in pairs]
        
def autodict( iDepth=None, funcDefault=None ):
    """ 
    Acts as a constructor; makes an avdict 
    Can terminate at a specified depth as a defaultdict with the specified constructor.
    Example1: x = funcAVD( 3, int ) allows x["foo"]["bar"]["net"] += 1 ( i.e., a counter )
    Example2: x = funcAVD( 2, list ) allows x["foo"]["bar"].append( "net" ) 
    """
    if iDepth is None:
        return defaultdict( lambda: autodict( None ) )
    elif iDepth >= 2:
        return defaultdict( lambda: autodict( iDepth - 1, funcDefault ) )
    elif iDepth == 1:
        return defaultdict( funcDefault )

# ---------------------------------------------------------------
# classes
# ---------------------------------------------------------------

class Ticker( ):
    def __init__( self, extent, step=None, pad=2 ):
        self.count = 0
        if type( extent ) in [float, int]:
            self.total = int( extent )
        else:
            self.total = len( extent )
        self.step = step if step is not None else int( self.total / 1000 )
        self.step = max( self.step, 1 )
        self.pad = " " * pad
    def tick( self ):
        self.count += 1
        if self.count % self.step == 0:
            self.report( )
    def report( self ):
        frac = self.count / float( self.total )
        sys.stderr.write( self.pad+"%.1f%%\r" % ( 100 * frac ) )

# ---------------------------------------------------------------
# tests
# ---------------------------------------------------------------

if __name__ == "__main__":
    pass
