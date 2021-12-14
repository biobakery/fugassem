#! /usr/bin/env python

# ---------------------------------------------------------------
# imports 
# ---------------------------------------------------------------

from __future__ import print_function

import os
import sys
import re
import copy
import csv
from itertools import chain

from scipy.stats import rankdata
from numpy import array

try:
    from fugassem.common.utils import try_open
except:
    sys.exit("fugassem is not installed!")


# ---------------------------------------------------------------
# constants 
# ---------------------------------------------------------------

c_strHeaders        = "headers"
c_strNA             = "#N/A"
c_iMaxPrintChoices  = 3
c_iMaxPeekChoices   = 5

# ---------------------------------------------------------------
# beginning of the table object
# ---------------------------------------------------------------

class table:

    """
    A class for representing and manipulating tabular data
    First row is unique col headers
    First col is unique row headers
    """

    # ---------------------------------------------------------------
    # initialize the table
    # ---------------------------------------------------------------
    
    def __init__( self, source=None, skip_remap=False, verbose=True ):
        """ Constructor can "open" ( i ) list of lists ( ii ) tab delimitted file/STDIN """
        # establish the table from: ( i ) list of lists or ( ii ) file ( handle )
        if isinstance( source, list ):
            self.data = [aRow[:] for aRow in source]
            self.source = "<list of lists>"
        else:
            with try_open( source ) if isinstance( source, str ) else sys.stdin as fh:
                self.data = [row for row in csv.reader( fh, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE )]
                self.source = source if source is not None else "<stdin>"
        # track transposition status
        self.istransposed=False
        # track verbosity status (used by the reporter)
        self.isverbose=verbose
        # attempt to set up the initial index for table
        if not skip_remap:
            self.remap()
            self.report( "new table with size", self.size() )
        else:
            self.report( "WARNING: table loaded in unmapped mode" )

    # ---------------------------------------------------------------
    # set up / update the indexing system
    # ---------------------------------------------------------------

    def remap( self ):
        """ rebuilding table indexing after instantiation or modification ( e.g. transpose ) """
        # convenience variables for knowing headers
        self.colheads = self.data[0][1:]
        self.rowheads = [row[0] for row in self.data[1:]]
        # convenience variables for mapping headers to indexes
        # note, the overall header ( 0,0 ) can be called by name or this default
        self.colmap = { c_strHeaders:0 }
        self.rowmap = { c_strHeaders:0 }
        # define colmap; check for collisions
        # **** elif is special case for tables that used c_strHeaders originally ****
        for i, header in enumerate( self.data[0] ):
            if header not in self.colmap:
                self.colmap[header] = i
            elif i == 0 and header == c_strHeaders:
                pass                
            else:
                self.report( "COL COLLISION:", header, "to be replaced with", header+"-dup" )
                header += "-dup"
                self.colmap[header] = i
                """
                self.report( "EXITING: col collision:", header, "at", i, "defined at", self.colmap[header] )
                sys.exit()
                """
        # define rowmap; check for collisions
        # **** elif is special case for tables that used c_strHeaders originally ****
        for i, row in enumerate( self.data ):
            header = row[0]
            if header not in self.rowmap:
                self.rowmap[header] = i
            elif i == 0 and header == c_strHeaders:
                pass
            else:
                self.report( "ROW COLLISION:", header, "to be replaced with", header+"-dup" )
                header += "-dup"
                self.rowmap[header] = i
                """
                self.report( "EXITING: row collision:", header, "at", i, "defined at", self.rowmap[header] )
                sys.exit()
                """
        # check that all rows have the same length
        aRowLens = [len( aRow ) for aRow in self.data]
        if len( set( aRowLens ) ) != 1:
            self.report( "EXITING: Not all rows have the same length", set( aRowLens ) )
            sys.exit()

    # ---------------------------------------------------------------
    # utilities
    # ---------------------------------------------------------------

    def __getitem__( self, value ):
        """ hack to allow [r][c] indexing """
        pass
    
    def rowdex( self, index ):
        """ Convert numerical or string rowhead to numerical rowhead index """
        return index if isinstance( index, int ) else self.rowmap[index]

    def coldex( self, index ):
        """ Convert numerical or string colhead to numerical colhead index """
        return index if isinstance( index, int ) else self.colmap[index]

    def report( self, *args ):
        """ generic reporter """
        if self.isverbose:
            print( self.source, ":", " ".join( [str( k ) for k in args] ), file=sys.stderr )
        
    def peek( self, r=c_iMaxPeekChoices, c=c_iMaxPeekChoices ):
        """ show part of the table """
        for row in self.data[0:min( len( self.data ), r+1 )]:
            print( "\t".join( [str( k ) for k in row[0:min( len( row ), c+1 )]] ), file=sys.stderr )
        
    def size( self ):
        """ returns size of table """
        return "<%s%d ROW x %d COL>" % ( "!transposed! " if self.istransposed else "", len( self.rowheads ), len( self.colheads ) )

    def copy( self ):
        """ Returns a DEEP copy of the table """
        return copy.deepcopy( self )

    def dump( self, output_file=None ):
        """ Print the table to a file """
        fh = try_open( output_file, "w" ) if output_file is not None else sys.stdout
        dumper = csv.writer( fh, delimiter="\t", quotechar="", quoting=csv.QUOTE_NONE )
        for row in self.data:
            dumper.writerow( row )
        fh.close()

    def rowsort( self, order=None ):
        """ alphasorts the rows based on rowheads """
        order = sorted( self.rowheads ) if order is None else order
        self.data = [self.data[0]] + [self.data[self.rowmap[rowhead]] for rowhead in order]
        self.remap()

    def colsort( self, order=None ):
        """ alphasorts the cols based on colheads ( via transpose rowsort ) """
        # no need for final remap, since it's included in transpose
        self.transpose()
        self.rowsort( order=order)
        self.transpose()

    # ---------------------------------------------------------------
    # slicing
    # ---------------------------------------------------------------

    def row( self, index, start=1 ):
        """ Returns requested row as a list; start=0 includes rowhead """
        index = self.rowdex( index )
        return self.data[index][start:]

    def col( self, index, start=1 ):
        """ Returns requested col as a list; start=0 includes colhead """
        index = self.coldex( index )
        return [row[index] for row in self.data[start:]]

    def rowdict( self, index ):
        """ Returns requested row as a [colhead]=value dictionary """
        temp_row = self.row( index )
        return {self.colheads[k]:temp_row[k] for k in range( len( self.colheads ) )}

    def coldict( self, index ):
        """ Returns requested col as a [rowhead]=value dictionary """
        temp_col = self.col( index )
        return {self.rowheads[k]:temp_col[k] for k in range( len( self.rowheads ) )}

    def entry( self, r, c ):
        """ Returns the ( r, c ) entry of the table; r and c can be int or named index """
        r, c = self.rowdex( r ), self.coldex( c )
        return self.data[r][c]

    # ---------------------------------------------------------------
    # appliers
    # ---------------------------------------------------------------

    def transpose( self ):
        """ transpose the table """
        self.data = [[row[index] for row in self.data] for index in range( len( self.data[0] ) )]
        self.istransposed = not self.istransposed
        self.remap()
        self.report( "transposed the table" )

    def promote( self, index ):
        """ Pick a row to become the new row[0]; i.e. make it the colhead row """
        self.data[0] = self.data.pop( self.rowdex( index ) )
        self.remap()

    def set( self, r, c, value ):
        """ Set the ( r, c ) entry of the table """
        r, c = self.rowdex( r ), self.coldex( c )
        self.data[r][c] = value

    def apply_rowheads( self, function ):
        """ applies a function to the rowheads and remaps ( to avoid collisions ) """
        for rowhead in self.rowheads:
            self.data[self.rowdex( rowhead )][0] = function( rowhead )
        self.remap()
                             
    def apply_colheads( self, function ):
        """ applies a function to the colheads and remaps ( to avoid collisions ) """
        for colhead in self.colheads:
            self.data[0][self.coldex( colhead )] = function( colhead )
        self.remap()

    def apply_entries( self, function ):
        """ applies a function to each entry in the table """
        # updated to be more efficient; now based on list comprehension not iteration
        for i in range( 1, len( self.data ) ):
            self.data[i] = [self.data[i][0]] + [function( k ) for k in self.data[i][1:]]

    def apply_row( self, index, function ):
        """ applies a function to each entry in a row """
        rowhead = self.rowdex( index )
        for colhead, value in self.rowdict( rowhead ).items( ):
            self.set( rowhead, colhead, function( value ) )

    def apply_col( self, index, function ):
        """ applies a function to each entry in a col """
        colhead = self.coldex( index )
        for rowhead, value in self.coldict( colhead ).items( ):
            self.set( rowhead, colhead, function( value ) )
                 
    # ---------------------------------------------------------------
    # generators
    # ---------------------------------------------------------------
    
    def iter_rows( self ):
        """ iterate over rows; yields rowhead, row """
        for rowhead in self.rowheads:
            yield rowhead, self.row( rowhead )

    def iter_cols( self ):
        """ iterate over cols; yields colhead, col """                         
        for colhead in self.colheads:
            yield colhead, self.col( colhead )

    def iter_entries( self ):
        """ iterate over entries; yields ( r )owhead, ( c )olhead, entry( r, c ) """
        for i, rowhead in enumerate( self.rowheads ):
            for j, colhead in enumerate( self.colheads ):
                # **** (row|col)heads don't include the origin position, so must offset i,j by 1 ****
                yield rowhead, colhead, self.data[i+1][j+1]

    # ---------------------------------------------------------------
    # reduce method, operates like python's filter on rows
    # ---------------------------------------------------------------

    def reduce( self, function, protect_headers=True, transposed=False, invert=False, in_place=True ):
        """ apply a function to the rows of the table and rebuild with or return true evals """
        # flip?
        if transposed: 
            self.transpose()
        # auto-pass headers?
        data2 = [self.data[0] if protect_headers else []]
        start = 1 if protect_headers else 0
        # apply test to rows
        if not invert:
            data2 += [row for row in self.data[start:] if function( row )]
        else:
            data2 += [row for row in self.data[start:] if not function( row )]
        # (1) rebuild table
        if in_place:
            self.data = data2
            if transposed: 
                self.transpose()
            self.remap()
            self.report( "--> reduced size is", self.size() )
            return None
        # (2) return rows that passed as new table
        else:
            if transposed: 
                self.transpose()
            # rebuild data2 to avoid referencing self.data rows
            data2 = [row[:] for row in data2]
            # this returns to the outer function, which must also return
            new_table = table( data2, verbose=self.isverbose )
            # if we were working on a transposed table, must also transpose selected rows
            if transposed: 
                new_table.transpose()
            return new_table

    # ---------------------------------------------------------------
    # methods that call reduce to do things
    # ---------------------------------------------------------------

    def grep( self, index, patterns, **kwargs ):
        """ restrict rows to those whose col[index] position matches a pattern """
        if isinstance( patterns, str ):
            patterns = [patterns]
        self.report( "applying grep", "index=", index, "patterns=", pretty_list( patterns ), kwargs )
        return self.reduce( lambda row: any( [re.search( k, row[self.coldex( index )] ) for k in patterns] ), **kwargs )

    def select( self, index, choices, **kwargs ):
        """ select rows whose col[index] entry is in choices """
        # this allows us to test if choices is a scalar or list; strange because strings are iterable?
        if isinstance( choices, str ):
            choices = [choices]
        self.report( "applying select", "index=", index, "choices=", pretty_list( choices ), kwargs )
        return self.reduce( lambda row: row[self.coldex( index )] in choices, **kwargs )

    def delete( self, index, choices, **kwargs ):
        """ delete rows whose col[index] entry is in choices """
        if isinstance( choices, str ):
            choices = [choices]
        # this allows us to test if choices is a scalar or list; strange because strings are iterable?
        self.report( "applying delete", "index=", index, "choices=", pretty_list( choices ), kwargs )
        return self.reduce( lambda row: row[self.coldex( index )] not in choices, **kwargs )
    
    def delete_row( self, index, **kwargs ):
        """ delete row with specific index """
        self.report( "applying raw", "index=", index, kwargs )
        return self.reduce( lambda row: self.rowdex( row[0] ) != self.rowdex( index ), **kwargs )

    def head( self, index, **kwargs ):
        """ keep only rows up to head; if inverted, delete rows up to head ( useful for stripping metadata ) """
        self.report( "applying head", "index=", index, kwargs )
        return self.reduce( lambda row: self.rowdex( row[0] ) <= self.rowdex( index ), **kwargs )

    def limit( self, index, operation, **kwargs ):
        """ """
        op, threshold = re.search( "([<>=]+)(.*)", operation ).groups()
        threshold = float( threshold )
        choices = { 
            "<" : lambda x: float( x ) <  threshold,
            "<=": lambda x: float( x ) <= threshold,
            ">" : lambda x: float( x ) >  threshold,
            ">=": lambda x: float( x ) >= threshold,
            }
        myfunc = choices[op]
        self.report( "applying limit, requiring field", index, "to be", op, threshold, kwargs )
        return self.reduce( lambda row: myfunc( row[self.coldex( index )] ), **kwargs )

    def unrarify( self, minlevel=0, mincount=1, **kwargs ):
        """ applies to numerical table: keep features ( row ) that exceed specified level in specified # of samples """
        self.report( "applying unrarify", "requiring at least", mincount, "row values exceeding", minlevel, kwargs )
        return self.reduce( lambda row: mincount <= sum( [1 if item > minlevel else 0 for item in row[1:]] ), **kwargs )

    # ---------------------------------------------------------------
    # groupby
    # ---------------------------------------------------------------

    def groupby( self, funcGrouper, funcSummarizer ):
        """ user grouper function on rowheads to cluster rows, then use summarizer function to combine values """ 
        groups = {} # map new rowheads to 1+ old rowheads
        for rowhead in self.rowheads:
            groupdict = groups.setdefault( funcGrouper( rowhead ), {} )
            groupdict[rowhead] = 1
        data2 = [self.data[0]] # colheads
        for group, dictRowheads in groups.items():
            aaRows = [self.row( rowhead ) for rowhead in dictRowheads]
            aaRowsTransposed = [[aRow[i] for aRow in aaRows] for i in range( len( aaRows[0] ) )]
            data2.append( [group] + [funcSummarizer( aRow ) for aRow in aaRowsTransposed] )
        self.data = data2
        self.remap()
        self.report( "applied groupby:", "rowheads now like <%s>" % ( self.rowheads[0] ), "; new size is", self.size() )

    # ---------------------------------------------------------------
    # stratification methods
    # ---------------------------------------------------------------

    def stratify( self, focus ):
        """ stratifies table on focal row value; returns dictionary of tables """
        self.transpose()
        dictDataArrays = {}
        for rowhead, row in self.iter_rows( ):
            stratum = self.entry( rowhead, focus )
            # the default is a copy of this table's rowheads (currently colheads)
            dictDataArrays.setdefault( stratum, [self.data[0][:]] ).append( [rowhead] + row )
        dictTables = {}
        for stratum, aaData in dictDataArrays.items():
            dictTables[stratum] = table( aaData, verbose=self.isverbose )
            dictTables[stratum].transpose( )
        self.transpose()
        self.report( "stratified on", focus, "into", len( dictTables ), "new tables" )
        return dictTables

    def funcify( self, func ):
        """ stratifies table on function return value """
        self.transpose( )
        dictDataArrays = {}
        for rowhead, row in self.iter_rows( ):
            ret = func( rowhead )
            if ret is not None:
                # the default is a copy of this table's rowheads (currently colheads)
                dictDataArrays.setdefault( ret, [self.data[0][:]] ).append( [rowhead] + row )
        dictTables = {}
        for stratum, aaData in dictDataArrays.items( ):
            dictTables[stratum] = table( aaData, verbose=self.isverbose )
            dictTables[stratum].transpose( )
        self.transpose( )
        self.report( "funcified into", len( dictTables ), "new tables" )
        return dictTables

    def groupify( self, funcGrouper ):
        """ stratifies table by return value of function on rowheads; returns dictionary of tables """
        dictDataArrays = {}
        for rowhead, row in self.iter_rows():
            group = funcGrouper( rowhead )
            # the default is a copy of this table's colheads
            dictDataArrays.setdefault( group, [self.data[0][:]] ).append( [rowhead] + row )
        dictTables = {}
        for group, aaData in dictDataArrays.items():
            dictTables[group] = table( aaData, verbose=self.isverbose )
        self.report( "groupified on f( rowhead ) into", len( dictTables ), "new tables" )
        return dictTables

    # ---------------------------------------------------------------
    # conversion methods
    # ---------------------------------------------------------------

    def table2tupledict( self ):
        """ return table as dictionary with tuple keys ~ ( rowheads[i], colheads[j] ) = data[i][j] """
        return {( rowhead, colhead ):value for rowhead, colhead, value in self.iter_entries()}

    def table2nesteddict( self ):
        """ return table as nested dictionary ~ [rowheads[i]][colheads[j]] = data[i][j] """
        return {rowhead:self.rowdict( rowhead ) for rowhead in self.rowheads}

    def table2array( self, last_metadata=None ):
        """ return just quant data as numpy 2d array """
        temp = []
        for rowhead, row in self.iter_rows( ):
            if last_metadata is None or self.rowmap[rowhead] > self.rowmap[last_metadata]:
                temp.append( map( float, row ) )
        return array( [array( row ) for row in temp] )
    
    # ---------------------------------------------------------------
    # methods for extending/combining tables
    # ---------------------------------------------------------------

    def insert( self, index, row ):
        """ inserts a pre-formatted row ( has header; proper order ) into the table before index """
        index = self.rowdex( index )
        self.data.insert( index, row )
        self.remap()

    def augment( self, table2 ):
        """ join table with table2 on colheads """
        setRowheadOverlap = set( self.rowheads ).__and__( set( table2.rowheads ) )
        if len( setRowheadOverlap ) > 0:
            self.report( "augmentee contains", len( setRowheadOverlap ), "duplicate rowheads (skip)" )
        # note: "colhead in dictColmap" much faster than "colhead in aColheads"
        self.data = self.data + [
            [rowhead2] + [table2.entry( rowhead2, colhead ) if colhead in table2.colmap else c_strNA for colhead in self.colheads] 
            for rowhead2 in table2.rowheads if rowhead2 not in self.rowmap
            ]
        self.remap()
        self.report( "augmented with rows from", table2.source, "new size is", self.size() )

    def extend( self, table2 ):
        """ join table with table2 on rowheads """
        setColheadOverlap = set( self.colheads ).__and__( set( table2.colheads ) )
        if len( setColheadOverlap ) > 0:
            self.report( "extendee contains", len( setColheadOverlap ), "duplicate colheads (skip)" )
        # find index for unique cols in table2
        aIndex = [j for j, colhead in enumerate( table2.colheads ) if colhead not in self.colmap]
        for i in range( len( self.data ) ):
            rowhead = self.data[i][0]
            # headers are a special case (may not align on rowhead if 0,0 entries differ)
            if i == 0:
                self.data[0] += [table2.colheads[j] for j in aIndex]
            # a shared row, add new col entries
            elif rowhead in table2.rowmap:
                row = table2.row( rowhead )
                self.data[i] += [row[j] for j in aIndex]
            # a self-specific row, add all NA entries
            else:
                self.data[i] += [c_strNA for j in aIndex]
        self.remap()
        self.report( "extended with cols from", table2.source, "new size is", self.size() )

    def merge( self, table2 ):
        """ slightly faster merge conducted as an extend and augment """
        self.extend( table2 )
        self.augment( table2 )

    # ---------------------------------------------------------------
    # methods for filling "empty" cells
    # ---------------------------------------------------------------

    def blank2na( self, value=c_strNA ):
        """ Anything that doesn't match a non-space char is filled with na """
        self.apply_entries( lambda k: k if re.search( "[^\s]", str( k ) ) else value )
        
    def na2zero( self, value=0 ):
        """ Convert NAs to zero values """
        self.apply_entries( lambda k: value if k == c_strNA else k )

    # ---------------------------------------------------------------
    # methods for working with metadata
    # ---------------------------------------------------------------

    def metasplit( self, last_metadata_row ):
        """ splits a table into top ( metadata ) and bottom ( features ); returns """
        tableMetadata = self.head( last_metadata_row, in_place=False )
        tableFeatures = self.head( last_metadata_row, in_place=False, invert=True )
        return tableMetadata, tableFeatures

    def metamerge( self, tableMetadata ):
        """ merges metadata ON TOP of current table """
        dataTemp = self.data[1:] # data rows
        self.head( 0 )
        self.augment( tableMetadata )
        self.data += dataTemp
        self.remap()
        self.report( "added metadata from", tableMetadata.source, "new size is", self.size() )

    # ---------------------------------------------------------------
    # methods based on math
    # ---------------------------------------------------------------

    def compress_zeroes( self ):
        """ replaces 0.0 with 0 for compression """
        self.apply_entries( lambda k: 0 if k == 0.0 else k )

    def float( self ):
        """ Attempt to float all non-header entries in the table """
        self.apply_entries( float )

    def unfloat( self ):
        """ convert 0.0 to 0 and reduce others to N sig figs for compression """
        self.apply_entries( lambda k: "0" if k == 0.0 else "%.6g" % ( k ) )

    def normalize_columns( self ):
        """ Normalizes the columns. Fails if there are non-numeric entries. """
        # transpose because working on rows is easier
        self.transpose()
        # iter over rows
        for i, header in enumerate( self.rowheads ):
            # note, rowhead 0 is row 1 of main table ***
            row = self.row( header )
            total = float( sum( row ) )
            if total > 0:
                row = [ k/total for k in row ]
            else:
                row = [ 0.0 for k in row ]
                self.report( "WARNING:", "sum of column", header, self.rowdex( header ), "is 0" )
            self.data[i+1] = [header] + row # ***
        # flip back
        self.transpose()
        # report
        self.report( "normalized the columns" )

    def rank_columns( self, normalize=False ):
        """ Ranks the columns. Fails if there are non-numeric entries. """
        # transpose because working on rows is easier
        self.float( )
        self.transpose()
        # iter over rows
        for i, header in enumerate( self.rowheads ):
            # note, rowhead 0 is row 1 of main table ***
            row = self.row( header )
            row = list( rankdata( row ) ) # default return is numpy array
            row = row if not normalize else [k/max( row ) for k in row]
            self.data[i+1] = [header] + row # ***
        # flip back
        self.transpose()

    # ---------------------------------------------------------------
    # end of the table object
    # ---------------------------------------------------------------

# ---------------------------------------------------------------
# table-related methods outside of the table object
# ---------------------------------------------------------------

def pretty_list( items ):
    """ simplifies printing long argument lists """
    items = list( items )
    if len( items ) <= c_iMaxPrintChoices:
        return items
    else:
        return items[0:c_iMaxPrintChoices] + ["{{ + %d others }}" % ( len( items ) - c_iMaxPrintChoices )]

def nesteddict2table( d, aRowheads=None, aColheads=None, origin=c_strHeaders, empty=c_strNA ):
    """ convert a nested dictionary [A][B]=value to a table with A as rows and B as columns;
    default ordering for A and B is alphabetical; can impose alternate order with aRowheads
    and aColheads; when encounting a missing cell, fill it with the empty value """
    if aRowheads is None:
        aRowheads = sorted( d.keys() )
    if aColheads is None:
        # for dictionary [A][B], get sorted list of unique Bs
        aColheads = sorted( set( chain.from_iterable( [d[k].keys() for k in d] ) ) )
    aaData = [ [origin] + aColheads ]
    for rowhead in aRowheads:
        inner = d.get( rowhead, {} )
        aRow = [ rowhead ] + [ inner.get( colhead, empty ) for colhead in aColheads ]
        aaData.append( aRow )
    return table( aaData, verbose=False )

def tupledict2table( d, aRowheads=None, aColheads=None, empty=c_strNA ):
    """ turn a dictionary [( r, c )] into a table with r as rows and c as cols;
    coerces first to [a][b]-style nested dictionary and then called nesteddict2table """
    dictdictTemp = {}
    for ( r, c ), value in d.items():
        if r not in dictdictTemp:
            dictdictTemp[r] = {}
        dictdictTemp[r][c] = value
    return nesteddict2table( dictdictTemp, aRowheads=aRowheads, aColheads=aColheads, empty=empty )
