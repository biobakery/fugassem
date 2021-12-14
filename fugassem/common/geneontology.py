#!/usr/bin/env python

"""
Module/script for interacting with Gene Ontology
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import sys
import re
import argparse
from collections import Counter
from fugassem.common.utils import path2name, warn, qw
from fugassem.common.dictation import polymap, col2dict

sys.setrecursionlimit(50000)

# ---------------------------------------------------------------
# constants
# ---------------------------------------------------------------

ALLOW_CROSSTALK = False
c_prop_pattern = r"^([^\s]*?): (.*)"
c_goid_pattern = r"(GO:[0-9]+)"

c_namespace_convert = qw( """
biological_process BP
molecular_function MF
cellular_component CC
""", as_dict=True )

c_relationship_subtypes = qw( """
part_of
#regulates
#positively_regulates
#negatively_regulates
""" )

c_singular_props = qw( """
id
name
namespace
""" )

# ---------------------------------------------------------------
# global tracking variables
# ---------------------------------------------------------------

parentage_types = Counter( )

# ---------------------------------------------------------------
# begin Stanza class
# ---------------------------------------------------------------

"""
Example of a Term Stanza from an OBO file
=========================================
[Term]
id: GO:0000070
name: mitotic sister chromatid segregation
namespace: biological_process
alt_id: GO:0016359
def: "The cell cycle process in which ... the mitotic cell." [GOC:ai, GOC:jl]
subset: goslim_pombe
synonym: "mitotic chromosome segregation" EXACT []
synonym: "mitotic sister-chromatid adhesion release" NARROW []
is_a: GO:0000819 ! sister chromatid segregation
is_a: GO:1903047 ! mitotic cell cycle process
relationship: part_of GO:0000278 ! mitotic cell cycle
relationship: part_of GO:0007067 ! mitotic nuclear division
"""

class Stanza:

    """query items from a stanza of the obo file"""
      
    def __init__( self ):
        self.lines = []
        self.props = {}

    def __getitem__( self, key ):
        if key in c_singular_props:
            assert len( self.props[key] ) == 1, \
                "not a singular relationship" + key + str( self.props[key] )
            return self.props[key][0]
        else:
            return self.props[key]

    def __contains__( self, key ):
        return key in self.props

    def get( self, key, default=None ):
        return self[key] if key in self else default
    
    def populate( self ):
        propname = None
        for line in self.lines:
            match = re.search( c_prop_pattern, line )
            if match:
                propname, line = match.groups( )
                self.props.setdefault( propname, [] ).append( line )

    def get_parent_goids( self ):
        parent_goids = []
        # extract is_a relationships
        # format: list of "goid ! definition"
        for line in self.props.get( "is_a", [] ):
            parent_goids.append( line.split( )[0] )
            parentage_types["is_a"] += 1
        # extract "relationship" relationships
        # format: list of "subtype goid ! definition"
        for line in self.props.get( "relationship", [] ):
            items = line.split( )
            if items[0] in c_relationship_subtypes:
                parent_goids.append( items[1] )
                parentage_types["relationship:"+items[0]] += 1
        # extract "intersection_of" relationships
        # format: list of "goid ! definition" || "subtype goid ! definition"
        for line in self.props.get( "intersection_of", [] ):
            items = line.split( )
            for k in [0, 1]:
                if "GO:" in items[k]:
                    parent_goids.append( items[k] )
                    parentage_types["intersection_of"] += 1
        # confirm parents look like GO ids
        return [k for k in parent_goids if re.search( c_goid_pattern, k )]

class OBOParser:

    def __init__( self, p_obo ):
        IN_TERMS = False
        self.stanzas = []
        with open( p_obo ) as fh:
            for line in fh:
                line = line.strip( )
                if line == "[Term]":
                    self.stanzas.append( Stanza( ) )
                    IN_TERMS = True
                elif IN_TERMS and line == "[Typedef]":
                    # signals end of the term stanzas
                    IN_TERMS = False
                elif IN_TERMS and line != "":
                    self.stanzas[-1].lines.append( line )
        # properties of a stanza are "populated" once all lines added
        [s.populate( ) for s in self.stanzas]

    def extract_terms( self ):
        for stanza in self.stanzas:
            yield Term( stanza )
                
# ---------------------------------------------------------------
# begin Term class
# ---------------------------------------------------------------

class Term:

    """Represents a Term in GO; Node in GO DAG"""

    def __init__( self, stanza ):        
        # set items as derived from the stanza
        self.goid            = stanza["id"]
        self.name            = stanza["name"]
        self.namespace       = stanza["namespace"]
        self.namespace_short = c_namespace_convert[self.namespace]
        self.is_obsolete     = True if "is_obsolete" in stanza else False
        self.replaced_by     = stanza.get( "replaced_by", None )
        self.alt_ids         = stanza.get( "alt_id", [] )
        self.parent_ids      = stanza.get_parent_goids( )
        # topological properties set by the ontology
        self.parents         = set( )
        self.children        = set( )
        self.genes           = set( )
        self.depth           = None
        self.is_root         = False
        self.is_leaf         = False
        # set by optional methods from the ontology
        self.progeny         = None
        self.progeny_genes   = None
        self.is_informative  = False
        self.is_pruned       = False
        self.is_acceptable   = True

    def __repr__( self ):
        return "{}: [{}] {}".format( self.goid, self.namespace_short, self.name )

    def add_gene( self, gene ):
        self.genes.add( gene )

    def get_progeny( self ):
        """ memoized via instance variable """
        if self.progeny is None:
            self.progeny = {self}
            for cterm in self.children:
                self.progeny.update( cterm.get_progeny( ) )
        return self.progeny

    def get_progeny_genes( self ):
        """ memoized via instance variable """
        if self.progeny_genes is None:
            self.progeny_genes = set( )
            self.progeny_genes.update( self.genes )
            for cterm in self.children:
                self.progeny_genes.update( cterm.get_progeny_genes( ) )
        return self.progeny_genes

# ---------------------------------------------------------------
# begin Ontology class
# ---------------------------------------------------------------

class Ontology:

    """Representation of the GO Ontology"""

    def __init__( self, p_obo ):

        # mapping from goid to term object
        self.terms = {}
        self.idmap = {}
        self.roots = []
        self.leaves = []
        self.attached_genes = set( )

        # populate terms and lookup
        for term in OBOParser( p_obo ).extract_terms( ):
            if not re.search( c_goid_pattern, term.goid ):
                continue
            elif term.is_obsolete:
                if term.replaced_by is not None:
                    self.idmap[term.goid] = term.replaced_by
            else:
                self.terms[term.goid] = term
                self.idmap[term.goid] = term.goid
                for alt_id in term.alt_ids:
                    self.idmap[alt_id] = term.goid

        # populate parent/child relationships
        for cterm in self.terms.values( ):
            for parent_id in cterm.parent_ids:
                pterm = self.terms[parent_id]
                if ALLOW_CROSSTALK or ( pterm.namespace == cterm.namespace ):
                    cterm.parents.add( pterm )
                    pterm.children.add( cterm )
                else:
                    warn( "parent->child spans namespaces (ignoring):", "\n",
                          "\t", pterm.__repr__( ), "\n",
                          "\t", cterm.__repr__( ), "\n", )
                    parentage_types["Cross-ontology [ignored]"] += 1

        # identify roots, leaves
        for termid, term in self.terms.items( ):
            if len( term.parents ) == 0:
                term.is_root = True
                self.roots.append( term )
            if len( term.children ) == 0:
                term.is_leaf = True
                self.leaves.append( term )

        # add depth information
        def recurse_depth( term ):
            step_down = term.depth + 1
            for cterm in term.children:
                if ( cterm.depth is None ) or \
                   ( cterm.depth > step_down ):
                    cterm.depth = step_down
                    recurse_depth( cterm )
        for root in self.roots:
            root.depth = 0
            recurse_depth( root )

    def iter_terms( self ):
        for goid, term in self.terms.items():
            yield term

    def attach_genes( self, polymap ):
        for gene, goids in polymap.items( ):
            self.attached_genes.add( gene )
            for goid in goids:
                # handle updated/aliased ids
                goid = self.idmap.get( goid, goid )
                if goid in self.terms:
                    self.terms[goid].add_gene( gene )

    def prune( self, goid ):
        def recurse_prune( term ):
            term.is_pruned = True
            for cterm in term.children:
                if not cterm.is_pruned:
                    recurse_prune( cterm )
        recurse_prune( self.terms[goid] )

    def set_informative( self, threshold ):
        def recurse_set_informative( term ):
            if len( term.get_progeny_genes( ) ) >= threshold:
                term.is_informative = True
                # check that no children will be called informative
                for cterm in term.children:
                    if len( cterm.get_progeny_genes( ) ) >= threshold:
                        term.is_informative = False
                        break
            else:
                for pterm in term.parents:
                    recurse_set_informative( pterm )
        for leaf in self.leaves:
            recurse_set_informative( leaf )

# ---------------------------------------------------------------
# utility functions
# ---------------------------------------------------------------

def get_args( ):
    """ """
    parser = argparse.ArgumentParser( )
    parser.add_argument( "obo", 
                         help="GO-provided obo file" )
    parser.add_argument( "--mapping",
                         default=None,
                         metavar="<path>",
                         help="mapping from genes to GO terms" )
    parser.add_argument( "--allowed-genes",
                         default=None,
                         metavar="<path>",
                         help="subset of allowed genes" )
    parser.add_argument( "--flip",
                         action="store_true",
                         help="mapping is from GO terms to genes" )
    parser.add_argument( "--depth",
                         type=int,
                         default=None,
                         metavar="<int>",
                         help="trimming option: only terms at this min distance from root" )
    parser.add_argument( "--grep",
                         default=None,
                         metavar="<regex>",
                         help="trimming option: only terms matching this regex" )
    parser.add_argument( "--prune",
                         default=None,
                         metavar="<term>",
                         help="trimming option: only terms descending from this term (GO id)" )
    parser.add_argument( "--namespace",
                         default=None,
                         choices=["BP", "MF", "CC"],
                         nargs="+", 
                         metavar="<BP/MF/CC>",
                         help="trimming option: only terms in this namespace" )
    parser.add_argument( "--informative",
                         default=None, 
                         metavar="<number OR fraction of genes>",
                         help="trimming option: only terms that are informative at a given level" )
    parser.add_argument( "--ignore-progeny",
                         action="store_true", 
                         help="do not let more specific annotations rise through the dag" )
    parser.add_argument( "--terms-only",
                         action="store_true", 
                         help="just output the terms list (no annotations)" )
    parser.add_argument( "--outfile",
                         default=None, 
                         help="output file" )
    args = parser.parse_args( )
    # warnings
    if args.ignore_progeny:
        warn( "Only considering direct annotations. Annotations will not rise through the DAG." )
    return args

# ---------------------------------------------------------------
# main
# ---------------------------------------------------------------

def main( ):

    args = get_args( )

    # load obo / report rel type
    obo = Ontology( args.obo )
    warn( "Summary of relationship types:" )
    for k in sorted( parentage_types ):
        warn( k, parentage_types[k] )  

    # attach genes
    if args.mapping is not None:
        mapping = polymap( args.mapping, reverse=args.flip )
        if args.allowed_genes is not None:
            allowed = col2dict( args.allowed_genes )
            mapping = {k:v for k, v in mapping.items( ) if k in allowed}
        obo.attach_genes( mapping )
        warn( "# of attached genes:", len( obo.attached_genes ) )

    # informative cut
    if args.informative is not None:
        threshold = float( args.informative )
        if threshold < 1:
            warn( "Intepretting informative cutoff as fraction of annotated genes" )
            threshold *= len( obo.attached_genes )
        threshold = int( threshold )
        obo.set_informative( threshold )
        for term in obo.iter_terms( ):
            if not term.is_informative:
                term.is_acceptable = False      

    # pruning cut
    if args.prune is not None:
        obo.prune( args.prune )
        for term in obo.iter_terms( ):
            if not term.is_pruned:
                term.is_acceptable = False

    # depth cut
    if args.depth is not None:
        for term in obo.iter_terms( ):
            if term.depth != args.depth:
                term.is_acceptable = False

    # grep cut
    if args.grep is not None:
        for term in obo.iter_terms( ):
            if not re.search( args.grep, term.name ):
                term.is_acceptable = False

    # namespace cut
    if args.namespace is not None:
        for term in obo.iter_terms( ):
            if term.namespace_short not in args.namespace:
                term.is_acceptable = False

    # output the new polymap
    fh = open( args.outfile, "w" ) if args.outfile is not None else sys.stdout
    for term in obo.iter_terms( ):
        if term.is_acceptable:
            outline = [str( term )]
            if not args.terms_only:
                outline += list( term.get_progeny_genes( ) if not args.ignore_progeny else term.genes )
            fh.write("\t".join( outline ) + "\n")
    fh.close( )

if __name__ == "__main__":
    main( )
