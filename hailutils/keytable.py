#!/usr/bin/python

import os
import sys
import re
import logging
import itertools
from hail import *
from pprint import pprint

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def join_lists(*args):
    return(list(itertools.chain(*args)))

def flatten_list_of_lists(l):
    return(list(itertools.chain.from_iterable(l)))

def change_type(str_to_set = None, set_to_str = None, hail_to_str = None, sep = ','):
    annotations = [] 
    if isinstance(str_to_set, list) and str_to_set:
        for col in str_to_set:
            annotations.append('{0} = {0}.split("{1}").toSet()'.format(col, sep))
    if isinstance(set_to_str, list) and set_to_str:
        for col in set_to_str:
            annotations.append('{0} = {0}.mkString("{1}")'.format(col, sep)) 
    if isinstance(hail_to_str, list) and hail_to_str:
        for col in hail_to_str:
            annotations.append('{0} = str({0})'.format(col))             
    return(annotations)

def set_type_missing(set_null = None, int_zero = None):
    annotations = []
    if isinstance(set_null, list) and set_null:
        for col in set_null:
            annotations.append('{0} = if (isDefined({0})) {0} else [""][:0].toSet()'.format(col))
    if isinstance(int_zero, list) and int_zero:
        for col in int_zero:
            annotations.append('{0} = if (isDefined({0})) {0} else 0'.format(col))
    return(annotations)

def index_before_join(kt, index, cols):
    return(kt.rename({ col:'{}{}'.format(col, index) for col in cols }))

def interval_as_variant(interval_col = 'region', variant_col = 'variant', ref = 'G', alt = 'C'):
    return(['{0} = Variant({1}.start.contig, {1}.start.position, "{2}", "{3}")'.format(variant_col, interval_col, ref, alt)])

def interval_string(interval_col = 'region', interval_str = 'region_str', replace = False):
    if replace:
        interval_str = interval_col
    return(['{0} = {1}.start.contig + ":" + str({1}.start.position) + "-" + str({1}.end.position)'.format(interval_str, interval_col)])

def interval_length(interval_col = 'region', length_col = 'length'):
    return(['{0} = {1}.end.position - {1}.start.position'.format(length_col, interval_col)])

def interval_to_bed(interval_col = 'region', chrom_col = 'chrom', start_col = 'start', end_col = 'end'):
    return([ '{1} = {0}.start.contig, {2} = {0}.start.position - 1, {3} = {0}.end.position - 1'.format(interval_col, chrom_col, start_col, end_col)])

def prettify_columns(columns, strip_match = None, strip_all = False):
    if strip_all:
        return({ col:re.split('\.', col)[-1] for col in columns })
    if strip_match:
        if not isinstance(strip_match, list):
            strip_match = [ strip_match ]
        regex = '|'.join(
            [ 
                '^(' + re.sub(r"\.", "\.", x) + '\.)' 
                 for x in 
                    sorted(
                        strip_match, 
                        reverse = True
                    ) 
            ]
        )
        return({ col:re.sub(regex, '', col) for col in columns })
    
def match_columns(columns, match = None, full_match = None):
    if full_match:
        if not isinstance(full_match, list):
            full_match = [ full_match ]
        full_match = [ x + '$' for x in full_match ]
    if match:
        if not isinstance(match, list):
            match = [ match ]
    if match and full_match:
        match = match + full_match
    elif not match and full_match:
        match = full_match
    regex = '|'.join(
        [ 
            '(' + re.sub(r"\.", "\.", x) + ')' 
             for x in 
                sorted(match, reverse = True) 
            
        ]
    )
    return([ col for col in columns if re.match(regex, col) ])   

def unflatten(kt, starts_with = 'va'):
    '''
    Unflattens KeyTable.
    '''
    exprs = []
    for col in kt.columns:
        if col.startswith(starts_with):
            exprs.append('{0} = `{0}`'.format(col))
    drop_cols = [ '{}'.format(col) for col in kt.columns if col.startswith(starts_with) ]
    return(kt.annotate(exprs).drop(drop_cols))

def flattened_to_vds(kt):
    '''
    Converts flattened KeyTable to VDS and unflattens structure
    '''
    exprs = []
    for col in kt.columns:
        if col.startswith('va'):
            exprs.append('{0} = va.`{0}`'.format(col))
    drop_cols = ', '.join([ '`{}`'.format(f) for f in kt.columns if f.startswith('va') ])
    return(
        VariantDataset
            .from_table(kt)
            .annotate_variants_expr(exprs)
            .annotate_variants_expr('va = drop(va, {})'.format(drop_cols))   
    )

def to_vds(kt):
    return(VariantDataset.from_table(kt))