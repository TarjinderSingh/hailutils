#!/usr/bin/python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *

def extract_region(vds, locus, start = None, end = None):
    interval_str = ''
    if locus and not start and not end:
        interval_str = '{}'.format(str(locus))
    if locus and start and end:
        interval_str = '{}:{}-{}'.format(locus, str(start), str(end))
    if locus and start and not end:
        interval_str = '{}:{}-{}'.format(locus, str(start), str(start + 1))
    return(vds.filter_intervals(Interval.parse(interval_str)))
    
def extract_region_file(vds, region_file, keep = True, 
                        bed = False, interval_list = False, locus = False, variant = False,
                        header = None):
    ext = re.sub('\.bgz$', '', os.path.splitext(region_file)[1])
    if bed or ext == '.bed':
        kt = KeyTable.import_bed(region_file)
    elif interval_list or ext == '.interval_list':
        kt = KeyTable.import_interval_list(region_file)
    elif locus or variant:
        if header:
            kt = hc.import_table(region_file, impute = True)
        else:
            kt = hc.import_table(region_file, no_header = True, impute = True)
            if locus:
                kt = kt.rename(['Locus'])
            else:
                kt = kt.rename(['Variant'])
        kt = kt.key_by(kt.columns[0])
    return(vds.filter_variants_table(kt, keep = keep))
    
def extract_variants(vds, variant_str_list):
    variant_str_list = [ variant_str_list ] if isinstance(variant_str_list, str) else variant_str_list
    return(vds.filter_variants_list([ Variant.parse(vstr) for vstr in variant_str_list ], keep = True))

def extract_genotypes(vds):
    return

def sites(vds):
    return(vds.drop_samples())

def subset_for_testing(vds, region = '22:15M-20M'):
    return(extract_region(vds, region))
    
def extract_inds(vds, sample_ids, keep = True):
    vds = vds.annotate_global('global.sample_ids', sample_ids, TSet(TString()))
    if keep:
        vds = vds.filter_samples_expr('global.sample_ids.contains(s)')
    else:
        vds = vds.filter_samples_expr('! global.sample_ids.contains(s)')
    return(vds)

def change_contig(outfile, infile, add = True):
    cmd = 'zless {} | modify_vcf_contig.py'.format(infile)
    cmd += ' -a' if add else ' -r'
    cmd += ' | bgzip -c > {}'.format(outfile)
    subprocess.call(cmd, shell = True)