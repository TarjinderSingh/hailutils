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
    
def sites(vds):
    return(vds.drop_samples())

def extract_inds(vds, keep = True):
    return