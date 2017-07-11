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

def extract_region_file(vds, region_file):
    """
    BED, INTERVAL list, variants
    """
    
def sites():

def extract_inds():
    
def remove_inds():
    
    
    
    