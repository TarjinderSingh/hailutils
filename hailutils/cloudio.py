#!/usr/bin/python

import os
import sys
import re
import logging
from pprint import pprint
import csv

from hail import *

# Basic IO from text files 
# Ported from glists.py in pyrunner

def parse(filepath, split=False, sep=None, cols=None, skip = None, strip = True):
    """ skip is 1-based, and represents the actual lines skipped. ie. 3 means skip lines 1,2,3 """
    with hadoop_read(filepath) as f:
        lines = f.readlines()
        if strip:
            lines = [ line.strip() for line in lines]
        if (isinstance(skip, int)) and (skip > 0):
            lines = lines[skip:]
        if split or sep or cols or (cols == 0):
            lines = [ line.split(sep) for line in lines ]
        if cols or (cols == 0):
            if isinstance(cols, int):
                lines = [ line[cols] for line in lines ]
            elif isinstance(cols, list):
                all_lines, lines = lines, []
                for line in all_lines:
                    lines.append([ line[i] for i in cols ])
    return(lines)

def read_list(filepath):
    return(parse(filepath))

def write_nested_list(outfile, nested_list):
    with hadoop_write(outfile) as f:   
        writer = csv.writer(f, delimiter='\t', lineterminator = '\n')
        for l in nested_list:
            writer.writerow(l)

def write_list(outfile, l):
    with hadoop_write(outfile) as f:
        writer = csv.writer(f, delimiter = '\n')
        writer.writerow(l)