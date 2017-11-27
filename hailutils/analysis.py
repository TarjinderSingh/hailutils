#!/usr/bin/python

import os
import sys
import re
import logging
from hail import *
from pprint import pprint

from keytable import *

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def log_case_count(vds):
    logger.info(
        '%s cases and %s controls observed in the data set.',
        vds.query_samples('samples.map(s => sa.pheno.isCase).map(x => x.toInt()).sum()'),
        vds.query_samples('samples.map(s => ! sa.pheno.isCase).map(x => x.toInt()).sum()')
    )
    
def log_sex_count(vds):
    logger.info(
        '%s females and %s males observed in the data set.',
        vds.query_samples('samples.map(s => sa.pheno.isFemale).map(x => x.toInt()).sum()'),
        vds.query_samples('samples.map(s => ! sa.pheno.isFemale).map(x => x.toInt()).sum()')
    )

def exact_keytable(vds, test = 'fisher', minCellCount = 5):
    logger.info('Annotate variants with counts.')
    vds = vds.annotate_variants_expr(
        [
            'va.exact.Xcase = gs.filter(g => sa.pheno.isCase).map(g => g.nNonRefAlleles()).sum()',
            'va.exact.Xctrl = gs.filter(g => !sa.pheno.isCase).map(g => g.nNonRefAlleles()).sum()',
            'va.exact.Ncase = gs.filter(g => sa.pheno.isCase).map(g => g.isCalled()).count()',
            'va.exact.Nctrl = gs.filter(g => !sa.pheno.isCase).map(g => g.isCalled()).count()',
            'va.exact.N = gs.map(g => g.isCalled()).count()',
            'va.exact.Xrcase = gs.filter(g => sa.pheno.isCase).map(g => 2 - g.nNonRefAlleles()).sum()',
            'va.exact.Xrctrl = gs.filter(g => !sa.pheno.isCase).map(g => 2 - g.nNonRefAlleles()).sum()',
        ]
    )

    logger.info('Run the following exact tests: %s', test)
    exprs = ['va.exact.AFcase = va.exact.Xcase/(2 * va.exact.Ncase)', 'va.exact.AFctrl = va.exact.Xctrl/(2 * va.exact.Nctrl)']
    if test == 'fisher' or test == 'both':
        exprs.append('va.exact.fet = fet(va.exact.Xcase, va.exact.Xctrl, va.exact.Xrcase, va.exact.Xrctrl)')
    if test == 'ctt' or test == 'both':
        exprs.append('va.exact.ctt = ctt(va.exact.Xcase, va.exact.Xctrl, va.exact.Xrcase, va.exact.Xrctrl, {})'.format(minCellCount))
    return(
        vds
            .annotate_variants_expr(exprs)
            .annotate_variants_expr('va = select(va, exact)')
            .variants_table()
    )

def pc_list(n, root = 'sa.covar'):
    return(['{}.PC{}'.format(root, i) for i in range(1, n + 1)])

def logreg_keytable(
    vds, method = 'wald', root = 'va.wald',
    vfilt_expr = None,
    case_label = 'sa.pheno.isCase', covar = pc_list(10)
):
    logger.info('Run logistic regression (method: %s).', method)
    if vfilt_expr:
        logger.info('Filter variants based on expression: %s', vfilt_expr)
        vds = vds.filter_variants_expr(vfilt_expr)
    kt = (
        vds
            .logreg(
                test = method, 
                y = case_label, 
                covariates = covar,
                root = root
            )
            .annotate_variants_expr('va = select(va, {})'.format(re.sub('va\.', '', root)))
            .variants_table()
    )
    return(kt)

def firth_keytable(vds, method = 'firth', root = 'va.firth', *args, **kwargs):
    return(logreg_keytable(vds, method = method, root = root, *args, **kwargs))
