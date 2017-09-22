#!/usr/bin/env python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def allele_metrics_exprs(root = "va.metrics", sample_filt_expr = ""):
    sample_filt_expr = "&& " + sample_filt_expr if sample_filt_expr else ""
    exprs = [
        '{}.nSample = gs.filter(g => true {}).count()',
        '{}.nCalled = gs.filter(g => g.isCalled() {}).count()',
        '{}.AC = gs.filter(g => g.isCalled() {}).map(g => g.nNonRefAlleles()).sum()',
        '{}.homrefGQ = gs.filter(g => g.isHomRef() {}).map(g => g.gq).stats()',
        '{}.hetGQ = gs.filter(g => g.isHet() {}).map(g => g.gq).stats()',
        '{}.homvarGQ = gs.filter(g => g.isHomVar() {}).map(g => g.dp).stats()',
        '{}.homrefDP = gs.filter(g => g.isHomRef() {}).map(g => g.dp).stats()',
        '{}.hetDP = gs.filter(g => g.isHet() {}).map(g => g.dp).stats()',
        '{}.homvarDP = gs.filter(g => g.isHomVar() {}).map(g => g.dp).stats()',
        '{}.homrefAB = gs.filter(g => g.isHomRef() {}).map(g => g.ad[1]/g.dp).stats()',
        '{}.hetAB = gs.filter(g => g.isHet() {}).map(g => g.ad[1]/g.dp).stats()',
        '{}.homvarAB = gs.filter(g => g.isHomVar() {}).map(g => g.ad[1]/g.dp).stats()'
    ]
    return([ x.format(root, sample_filt_expr) for x in exprs ])

def allele_metrics_exprs2(root = "va.metrics", sample_filt_expr = ""):
    sample_filt_expr = "&& " + sample_filt_expr if sample_filt_expr else ""
    exprs = [
        '{0}.callRate = {0}.nCalled/{0}.nSample',
        '{0}.AN = {0}.nCalled * 2'
    ]
    return([ x.format(root, sample_filt_expr) for x in exprs ])

def genotype_filter_expr(
    minDP = 10, homrefAB = 0.1, minhetAB = 0.2, 
    minhomrefGQ = 20, minSNPGQ = 20, minindelGQ = 90):
    expr = '''
        let ab = g.ad[1] / g.ad.sum in
        g.ad.sum >= {minDP} 
        && 
        (
            (g.isHomRef && ab <= {homrefAB}) || 
            (g.isHet && ab >= {minhetAB} && ab <= 1 - {minhetAB}) || 
            (g.isHomVar && ab >= 1 - {homrefAB})
        )
        &&
        (
            (g.isHomRef && g.gq >= {minhomrefGQ}) ||
            (
                v.altAllele.isSNP && 
                ( (g.isHet && g.gq >= {minSNPGQ}) || (g.isHomVar && g.gq >= {minSNPGQ}) )
            ) ||
            (
                v.altAllele.isIndel && 
                ( (g.isHet && g.gq >= {minindelGQ}) || (g.isHomVar && g.gq >= {minindelGQ}) )
            )
        )
        '''.format(minDP = minDP, homrefAB = homrefAB, minhetAB = minhetAB, 
                   minhomrefGQ = minhomrefGQ, minSNPGQ = minSNPGQ, minindelGQ = minindelGQ)
    logger.info('We use the following expression to filter genotypes: %s', expr)
    return(expr)