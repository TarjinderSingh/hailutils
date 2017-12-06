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

def allele_metrics_exprs(root = "va.metrics", sample_filt_expr = "", hardcall = False):
    sample_filt_expr = "&& " + sample_filt_expr if sample_filt_expr else ""
    exprs = [
        '{}.nSample = gs.filter(g => true {}).count()',
        '{}.nCalled = gs.filter(g => g.isCalled() {}).count()',
        '{}.AC = gs.filter(g => g.isCalled() {}).map(g => g.nNonRefAlleles()).sum()'
    ]
    if not hardcall:
        exprs.extend(
            [
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
        )
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
    logger.info('The following expression is used to filter genotypes: %s', expr)
    return(expr)

def variant_filter_expr(
    min_gqmean = 20, 
    min_dpmean = 10, 
    pre_pHWE = 1e-6, 
    post_pHWE = 1e-6, 
    min_hetAB = 0.3, 
    post_callrate = 0.97, 
    min_AC = 0,
    empty_filters = True
):
    expr_list = [
        'va.prefilt.vqc.gqMean >= {}'.format(min_gqmean),
        'va.prefilt.vqc.dpMean >= {}'.format(min_dpmean),
        'va.prefilt.hetAB.mean >= {}'.format(min_hetAB),           
        'va.prefilt.hetAB.mean <= {}'.format(1 - min_hetAB),
        'va.prefilt.vqc.pHWE >= {}'.format(pre_pHWE),
        'va.postfilt.vqc.pHWE >= 1e-6'.format(post_pHWE),
        'va.postfilt.callRate >= {}'.format(post_callrate),
        'va.postfilt.AC >= {}'.format(min_AC)
    ]

    if empty_filters:
        expr_list.append('va.filters.isEmpty()')

    expr = ' && '.join(expr_list)
    logger.info('The following expression is used to filter variants: %s', expr)
    return(expr)

def basic_variant_qc(vds, root = 'va.basic_qc', sample_filt_expr = ""):
    logger.info('Perform basic variant QC.')
    if sample_filt_expr != "":
        logger.info('Restricting to samples using the following expression: %s', sample_filt_expr)
    vds = vds.annotate_variants_expr(allele_metrics_exprs(root = root, sample_filt_expr = sample_filt_expr, hardcall = True))
    return(vds)

def exclude_monoallelic(vds):
    logger.info('Calculate AC (ignoring multiallelic sites).')
    vds = vds.annotate_variants_expr(allele_metrics_exprs(hardcall = True))

    logger.info('Exclude monoallelic sites.')
    vds = vds.filter_variants_expr('va.metrics.AC > 0')
    return(vds)

def variant_qc_subset(vds, sample_expr, root = 'va.qc'):
    logger.info('')
    kt = (
        vds
        .filter_samples_expr(sample_expr)
        .variant_qc(root)
        .annotate_variants_expr('va = select(va, `{}`)'.format(root.split('.')[1]))
        .variants_table()
    )
    vds = vds.annotate_variants_table(kt, expr = 'va = merge(va, table)')
    return(vds)