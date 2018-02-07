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
        '{}.nSample = gs.filter(g => true {}).count().toInt()',
        '{}.nCalled = gs.filter(g => g.isCalled() {}).count().toInt()',
        '{}.nHomRef = gs.filter(g => g.isHomRef() {}).count().toInt()',
        '{}.nHet = gs.filter(g => g.isHet() {}).count().toInt()',
        '{}.nHomVar = gs.filter(g => g.isHomVar() {}).count().toInt()',
        '{}.AC = gs.filter(g => g.isCalled() {}).map(g => g.nNonRefAlleles()).sum().toInt()'
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

def sex_aware_allele_metrics_exprs2(root = "va.metrics", root_m = "va.metrics.male", root_f = "va.metrics.female"):
    exprs = [
        '''
        {root}.callRate = 
            if (! v.inYNonPar())
                ({root_m}.nCalled + {root_f}.nCalled)/({root_m}.nSample + {root_f}.nSample)
            else if (v.inYNonPar())
                {root_m}.nCalled/{root_m}.nSample
            else
                NA: Float
        ''',
        '''
        {root}.AC = 
            if (! v.inXNonPar() && ! v.inYNonPar())
                {root_m}.AC + {root_f}.AC
            else if (v.inXNonPar())
                {root_m}.nHomVar + {root_f}.AC
            else if (v.inYNonPar())
                {root_m}.nHomVar
            else
                NA: Int
        ''',
        '''
        {root}.AN = 
            if (! v.inXNonPar() && ! v.inYNonPar())
                ({root_m}.nCalled + {root_f}.nCalled) * 2
            else if (v.inXNonPar())
                {root_m}.nCalled + {root_f}.nCalled * 2
            else if (v.inYNonPar())
                {root_m}.nCalled
            else
                NA: Int
        '''
    ]
    return([ x.format(root = root, root_m = root_m, root_f = root_f) for x in exprs ])    

def annotate_sex_aware_allele_metrics(vds, root = 'va.metrics', sex_boolean = 'sa.isFemale', hardcall = True):
    logger.info('Calculate and annotate with sex-aware variant metrics.')
    exprs1 = (
        allele_metrics_exprs("{}.male".format(root), "! {}".format(sex_boolean), hardcall) + 
        allele_metrics_exprs("{}.female".format(root), "{}".format(sex_boolean), hardcall)
    )
    exprs2 = sex_aware_allele_metrics_exprs2(root, "{}.male".format(root), "{}.female".format(root))
    return(vds.annotate_variants_expr(exprs1).annotate_variants_expr(exprs2))

def hwe_expr(root = "va.metrics"):
    return('{0}.hwe = hwe({0}.nHomRef, {0}.nHet, {0}.nHomVar)'.format(root))

def sex_aware_hwe_expr(root = "va.metrics", root_m = "va.metrics.male", root_f = "va.metrics.female"):
    return(
        '''
        {root}.hwe = 
            if (v.isAutosomal() || v.inXPar())
                let x = hwe({root_m}.nHomRef + {root_f}.nHomRef, {root_m}.nHet + {root_f}.nHet, {root_m}.nHomVar + {root_f}.nHomVar) in {{ rExpectedHetFrequency: x.rExpectedHetFrequency, pHWE: x.pHWE }}
            else if (v.inXNonPar())
                let x = hwe({root_f}.nHomRef, {root_f}.nHet, {root_f}.nHomVar) in {{ rExpectedHetFrequency: x.rExpectedHetFrequency, pHWE: x.pHWE }}
            else
                NA: Struct{{rExpectedHetFrequency:Double, pHWE:Double}}
        '''.format(root = root, root_m = root_m, root_f = root_f)
    )

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

def sex_aware_genotype_filter_expr(
    isFemale = 'sa.pheno.isFemale', 
    minDP = 10, minDP_male_PAR = 5, 
    homrefAB = 0.1, minhetAB = 0.25,
    minhomrefGQ = 20, minSNPGQ = 20, minindelGQ = 20):
    expr = '''
        let ab = g.ad[1] / g.dp in
        (
            (! v.inXNonPar() && ! v.inYNonPar() && g.dp >= {minDP}) ||
            ( (v.inXNonPar() || v.inYNonPar()) && 
                (
                    ( isMissing({isFemale}) && g.dp >= {minDP} ) ||
                    ( {isFemale} && g.dp >= {minDP}) ||
                    ( ! {isFemale} && ((! g.isHet && g.dp >= {minDP_male_PAR}) || (g.isHet && g.dp >= {minDP})))
                )
            )
        )
        && 
        (
            (g.isHomRef && ab <= {homrefAB}) || 
            (g.isHet && ab >= {minhetAB} && ab < 1 - {homrefAB}) || 
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
        '''.format(isFemale = isFemale, 
                   minDP = minDP, minDP_male_PAR = minDP_male_PAR,
                   homrefAB = homrefAB, minhetAB = minhetAB, 
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

def gnomad_genotype_filter_expr(ADJ_GQ = 20, ADJ_DP = 10, ADJ_AB = 0.2):
    ADJ_CRITERIA = 'g.gq >= %(gq)s && g.dp >= %(dp)s && (' \
               '!g.isHet || ' \
               '(g.gtj == 0 && g.ad[g.gtk]/g.dp >= %(ab)s) || ' \
               '(g.gtj > 0 && g.ad[g.gtj]/g.dp >= %(ab)s && g.ad[g.gtk]/g.dp >= %(ab)s)' \
               ')' % {'gq': ADJ_GQ, 'dp': ADJ_DP, 'ab': ADJ_AB}
    return(ADJ_CRITERIA)
