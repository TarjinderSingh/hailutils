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
    """Log case counts in VDS.
    
    Parameters
    ----------
    vds : Hail variant data set
    """
    logger.info(
        '%s cases and %s controls observed in the data set.',
        vds.query_samples('samples.map(s => sa.pheno.isCase).map(x => x.toInt()).sum()'),
        vds.query_samples('samples.map(s => ! sa.pheno.isCase).map(x => x.toInt()).sum()')
    )
    
def log_sex_count(vds):
    """Log male and female counts in VariantDataset.
    
    Parameters
    ----------
    vds : Variant Data Set format from Hail.
    """ 
    logger.info(
        '%s females and %s males observed in the data set.',
        vds.query_samples('samples.map(s => sa.pheno.isFemale).map(x => x.toInt()).sum()'),
        vds.query_samples('samples.map(s => ! sa.pheno.isFemale).map(x => x.toInt()).sum()')
    )

def exact_keytable(vds, test = 'fisher', minCellCount = 5):
    """Calculates case-control counts for each variant and a P-value using the exact test or Chi-sq test.
    
    Parameters
    ----------
    vds : TYPE
        Description
    test : str, optional
        'fisher', 'ctt', or 'both' depending method to be used.
    minCellCount : int, optional
        Minimum cell counts before Chi-sq test is used over the Exact test.
    
    Returns
    -------
    VariantDataset
        Annotated with va.exact in the variant schema.
    """
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

def get_results_table(
    vds, 
    filter_expr = None,
    keep_cols = [
        'Variant',
        'va.rsid',
        'va.info.VQSLOD',
        'va.prefilt.vqc.callRate',
        'va.prefilt.vqc.AC',
        'va.prefilt.vqc.dpMean',
        'va.prefilt.vqc.pHWE',
        'va.prefilt.hetGQ.mean',
        'va.prefilt.hetAB.mean',
        'va.prefilt.hetDP.mean',
        'va.postfilt.vqc.callRate',
        'va.firth.beta',
        'va.firth.chi2',
        'va.firth.pval',
        'va.firth.fit.converged',
        'va.wald.beta',
        'va.wald.se',
        'va.wald.zstat',
        'va.wald.pval',
        'va.wald.fit.converged',
        'va.exact.Xcase',
        'va.exact.Xctrl',
        'va.exact.Ncase',
        'va.exact.Nctrl',
        'va.exact.AFcase',
        'va.exact.AFctrl',
        'va.exact.ctt.pValue',
        'va.exact.ctt.oddsRatio'
    ],
    prettify_headers = ['va', 'va.exact']):
    
    if filter_expr:
        logger.info('Filter results VDS based on expression: %s', filter_expr)
        vds = vds.filter_variants_expr(filter_expr)
    
    logger.info('Return results as KeyTable.')
    kt = (
        vds 
            .variants_table()
            .flatten()
            .annotate('Variant = str(v)')
            .select(keep_cols)
    )
    kt = prettify_columns(kt, strip_match = prettify_headers)
    return(kt)

# filt_expr = '`va.postfilt_metrics.AC` <= 5 && (isMissing(`va.nonpsych_gnomad.AC`) || `va.nonpsych_gnomad.AC` <= 5)'

def nonref_genotypes_table(vds, filt_expr = None, explode_col = None):
    '''
    filt_expr: '`va.postfilt_metrics.AC` <= 5 && (isMissing(`va.nonpsych_gnomad.AC`) || `va.nonpsych_gnomad.AC` <= 5)'
    explode_col: 'va.ann.canonical.gene_id'
    '''
    kt = (
        vds
            .genotypes_table()
            .filter('g.gt > 0')
            .flatten()
    )
    if filt_expr:
        kt = kt.filter(filt_expr)
    if explode_col:
        kt = kt.explode(explode_col)
    return(kt)

def sample_burden(kt, sample_col = 's', phe_col = 'sa.isCase', group_col = 'va.ann.canonical.gene_id', batch_list = None, csq_col = 'csq2', func = 'binary'):
    '''
    fun: binary or sum
    '''
    expr = None
    if func == 'binary':
        agg_expr = 'X = (g.map(g => g.gt).sum() > 0).toInt()'
    elif func == 'sum':
        agg_expr = 'X = g.map(g => g.gt).sum().toInt()'
    else:
        return(None)
    
    keys_expr = [
        'Sample = `{}`'.format(sample_col),
        'phe = `{}`'.format(phe_col),
        'group = `{}`'.format(group_col),
        'csq = `{}`'.format(csq_col)
    ]
    
    if batch_list: 
        keys_expr.extend([ '`{0}` = `{0}`'.format(b) for b in batch_list ])
    
    return(kt.aggregate_by_key(keys_expr, agg_expr))

def grouped_burden(kt, keys_list = [ 'group', 'phe', 'csq', 'batch' ], agg_expr = 'X = X.sum()'):
    keys_expr = [ '`{0}` = `{0}`'.format(k) for k in keys_list ]
    return(kt.aggregate_by_key(keys_expr, agg_expr))

def quo(x):
    return('`{}`'.format(x))

def annotate_csq2_expr(new_col = 'va.csq2', csq_col = 'va.ann.canonical.csq', mpc_col = 'va.mpc.MPC', is_kt = True):
    if is_kt:
        new_col = quo(new_col)
        csq_col = quo(csq_col)
        mpc_col = quo(mpc_col)
    return('{0} = if (isDefined({1}) && {1} >= 2) "mis2" else {2}'.format(new_col, mpc_col, csq_col))

def annotate_mpc_expr(new_col = 'va.csq2', csq_col = 'va.ann.canonical.csq', mpc_col = 'va.mpc.MPC', is_kt = True):
    if is_kt:
        new_col = quo(new_col)
        csq_col = quo(csq_col)
        mpc_col = quo(mpc_col)
    expr = '''
    {0} = 
        if (isDefined({1}) && {1} >= 3 && {2} == "mis") 
            "mis3"
        else if (isDefined({1}) && {1} >= 2 && {2} == "mis") 
            "mis2"
        else if (isDefined({1}) && {1} >= 1.5 && {2} == "mis") 
            "mis15"
        else if (isDefined({1}) && {1} >= 1 && {2} == "mis") 
            "mis1"    
        else 
            {2}
    '''.format(new_col, mpc_col, csq_col)
    return(expr)

def annotate_startlost_expr(new_col = 'va.csq2', csq_col = 'va.ann.canonical.csq'):
    if is_kt:
        new_col = quo(new_col)
        csq_col = quo(csq_col)
    return('{0} = if (isDefined({1}) && {1} == "startlost") "mis" else {1}'.format(new_col, csq_col))

def annotate_loftee_expr(new_col = 'csq2', csq_col = 'va.ann.canonical.csq', loftee_col = 'va.ann.canonical.loftee', is_kt = True):
    if is_kt:
        new_col = quo(new_col)
        csq_col = quo(csq_col)
        loftee_col = quo(loftee_col)
    expr = '''
    {0} = 
        if (isDefined({1}) && {1} == "HC") 
            "lofHC"
        else 
            {2}
    '''.format(new_col, loftee_col, csq_col)
    return(expr)

def annotate_analysis_consequences_keytable(kt, anns = [ 'canonical', 'basic', 'all' ]):
    for ann in anns:
        expr = annotate_mpc_expr(new_col = 'va.csq.{}'.format(ann), csq_col = 'va.ann.{}.csq'.format(ann), mpc_col = 'va.mpc.MPC', is_kt = True)
        logger.info('Annotate kt with expression: %s', expr)
        kt = kt.annotate(expr)
        expr = annotate_loftee_expr(new_col = 'va.csq.{}'.format(ann), csq_col = 'va.csq.{}'.format(ann), loftee_col = 'va.ann.{}.loftee'.format(ann), is_kt = True)
        logger.info('Annotate kt with expression: %s', expr)
        kt = kt.annotate(expr)
    return(kt)

def annotate_analysis_consequences_vds(vds, anns = [ 'canonical', 'basic', 'all' ]):
    for ann in anns:
        expr = annotate_mpc_expr(new_col = 'va.csq.{}'.format(ann), csq_col = 'va.ann.{}.csq'.format(ann), mpc_col = 'va.mpc.MPC', is_kt = False)
        logger.info('Annotate VDS with expression: %s', expr)
        vds = vds.annotate_variants_expr(expr)
        expr = annotate_loftee_expr(new_col = 'va.csq.{}'.format(ann), csq_col = 'va.csq.{}'.format(ann), loftee_col = 'va.ann.{}.loftee'.format(ann), is_kt = False)
        logger.info('Annotate VDS with expression: %s', expr)
        vds = vds.annotate_variants_expr(expr)
    return(vds)

def annotate_noncoding_expr(new_col = 'va.csq2', score_col = 'va.ann.canonical.minscore', csq_col = 'va.ann.canonical.csq', is_kt = True):
    if is_kt:
        new_col = quo(new_col)
        csq_col = quo(csq_col)
        score_col = quo(score_col)
    expr = '''
    {0} = 
        if (isDefined({1}) && isMissing({2}) && ({1} == 18 || {1} == 19))
            "utr"
        else if (isDefined({1}) && isMissing({2}) && {1} == 21)
            "intron"
        else if (isDefined({1}) && isMissing({2}) && ({1} == 24 || {1} == 25))
            "flank"
        else
            {2}
    '''.format(new_col, score_col, csq_col)
    return(expr)