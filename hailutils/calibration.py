#!/usr/bin/env python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *
import pandas as pd

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def annotate_gatk_bundle(vds):
    logger.info('Annotate with GATK bundle.')
    logger.info('Annotate with HapMap sites.')
    vds = vds.annotate_variants_vds(
        vds.hc
            .read('gs://exome-qc/resources/gatkbundle_b37/hapmap_3.3.b37.vds')
            .annotate_variants_expr('va.inHapMap.isTrue = 1, va.inHapMap.AF = va.info.AF'),
        expr = 'va.inHapMap = vds.inHapMap'
    )
    logger.info('Annotate with Omni2.5 sites.')
    vds = vds.annotate_variants_vds(
        vds.hc
            .read('gs://exome-qc/resources/gatkbundle_b37/1000G_omni2.5.b37.vds')
            .annotate_variants_expr('va.inOmni25 = 1'),
        expr = 'va.inOmni25 = vds.inOmni25'
    )
    logger.info('Annotate with 1000 Genomes sites.')
    vds = vds.annotate_variants_vds(
        vds.hc
            .read('gs://exome-qc/resources/gatkbundle_b37/1000G_phase3_v4_20130502.snvs.sites.vds')
            .annotate_variants_expr('va.in1000G.isTrue = 1, va.in1000G.AF = va.info.AF'),
        expr = 'va.in1000G = vds.in1000G'
    )
    logger.info('Annotate with Mills sites.')
    vds = vds.annotate_variants_vds(
        vds.hc
            .read('gs://exome-qc/resources/gatkbundle_b37/Mills_and_1000G_gold_standard.indels.b37.vds')
            .annotate_variants_expr('va.inMills = 1'),
        expr = 'va.inMills = vds.inMills'
    )
    return(vds)

def gatk_feature_exprs(root = "va.features"):
    return([ 
        'va.{0}.{1} = va.info.{1}'.format(root, ann) for ann in 
        [ 'InbreedingCoeff', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum' ] 
    ])

def allele_feature_exprs(root = "va.features", sample_filt_expr = ""):
    sample_filt_expr = "&& " + sample_filt_expr if sample_filt_expr else ""
    exprs = [
        '{}.nrGQ = gs.filter(g => g.isCalledNonRef {}).map(g => g.gq).stats()',
        '{}.nrDP = gs.filter(g => g.isCalledNonRef {}).map(g => g.dp).stats()',
        '{}.nrq = gs.filter(g => g.isCalledNonRef {}).map(g => -log10(g.gp[0])).stats()',
        '{}.AB = gs.filter(g => g.isHet {}).map(g => g.ad[1]/g.dp).stats()',
        '{}.optAB = gs.filter(g => g.isHet {}).map(g => abs((g.ad[1]/g.dp) - 0.5)).min()',
        '{}.pAB = gs.filter(g => g.isHet {}).map(g => g.pAB()).stats()',
        '{}.nrDPsum = gs.filter(g => g.isCalledNonRef {}).map(g => g.dp).sum()',
        '{}.QUAL = -10 * gs.filter(g => g.isCalledNonRef {}).map(g => if (g.pl[0] > 3000) -300 else log10(g.gp[0])).sum()'
    ]
    return([ x.format(root, sample_filt_expr) for x in exprs ])

def allele_feature_exprs2(root = "va.features"):
    exprs = [
        '{0}.QD = {0}.QUAL/{0}.nrDPsum',
        '{0}.nrGQmean = {0}.nrGQ.mean',
        '{0}.ABmean = {0}.AB.mean'
    ]
    return([ x.format(root) for x in exprs ])

def site_feature_exprs(root = "va.features"):
    exprs = [
        '{}.isSNP = v.altAlleles.map(a => (! a.isSNP()).toInt()).sum() == 0',
        '{}.nAlt = v.nAltAlleles()',
        '{}.nIndels = v.altAlleles.map(a => a.isIndel().toInt()).sum()',
        '{}.mixed = (v.altAlleles.map(a => a.isSNP().toInt()).sum() > 0) && (v.altAlleles.map(a => a.isIndel().toInt()).sum() > 0)',
        '{}.inStar = v.altAlleles.map(a => a.isStar().toInt()).sum() > 0'
    ]
    return([ x.format(root) for x in exprs ])

def hard_filter_exprs(root = 'va.prior'):
    logger.info('Identify low-quality sites based on GATK hard filters.')
    return([ 
        '''
        {}.failed_gatk_hard_filters = 
            if (va.info.QD < 2 || va.info.FS > 60 || va.info.MQ < 30)
                true
            else
                false
        '''.format(root)
    ])

rf_features = [ 
    'va.features.nrGQmean', 'va.features.ABmean', 'va.features.QD',
    'va.features.isSNP', 'va.features.nAlt', 'va.features.nIndels',
    'va.features.mixed', 'va.features.inStar',
    'va.info.InbreedingCoeff', 'va.info.FS', 'va.info.MQRankSum',
    'va.info.ReadPosRankSum', 'va.info.MQ'
]

def extract_rf_keytable(vds):
    return(
        vds
            .variants_table()
            .flatten()
            .annotate('v = str(v)')
            .select([ 'v', 'va.in1000G.AF', 'va.inMills', 'va.info.QD' ] + rf_features)
    )

def prepare_vds_for_calibration(vds):
    vds = annotate_gatk_bundle(vds)
    logger.info('Annotate with site and allele features for variant calibration.')
    return(
        vds
            .annotate_variants_expr(site_feature_exprs())
            .split_multi()
            .annotate_variants_expr(allele_feature_exprs() + allele_metrics_exprs())
            .annotate_variants_expr(allele_feature_exprs2() + allele_metrics_exprs2())
            .annotate_variants_expr(hard_filter_exprs())
    )

def concordance(left_vds, right_vds):
    summary, samples, variants = left_vds.concordance(right_vds)
    summary = summarize_concordance(summary)
    return(summary, samples, variants)

def summarize_concordance(summary):
    concordance_table = pd.DataFrame(summary)
    concordance_table.columns = [
        'No Data (missing variant)',
        'No Call (missing genotype call)',
        'Hom Ref',
        'Heterozygous',
        'Hom Var'
    ]
    concordance_table.index = [
        'No Data (missing variant)',
        'No Call (missing genotype call)',
        'Hom Ref',
        'Heterozygous',
        'Hom Var'
    ]
    return(concordance_table)

def calculate_sample_concordance(skt):
    return(
        skt
            .rename({'concordance': 'm'})
            .annotate([
                'TP = m.map(x => x[3])[3] + m.map(x => x[4])[4]', # true positives
                'P = m.map(x => x[3]).sum() + m.map(x => x[4]).sum()', # total true variants, TP + FN
                'FP = m.map(x => x[2])[3] + m.map(x => x[2])[4] + m.map(x => x[3])[4]', # false positives
                'N = m.map(x => x[2]).sum()', # total true variants, TP + FN

            ])
            .annotate('TPR = TP/P, FPR = FP/N')
            .select(['s', 'nDiscordant', 'TP', 'P', 'TPR', 'FP', 'N', 'FPR'])
    )

def calculate_variant_concordance(vkt):
    return(
        vkt
            .rename({'concordance': 'm'})
            .annotate(
                [
                    'nConcordant = m.map(x => x[2])[2] + m.map(x => x[3])[3] + m.map(x => x[4])[4]',
                    'nConcordantNonRef = m.map(x => x[3])[3] + m.map(x => x[4])[4]',
                    'nLRefRHet = m.map(x => x[3])[2]',
                    'nLHetRRef = m.map(x => x[2])[3]',
                    'nLCalled = ' + ' + '.join(
                        [ 'm.map(x => x[{}])[{}]'.format(i, j) for i in range(0, 4 + 1) for j in range(2, 4 + 1) ]),
                    'nLHet = ' + ' + '.join([ 'm.map(x => x[{}])[3]'.format(i) for i in range(0, 4 + 1) ]),
                    'nLHomVar = ' + ' + '.join([ 'm.map(x => x[{}])[4]'.format(i) for i in range(0, 4 + 1) ]),
                    'nRCalled = m.map(x => x[2]).sum() + m.map(x => x[3]).sum() +  m.map(x => x[4]).sum()',
                    'nRHet = m.map(x => x[3]).sum()',
                    'nRHomVar = m.map(x => x[4]).sum()'
                ]
            )
            .drop('m')
    )