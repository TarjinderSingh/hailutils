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

def annotate_gatk_bundle(vds):
    logger.info('Annotate with GATK bundle.')
    logger.info('Annotate with HapMap sites.')
    vds = vds.annotate_variants_vds(
        hc
            .read('gs://exome-qc/resources/gatkbundle_b37/hapmap_3.3.b37.vds')
            .annotate_variants_expr('va.inHapMap.isTrue = 1, va.inHapMap.AF = va.info.AF'),
        expr = 'va.inHapMap = vds.inHapMap'
    )
    logger.info('Annotate with Omni2.5 sites.')
    vds = vds.annotate_variants_vds(
        hc
            .read('gs://exome-qc/resources/gatkbundle_b37/1000G_omni2.5.b37.vds')
            .annotate_variants_expr('va.inOmni25 = 1'),
        expr = 'va.inOmni25 = vds.inOmni25'
    )
    logger.info('Annotate with 1000 Genomes sites.')
    vds = vds.annotate_variants_vds(
        hc
            .read('gs://exome-qc/resources/gatkbundle_b37/1000G_phase3_v4_20130502.snvs.sites.vds')
            .annotate_variants_expr('va.in1000G.isTrue = 1, va.in1000G.AF = va.info.AF'),
        expr = 'va.in1000G = vds.in1000G'
    )
    logger.info('Annotate with Mills sites.')
    vds = vds.annotate_variants_vds(
        hc
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

def allele_metrics_exprs_2(root = "va.metrics", sample_filt_expr = ""):
    sample_filt_expr = "&& " + sample_filt_expr if sample_filt_expr else ""
    exprs = [
        '{0}.callRate = {0}.nCalled/{0}.nSample',
        '{0}.AN = {0}.nCalled * 2'
    ]
    return([ x.format(root, sample_filt_expr) for x in exprs ])

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