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

os.environ['CLOUD'] = 'dataproc'

def ld_prune(vds, r2 = 0.1, window = 1000000, memory_per_core = 512, num_cores = 90):    
    logger.info('Exclude variants in high LD regions.')
    kt = KeyTable.import_interval_list(
        '/user/tsingh/resources/price_2010-high_ld_regions.tsv' 
        if os.environ['CLOUD'] == 'local' else 'gs://exome-qc/resources/high-ld-regions/price_2010-high_ld_regions.tsv'
    )
    vds = vds.filter_variants_table(kt, keep = False)
    logger.info('After excluding high LD regions, %s samples and %s variants observed.', *vds.count())
    logger.info('Begin LD pruning.')
    vds = vds.ld_prune(r2 = r2, window = window, memory_per_core = memory_per_core, num_cores = num_cores)
    logger.info('After LD pruning, %s samples and %s variants observed.', *vds.count())  
    return(vds)

def filter_and_ld_prune(
    vds, sample_ids = None, 
    pHWE_threshold = 1e-05, callrate_threshold = 0.99,
    af_threshold = 0.05, meandp_threshold = 7,
    **kwargs
):
    if sample_ids:
        logger.info('Filter to %s samples.', len(sample_ids))  
        vds = vds.filter_samples_list(sample_ids, keep = True)

    logger.info('Run sample and variant QC.')
    vds = vds.sample_qc('sa.rqc').variant_qc('va.rqc')

    logger.info(
        'Restrict analysis to variants with: pHWE >= %s, callRate >= %s, AF >= %s, dpMean >= %s.',
        pHWE_threshold, callrate_threshold, af_threshold, meandp_threshold) 

    vds = (
        vds.filter_variants_expr(
            '''
            v.altAllele.isSNP() &&
            va.rqc.pHWE >= {} && 
            va.rqc.callRate >= {} && 
            va.rqc.AF >= {} && 
            va.qc.dpMean >= {}
            '''.format(pHWE_threshold, callrate_threshold, af_threshold, meandp_threshold)
        )
    )
    logger.info('After variant filtering, we observe %s samples and %s variants', *vds.count())
    vds = ld_prune(vds, **kwargs)    
    return(vds)

def pc_relate_subset(
    vds, sample_ids, 
    min_kin = 0.04,
    k = 5, maf = 0.05, block_size = 1024,
    **kwargs
):
    vds = filter_and_ld_prune(vds, sample_ids, **kwargs)
    logger.info(
        'Run PC-Relate with the following parameters: k = %s, maf = %s, and blocksize = %s.', 
        k, maf, block_size)
    return(
       vds
            .pc_relate(k = k, maf = maf, block_size = block_size)
            .filter('kin >= {}'.format(min_kin))
    )

def pc_project(vds, pc_vds, pca_loadings_root = 'va.pca_loadings'):
    """
    Projects samples in `vds` on PCs computed in `pc_vds`
    :param vds: VDS containing the samples to project
    :param pc_vds: VDS containing the PC loadings for the variants
    :param pca_loadings_root: Annotation root for the loadings. Can be either an Array[Double] or a Struct{ PC1: Double, PC2: Double, ...}
    :return: VDS with
    """
    pc_vds = pc_vds.annotate_variants_expr('va.pca.calldata = gs.callStats(g => v)')

    pcs_struct_to_array = ",".join(['vds.pca_loadings.PC%d' % x for x in range(1, 21)])
    arr_to_struct_expr = ",".join(['PC%d: sa.pca[%d - 1]' % (x, x) for x in range(1, 21)])

    vds = (
        vds
            .filter_multi()
            .annotate_variants_vds(
                pc_vds, 
                expr = 'va.pca_loadings = [%s], va.pca_af = vds.pca.calldata.AF[1]' % pcs_struct_to_array
            )
            .filter_variants_expr('!isMissing(va.pca_loadings) && !isMissing(va.pca_af)')
    )
    n_variants = vds.query_variants(['variants.count()'])[0]

    return(
        vds
            .annotate_samples_expr(
                'sa.pca = gs.filter(g => g.isCalled && va.pca_af > 0.0 && va.pca_af < 1.0).map(g => let p = va.pca_af in (g.gt - 2 * p) / sqrt(%d * 2 * p * (1 - p)) * va.pca_loadings).sum()' % n_variants
            )
           .annotate_samples_expr('sa.pca = {%s}' % arr_to_struct_expr)
    )

def autopca(vds, sample_ids = None, k = 20, **kwargs):
    return(
        filter_and_ld_prune(vds, sample_ids, **kwargs)
            .pca(scores = 'sa.scores', loadings = 'va.pca_loadings', k = k)
    )
        
