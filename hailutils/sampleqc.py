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
    return(vds.pc_relate(k = k, maf = maf, block_size = block_size, min_kinship = min_kin))

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

def prune_samples(kinship_kt, sample_kt = None, tiebreak_expr = None, min_k = None):
    if min_k:
        logger.info('Retain only pairwise relationships with kin >= %s', min_k)
        kinship_kt = kinship_kt.filter("kin >= {}".format(min_k))
    
    logger.info('Define pairwise samples as sets.')
    related_samples = kinship_kt.query('i.flatMap(i => [i,j]).collectAsSet()')
    
    if (tiebreak_expr == 'case'):
        ann = 'isCase: `sa.isCase`'
        tiebreak_expr = '''
            if (l.isCase && !r.isCase) 
                -1
            else if (!l.isCase && r.isCase) 
                1 
            else 
                0
            '''
    elif (tiebreak_expr == 'score'):
        ann = 'isCase: `sa.isCase`, score: `sa.score`'
        tiebreak_expr = '''
            if (l.isCase && !r.isCase) 
                -1
            else if (!l.isCase && r.isCase) 
                1 
            else if (l.score > r.score)
                -1
            else if (l.score < r.score)
                1
            else
                0
            '''
    
    logger.info('Identify maximum independent set of samples.')
    if (tiebreak_expr):
        logger.info('Tiebreak expression. %s', tiebreak_expr)
        related_samples_to_keep = (
            kinship_kt
                .key_by("i")
                .join(sample_kt)
                .annotate('iAndCase = { id: i, %s }' % ann)
                .select(['j', 'iAndCase'])
                .key_by("j")
                .join(sample_kt)
                .annotate('jAndCase = { id: j, %s }' % ann)
                .select(['iAndCase', 'jAndCase'])
                .maximal_independent_set(
                    "iAndCase", 
                    "jAndCase",
                    tie_breaker = tiebreak_expr
                )
        )
        related_samples_to_remove = related_samples - {x.id for x in related_samples_to_keep}
    else:
        related_samples_to_keep = kinship_kt.maximal_independent_set("i", "j")
        related_samples_to_remove = related_samples - set(related_samples_to_keep)
        
    logger.info('Of the %s samples observed in the kinship matrix, we retain %s and exclude %s.', len(related_samples), len(set(related_samples_to_keep)), len(related_samples_to_remove))

    return(related_samples_to_remove)

def sex_check(vds):
    logger.info('Read in 1000 Genomes variants (with sex chromosomes).')
    kt = (
        vds.hc
            .read_table('gs://exome-qc/resources/1000-genomes-phase3/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.snps.maf005.sites.kt')
            .filter('v.contig == "22" || v.contig == "X" || v.contig == "Y"')
    )

    logger.info('Filter to 1000 Genomes variants.')
    vds = vds.split_multi().filter_variants_table(kt, keep = True).repartition(500)
     
    logger.info('Variants in each chromosome: %s', vds.query_variants('variants.map(v => v.contig).counter()'))

    logger.info('Apply variant filters to chr22, chrX, and chrY variants.')
    vds = (
        vds
            .variant_qc()
            .filter_variants_expr(
                '''
                (
                    v.contig == "22" &&
                    va.qc.callRate >= 0.97 && va.qc.AF >= 0.005 && 
                    va.qc.pHWE >= 1e-05 && 
                    va.qc.dpMean >= 7 && va.qc.dpMean <= 100 && 
                    va.qc.dpMean/va.qc.dpStDev >= 0.5
                ) || 
                (
                    v.contig == "X" &&
                    va.qc.callRate >= 0.97 && va.qc.AF >= 0.005 && 
                    va.qc.dpMean >= 7 && va.qc.dpMean <= 100
                ) ||
                (
                    v.contig == "Y" &&
                    va.qc.dpMean >= 3.5
                )   
                '''     
            )
    )
    logger.info('Variants in each chromosome: %s', vds.query_variants('variants.map(v => v.contig).counter()'))

    logger.info('Calculate depth metrics for chr22, chrX, and chrY.')
    vds = (
        vds.annotate_samples_expr(
            [
                'sa.dpMean_chr22 = gs.filter(g => v.contig == "22").map(g => g.dp).stats().mean',
                'sa.dpMean_chrXPar = gs.filter(g => v.inXPar()).map(g => g.dp).stats().mean',
                'sa.dpMean_chrYPar = gs.filter(g => v.inYPar()).map(g => g.dp).stats().mean',
                'sa.dpMean_chrXNonPar = gs.filter(g => v.inXNonPar()).map(g => g.dp).stats().mean',
                'sa.dpMean_chrYNonPar = gs.filter(g => v.inYNonPar()).map(g => g.dp).stats().mean',
                'sa.dpZscoreMean_chr22 = gs.filter(g => v.contig == "22").map(g => (g.dp - va.qc.dpMean)/va.qc.dpStDev).stats().mean',
                'sa.dpZscoreMean_chrX = gs.filter(g => v.inXNonPar()).map(g => (g.dp - va.qc.dpMean)/va.qc.dpStDev).stats().mean',
                'sa.dpZscoreMean_chrY = gs.filter(g => v.inYNonPar()).map(g => (g.dp - va.qc.dpMean)/va.qc.dpStDev).stats().mean'
            ]
        )
    )

    logger.info('Impute sex.')
    vds = vds.impute_sex()

    logger.info('Return sample keytable.')
    kt = vds.samples_table().flatten()

    kt = prettify_columns(kt, ['sa', 'sa.imputesex'])
    return(kt)