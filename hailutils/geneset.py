#!/usr/bin/python

import os
import sys
import re
import logging
from hail import *
from pprint import pprint

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def annotate_geneset(vds, name, gene_ids = None, gene_key = 'va.ann.canonical.gene_id', ann_label = 'canonical'):
    if isinstance(gene_ids, set):
        vds = vds.annotate_global('global.{}'.format(name), gene_ids, TSet(TString()))
    vds = (
        vds
            .annotate_variants_expr(
                '''
                va.genesets.{0}.{1} = 
                    if (! {2}.filter(x => global.{1}.contains(x)).isEmpty())
                        true
                    else
                        false
                '''.format(ann_label, name, gene_key)
        )
    )
    return(vds)

def annotate_genesets_dict(vds, geneset_dict, gene_key = 'va.ann.canonical.gene_id', ann_label = 'canonical'):
    for key in geneset_dict:
        vds = vds.annotate_global('global.{}'.format(key), geneset_dict[key], TSet(TString()))
    exprs = [
                '''
                va.genesets.{0}.{1} = 
                    if (! {2}.filter(x => global.{1}.contains(x)).isEmpty())
                        true
                    else
                        false
                '''.format('canonical', key, 'va.ann.canonical.gene_id')
                for key in geneset_dict
        
    ]
    vds = vds.annotate_variants_expr(exprs)
    return(vds)

def annotate_nonpsych_lof_intolerant_genes(vds, gene_key = 'va.ann.canonical.gene_id', ann_label = 'canonical'):
    gene_ids = (
        vds.hc
            .read_table('gs://exome-qc/resources/exac_release0.3.1/constraint/exac_constraint_nonpsych.kt')
            .filter('pLI >= 0.9')
            .query('gene_id.collect().toSet()')
    )
    
    vds = annotate_geneset(vds, name = 'nonpsych_lof_intolerant_genes', gene_ids = gene_ids, gene_key = gene_key, ann_label = ann_label)
    return(vds)

def annotate_lof_intolerant_genes(vds, gene_key = 'va.ann.canonical.gene_id', ann_label = 'canonical'):
    gene_ids = (
        vds.hc
            .read_table('gs://exome-qc/resources/exac_release0.3.1/constraint/exac_constraint.kt')
            .filter('pLI >= 0.9')
            .query('gene_id.collect().toSet()')
    )
    
    vds = annotate_geneset(vds, name = 'lof_intolerant_genes', gene_ids = gene_ids, gene_key = gene_key, ann_label = ann_label)
    return(vds)

def annotate_geneset_all_annotations(vds, name, gene_ids = None):
    vds = annotate_geneset(vds, name, gene_ids = gene_ids, gene_key = 'va.ann.canonical.gene_id', ann_label = 'canonical')
    vds = annotate_geneset(vds, name, gene_ids = gene_ids, gene_key = 'va.ann.basic.gene_id', ann_label = 'basic')
    vds = annotate_geneset(vds, name, gene_ids = gene_ids, gene_key = 'va.ann.all.gene_id', ann_label = 'all')
    return(vds)
