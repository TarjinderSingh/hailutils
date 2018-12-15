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

def read_genesets_table(path = '/psych/genetics_data/tsingh/projects/sczexomes/analysis/release_v2.0.5/2018-07-19_merged-genesets-for-analysis.tsv'):
    '''Read gene sets from the standard gene set format into a Pandas dataframe.
        
    Parameters
    ----------
    path : {str}, optional
        Location of the gene set file.(the default is '/psych/genetics_data/tsingh/projects/sczexomes/analysis/release_v2.0.5/2018-07-19_merged-genesets-for-analysis.tsv', which [default_description])
    '''
    
    import pandas as pd

    logger.info('Read gene sets.')
    data = pd.read_table(path)
    logger.info('Clean up table.')
    data = data.drop('gene_name', axis = 1)
    data['name'] = data['name'].str.replace(':', '_').str.replace('.', '_').str.replace('-', '_')
    data = data[data.gene_id.notnull()]
    logger.info('Groups in gene set table: %s', set(data.group))
    return(data)

def generate_genesets_dict(data):
    '''Convert Pandas dataframe containing gene sets into a dictionary.
    
    The dictionary returned has the gene set name as key and a set of genes as values.
    
    Parameters
    ----------
    data : {[type]}
        [description]
    '''
    import pandas as pd
    from pyrunner.pycore import ProgressBar
    names = list(set(data.name))
    bar = ProgressBar(names, 0.01, 50)
    logger.info('Convert to dictionary.')
    geneset_dict = {}
    for name in names:
        geneset_dict[name] = set(data.query('name == "{}"'.format(name))['gene_id'].tolist())
        bar.update()
    return(geneset_dict)

def conditional_genesets_table(data, cond_name = 'gnomAD_constraint_category__pLI', label = 'pLI'):
    '''Generate conditional gene set analysis Pandas table from standard gene set input.
    
    Parameters
    ----------
    data : {[dataframe]}
        Pandas dataframe
    cond_name : {str}, optional
        Gene set used for conditional analysis (the default is 'gnomAD_constraint_category__pLI', which [default_description])
    label : {str}, optional
        String for indicating a gene set has been conditioned. (the default is 'pLI', which [default_description])
    '''
    import pandas as pd
    cond_genes = list(set(data[data.name == cond_name].gene_id.tolist()))
    
    df = data[data.name == cond_name]
    df['group'] = 'conditional'
    df['name'] = label
    
    data['in_cond'] = False
    data.loc[data.gene_id.isin(cond_genes), 'in_cond'] = True

    data['old_name'] = data.name
    data.loc[data.in_cond, 'name'] = data.loc[data.in_cond, 'name'] + '__in__' + label
    data.loc[~data.in_cond, 'name'] = data.loc[~data.in_cond, 'name'] + '__ex__' + label
    
    data = pd.concat([data, df], axis = 0)
    
    return(data)


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
                '''.format(ann_label, key, gene_key)
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

def annotate_variants_geneset(vds, geneset_dict, anns = [ 'all', 'canonical', 'basic']):
    '''Annotate variants in a VDS with gene set membership based on a hierarchy of annotations.
        
    Parameters
    ----------
    vds : {[hail.VariantDataset]}
        Variant Dataset from Hail
    geneset_dict : {[dict]}
        Dictionary of gene sets
    anns : {list}, optional
        Variant annotations for which to annotate variants with (the default is [ 'all', 'canonical', 'basic'])
    '''
    logger.info('Annotate with genesets.')
    for ann in anns:
        vds = annotate_genesets_dict(vds, geneset_dict, gene_key = 'va.ann.{}.gene_id'.format(ann), ann_label = ann)
    return(vds)

def geneset_annotated_genotypes_table(vds, freq_expr = '(isMissing(va.AC_filt_all) || va.AC_filt_all <= 1) && (isMissing(va.AC_internal) || va.AC_internal == 1)'):
    '''Export per-sample gene set genotypes table for aggregation.
        
    Parameters
    ----------
    vds : {[Hail.VariantDataset]}
        Variantdataset
    freq_expr : {str}, optional
        Hail expression language string for filtering variants (the default is '(isMissing(va.AC_filt_all) || va.AC_filt_all <, which [default_description])
    '''
    logger.info('Generate KeyTable of singletons.')
    gkt = (
        vds
            .filter_variants_expr(freq_expr)
            .annotate_variants_expr('va = select(va, s, csq, gt, genesets)') # , expression_proportions
            .variants_table()
            .flatten()
    )
    return(gkt)

def geneset_counts_table(
    gkt, 
    anns = [ 'all', 'canonical', 'basic',  'cseven', 'pfour', 'pthree', 'pseven', 'primate' ],  # 'all_expr', 
    genesets = [ 'gnomAD_constraint_category__pLI', 'gnomAD_constraint_category__pNull', 'gnomAD_constraint_category__All' ]
):
    '''
    Export gene set counts table from gene set genotype table.
    '''
    logger.info('Reduce to 100 partitions.')
    gkt = gkt.repartition(100, shuffle = False)

    logger.info('Define gene set loop.')
    exprs = [ '`va.genesets.{}.{}`'.format(ann, g) for ann in anns for g in genesets  ]
    
    logger.info('Generate non-ref counts for gene sets.')
    kt_list = []
    for expr in exprs:
        ann = re.search('va.genesets\.([a-z_0-9]+)\.', expr).group(1)
        name = re.search('va.genesets\.[a-z_0-9]+\.(.*)\`', expr).group(1)
        #logger.info('Generate non-ref counts for gene set: %s with analysis consequences: %s', name, ann)
        kt_list.append(
            gkt
                .filter(expr)
                .aggregate_by_key(
                    [ 's = `va.s`', 'csq = `va.csq.{}`'.format(ann) ],
                    'X = `va.gt`.sum().toInt()'
                )
                .annotate('name = "{}", ann = "{}"'.format(name, ann))
        )
    logger.info('Join results keytables and repartition to 500.')
    kt = union_keytables(kt_list).repartition(500, shuffle = False)
    return(kt)