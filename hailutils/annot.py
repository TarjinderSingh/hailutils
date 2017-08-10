#!/usr/bin/python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("annot")
logger.setLevel(logging.INFO)

def get_ann_type(annotation, schema, root='va'):
    return get_ann_field(annotation, schema, root).typ

def get_ann_field(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]

def flatten_struct(struct, root='', leaf_only=True):
    result = {}
    for f in struct.fields:
        path = '%s.%s' % (root, f.name)
        if isinstance(f.typ, TStruct):
            result.update(flatten_struct(f.typ, path))
            if not leaf_only:
                result[path] = f
        else:
            result[path] = f
    return result

def get_numbered_annotations(vds, root='va.info'):
    """
    Get all 1-, A-, G- numbered annotations from a VDS based on the Number va attributes.
    In addition returns arrays with no Number or Number=. va attribute separately
    :param vds: Input VDS
    :param root: Place to find annotations (defaults to va.info)
    :return: annotations, a_annotations, g_annotations, dot_annotations as list[Field]
    """
    a_annotations = []
    g_annotations = []
    dot_annotations = []
    annotations = []

    release_info = get_ann_type(root, vds.variant_schema)
    for field in release_info.fields:
        if isinstance(field.typ, TArray):
            if 'Number' in field.attributes:
                number = field.attributes['Number']
                if number == "A":
                    a_annotations.append(field)
                elif number == "G":
                    g_annotations.append(field)
                else:
                    dot_annotations.append(field)
        else:
            annotations.append(field)

    logger.info("Found the following fields:")
    logger.info("1-based annotations: " + ",".join([x.name for x in annotations]))
    logger.info("A-based annotations: " + ",".join([x.name for x in a_annotations]))
    logger.info("G-based annotations: " + ",".join([x.name for x in g_annotations]))
    logger.info("dot annotations: " + ",".join([x.name for x in dot_annotations]))

    return annotations, a_annotations, g_annotations, dot_annotations

def index_into_arrays(a_based_annotations=None, r_based_annotations=None, vep_root=None, drop_ref_ann = False):
    """
    Creates annotation expressions to get the correct values when splitting multi-allelics
    :param list of str a_based_annotations: A-based annotations
    :param list of str r_based_annotations: R-based annotations
    :param str vep_root: Root of the vep annotation
    :param bool drop_ref_ann: If set to True, then the reference value of R-based annotations is removed (effectively converting them in A-based annotations)
    :return: Annotation expressions
    :rtype: list of str
    """
    annotations = []
    if a_based_annotations:
        for ann in a_based_annotations:
            annotations.append('{0} = {0}[va.aIndex - 1]'.format(ann))
    if r_based_annotations:
        expr = '{0} = {0}[va.aIndex]' if drop_ref_ann else '{0} = [{0}[0], {0}[va.aIndex]]'
        for ann in r_based_annotations:
            annotations.append(expr.format(ann))
    if vep_root:
        sub_fields = ['transcript_consequences', 'intergenic_consequences', 'motif_feature_consequences', 'regulatory_feature_consequences']
        annotations.extend(['{0}.{1} = {0}.{1}.filter(x => x.allele_num == va.aIndex)'.format(vep_root, sub_field) for sub_field in sub_fields])

    return annotations


####
# Map VEP consequences to integer

SO = ['transcript_ablation', # HIGH 0
      'frameshift_variant', # HIGH 1
      'stop_gained', # HIGH 2
      'start_lost', # 'initiator_codon_variant', # HIGH 3
      'splice_acceptor_variant', # HIGH 4
      'splice_donor_variant', # HIGH 5
      'stop_lost', # LOW 6
      'transcript_amplification', #HIGH
      'inframe_insertion', # MODERATE, NS
      'inframe_deletion', # MODERATE, NS
      'missense_variant', # MODERATE, NS, 10
      'protein_altering_variant', # MODERATE
      'splice_region_variant', # LOW
      'incomplete_terminal_codon_variant', # LOW
      'stop_retained_variant', # LOW
      'synonymous_variant', # LOW 
      'coding_sequence_variant', #MODIFIER 
      'mature_miRNA_variant', 
      '5_prime_UTR_variant',
      '3_prime_UTR_variant',
      'non_coding_transcript_exon_variant',
      'intron_variant',
      'NMD_transcript_variant',
      'non_coding_transcript_variant',    
      'upstream_gene_variant',
      'downstream_gene_variant',
      'TFBS_ablation', # MODERATE
      'TFBS_amplification',
      'TF_binding_site_variant',
      'regulatory_region_variant',
      'regulatory_region_ablation',
      'regulatory_region_amplification',
      'feature_elongation',
      'regulatory_region_variant',
      'feature_truncation',
      'intergenic_variant']


def run_vep_pipeline(vds):
    logger.info('Annotating with VEP.')
    vds = vds.vep(config = '/vep/vep-gcloud.properties')

    logger.info('Parsing VEP output.')
    vds = (
        parse_loftee(
            parse_vep(
                annotate_vep_global(vds)
            )
         )
    )
    return(vds)
    
def annotate_vep_global(vds):  
    ####
    # Extract protein-coding Ensembl IDs
    gene_ids = (
        vds
            .query_variants(
                """
                variants
                    .map(
                        v => va.vep.transcript_consequences
                            .filter(tc => tc.biotype == "protein_coding")
                            .map(x => x.gene_id)
                    )
                    .collect()
                    .flatten()
                    .toSet()
                    .toArray()
                """
            )
    )

    vds = vds.annotate_global('global.gene_ids', gene_ids, TArray(TString()))

    ####
    # Extract possible biotypes 
    biotypes = (
        vds
            .query_variants(
                """
                variants
                    .map(
                        v => va.vep.transcript_consequences
                            .map(tc => tc.biotype)
                        )
                    .collect()
                    .flatten()
                    .toSet()
                    .toArray()
                """
            )
    )

    vds = vds.annotate_global('global.biotypes', biotypes, TArray(TString()))

    ####
    # Extract consequences for protein-coding genes
    csqs = (   
        vds
            .query_variants(
                """
                variants
                    .map(
                        v => va.vep.transcript_consequences
                            .filter(tc => tc.biotype == "protein_coding")
                            .map(tc => tc.consequence_terms)
                        )
                    .collect()
                    .flatten()
                    .flatten()
                    .toSet()
                    .toArray()
                """
            )
    )

    vds = vds.annotate_global('global.csqs', csqs, TArray(TString()))

    ####
    # Annotate with sequence ontology map (essential for parse_vep)
    so_dict = { s: i for i, s in enumerate(SO) }
    vds = vds.annotate_global('global.somap', so_dict, TDict(TString(), TInt()))
    
    return(vds)

def parse_vep(vds):
    vds = (
        vds
            # Identify only protein-coding transcripts
            .annotate_variants_expr(
                """
                va.tcsq = 
                    va.vep.transcript_consequences
                        .filter(tc => tc.biotype == "protein_coding")
                """
            )
            # Find worst consequence sequence ontology score across transcripts
            .annotate_variants_expr(
                """
                va.minscore =
                    va.tcsq
                        .map(
                            tc => tc.consequence_terms
                                .map(x => global.somap[x])
                            )
                        .flatten()
                        .min()
                """
            )
            .annotate_variants_expr(
                """
                va.csq = 
                    if (va.minscore <= 5)
                        "lof"
                    else if (va.minscore == 8 || va.minscore == 9)
                        "ns"
                    else if (va.minscore == 10)
                        "mis"
                    else if (va.minscore == 12)
                        "splice"
                    else if (va.minscore == 15)
                        "syn"
                    else
                        NA: String
                """
            )
            # Find worst transcript consequence struct(s)
            .annotate_variants_expr(
                """
                va.min_tcsq = 
                    va.tcsq
                        .filter(
                            tc => tc.consequence_terms
                                .map(x => global.somap[x]).min() == va.minscore
                        )
                """
            )
            # Find gene_ids for the worst consequence
            .annotate_variants_expr(
                """
                va.gene_id = 
                    va.min_tcsq
                        .map(tc => tc.gene_id)
                        .toSet()
                        .toArray()
                """
            )
            # Find transcripts affected by the worst functional consequence
            .annotate_variants_expr(
                """
                va.transcript_ids = 
                    let gene2transcripts_structlist = 
                        va.gene_id
                            .map(
                                x => {gene_id: x, 
                                      transcript_ids: va.min_tcsq
                                          .filter(tc => tc.gene_id == x)
                                          .map(tc => tc.transcript_id)
                                          .toSet()
                                          .toArray()
                                     }
                            ) in 
                        index(gene2transcripts_structlist, gene_id)
                            .mapValues(x => x.transcript_ids)
                """            
            )
            # Annotate with amino acid change
            .annotate_variants_expr(
                """
                va.aachange =       
                    let aachange_structlist =               
                        va.transcript_ids
                            .values
                            .flatten()
                            .map(x => {
                                transcript_id: x,
                                aachange: va.min_tcsq
                                    .filter(tc => tc.transcript_id == x)
                                    .map(tc => tc.hgvsp)
                                    .toSet()
                                    .toArray()
                                }
                            ) in                    
                        index(aachange_structlist, transcript_id)
                            .mapValues(x => x.aachange)  
                """      
            )
    )
    return(vds)

def parse_loftee(vds):
    vds = (
        vds
            # Annotate with genes affected by a HC loftee variant
            .annotate_variants_expr(
                """
                va.loftee = 
                    if (va.min_tcsq.exists(tc => tc.lof == "HC"))
                        va.min_tcsq
                            .filter(tc => tc.lof == "HC")                     
                            .map(
                                tc => tc.gene_id
                            )
                            .toSet()
                            .toArray()
                    else if (va.min_tcsq.exists(tc => tc.lof == "LC"))
                        [ "LC" ]
                    else 
                        NA: Array[String]
                """ 
            )
            # Annotate with a loftee flag, which indicates additional problems in addition to HC/LC
            .annotate_variants_expr(
                """
                va.loftee_flags = 
                    if (! va.loftee.isEmpty() || va.loftee != ["LC"])
                        let flags = 
                            va.loftee
                                .map(x => va.min_tcsq.filter(tc => tc.gene_id == x).map(tc => tc.lof_flags))
                                .flatten()
                                .toSet()
                        in 
                            if (flags.contains(""))
                                "HC"
                            else
                                flags.mkString(",")
                    else
                        NA: String    
                """
            )
    )
    return(vds)
           
def annotate_mpc(vds):
    logger.info('Annotate with MPC scores.')
    kt = (
        vds.hc.read_table(
            'gs://exome-qc/resources/missense_constraint/constraint_official/fordist_constraint_official_mpc_values.kt'
        )

        .key_by('variant')

        .select(
            [
                'variant',
                'ENSG',
                'gene_name',
                'SIFT',
                'PolyPhen',
                'obs_exp',
                'mis_badness',
                'fitted_score',
                'MPC'
            ]
        )
    )         
    vds = vds.annotate_variants_table(kt, root = 'va.mpc')
    return(vds)

def annotate_splice(vds):
    logger.info('Annotate with splice scores.')
    kt = (
        vds.hc.read_table(
            'gs://exome-qc/resources/splice_region_constraint/splice_constrained_variants.kt'
        )
    
        .key_by('region')

        .select(
            [
                'region',
                'chrom',
                'pos',
                'ref',
                'splice_position',
                'gene_id',
                'gene_name',
                'splice_region_damaging'
            ]
        )
    )
    vds = vds.annotate_variants_table(kt, root = 'va.splice')
    return(vds)
           
def annotate_gnomad_frequencies(vds):
    logger.info('Annotate with gnomAD frequencies.')
    gnomad_vds = (
        vds.hc
            .read('gs://sczmeta_exomes/data/gnomad_reference/release_v1.1/gnomad_merged.reduced.vep.r2.0.1.nonpsych.sites.vds')
            .annotate_variants_expr('va = select(va, gnomad)')
     )
    vds = vds.annotate_variants_vds(gnomad_vds, expr = 'va = merge(va, vds)')
    return(vds)
    
def annotate_nonpsychexac_frequencies(vds):
    logger.info('Annotate with non-psych ExAC frequencies.')
    exac_vds = vds.hc.read('gs://exome-qc/resources/exac_release0.3.1/ExAC.r0.3.nonpsych.sites.vds')
    vds = vds.annotate_variants_vds(exac_vds, 'va.in_nonpsych_ExAC = vds.in_nonpsych_ExAC')
    return(vds)
    
def annotate_discovEHR_frequencies(vds):
    logger.info('Annotate with non-psych ExAC frequencies.')
    discov_vds = vds.hc.read('gs://exome-qc/resources/discovEHR/discovEHR_freeze_50.vds')
    vds = vds.annotate_variants_vds(discov_vds, 'va.in_discovEHR = vds.in_discovEHR')
    return(vds)

def annotate_cadd13(vds):
    logger.info('Annotate with CADD scores.')
    cadd_kt = (
        vds.hc
            .read_table('gs://exome-qc/resources/cadd/cadd1.3_whole_exome_SNVs.kt')
            .key_by('variant')
            .select([ 'variant', 'phred', 'rawscore'])
    )
    vds = vds.annotate_variants_table(cadd_kt, root = 'va.cadd13')
    return(vds)

def annotate_cadd10(vds):
    logger.info('Annotate with CADD scores.')
    cadd_kt = (
        vds.hc
            .read_table('gs://exome-qc/resources/cadd/cadd1.0_whole_exome_SNVs.kt')
            .key_by('variant')
            .select([ 'variant', 'phred', 'rawscore'])
    )
    vds = vds.annotate_variants_table(cadd_kt, root = 'va.cadd10')
    return(vds)

def annotate_constraint(vds):
    logger.info('Import missense constrained regions.')
    mis_regions = (
        vds.hc
            .read_table('gs://exome-qc/resources/missense_constraint/constraint_official/missense_constrained_subregions.kt')
            .key_by('region')
            .select(
                [
                    'region',
                    'region_name',
                    'gene_id',
                    'gene_name',
                    'obs_exp',
                    'chisq_diff_null',
                    'lambda_constrained',
                    'multiregional',
                    'mis_z',
                    'n_regions'
                ]
            )
    )
    
    logger.info('Annotate VDS with missense regions.')
    vds = vds.annotate_variants_table(
        mis_regions.select(
            [
                'region',
                'gene_id', 
                'gene_name',
                'obs_exp',
                'lambda_constrained'
            ]
        ), 
        root = 'va.constraint'
    )

    logger.info('Annotate with loss-of-function intolerant genes.')
    lof_intolerant_genes = set(mis_regions.filter('pLI >= 0.9').query('gene_id.collect().toSet()'))
    vds = vds.annotate_global('global.lof_intolerant_genes', lof_intolerant_genes, TSet(TString()))
    
    logger.info('Annotate with missense intolerant genes.')
    mis_intolerant_genes = set(mis_regions.filter('mis_z >= 3.09').query('gene_id.collect().toSet()'))
    vds = vds.annotate_global('global.mis_intolerant_genes', mis_intolerant_genes, TSet(TString()))
    
    logger.info('Annotate with genes that contain constrained regions.')
    regional_intolerant_genes = set(mis_regions.filter('obs_exp <= 0.6').query('gene_id.collect().toSet()'))
    vds = vds.annotate_global('global.regional_intolerant_genes', regional_intolerant_genes, TSet(TString()))

    return(vds)