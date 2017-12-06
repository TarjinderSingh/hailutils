#!/usr/bin/python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *

from cloudio import *
from vds import *
from keytable import *

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
            parse_vep_transcript_consequences(
                annotate_vep_global(vds)
            )
         )
    )
    return(vds)
    
def run_vep_pipeline2(vds):
    logger.info('Annotating with VEP.')
    vds = vds.vep(config = '/vep/vep-gcloud.properties')

    logger.info('Parsing VEP output.')
    vds = (
        parse_canonical_transcript_consequences(
            parse_vep_transcript_consequences(
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

def parse_vep_transcript_consequences(vds):
    logger.info('Parse VEP transcript consequences.')
    return(
        vds
            .annotate_variants_expr(
                """
                va.tcsq = 
                    va.vep.transcript_consequences
                        .filter(tc => tc.biotype == "protein_coding") 
                        .map(
                            x => {
                                transcript_id: x.transcript_id,
                                canonical: x.canonical == 1,
                                consequence_term: x.consequence_terms.mkString(","),
                                minscore: x.consequence_terms.map(x => global.somap[x]).min(),
                                csq: let minscore = x.consequence_terms.map(x => global.somap[x]).min() in 
                                    if ((minscore <= 5) && (minscore != 3))
                                        "lof"
                                    else if (minscore == 3)
                                        "startlost"
                                    else if (minscore == 6)
                                        "stoplost"
                                    else if (minscore == 8 || minscore == 9)
                                        "ns"
                                    else if (minscore == 10)
                                        "mis"
                                    else if (minscore == 12)
                                        "splice"
                                    else if (minscore == 15)
                                        "syn"
                                    else
                                        NA: String,
                                gene_id: x.gene_id,
                                aachange: x.hgvsp,
                                lof: x.lof,
                                lof_filter: x.lof_filter,
                                lof_flags: x.lof_flags,
                                polyphen: x.polyphen_prediction,
                                sift: x.sift_prediction,
                                loftee: 
                                    if ((x.lof == "HC") && (x.lof_flags == ""))
                                        "HC"
                                    else if ((x.lof == "HC") && (x.lof_flags != ""))
                                        "HCflagged"
                                    else if (x.lof == "LC")
                                        "LC"
                                    else 
                                        NA: String
                            }
                        )
                """
            )
    )

def parse_canonical_transcript_consequences(vds):
    logger.info(
        'Determine minscores of canonical transcripts and all coding transcripts '
        'and identify all corresponding transcripts with those minscore values.'
    )
               
    vds = (
        vds
            .annotate_variants_expr(
                '''
                va.ann.canonical.minscore = va.tcsq.filter(tc => tc.canonical).map(tc => tc.minscore).min(),
                va.ann.all.minscore = va.tcsq.map(tc => tc.minscore).min()
                '''
            )
            .annotate_variants_expr(
                '''
                va.ann.canonical.min_tcsq = va.tcsq.filter(tc => tc.canonical && tc.minscore == va.ann.canonical.minscore),
                va.ann.all.min_tcsq = va.tcsq.filter(tc => tc.minscore == va.ann.all.minscore)
                '''
            )
    )
    
    logger.info('Annotate the genes and predicted consequences of all minscore transcripts.')
    exprs = [
        '''
        {0}.gene_id = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.gene_id).toSet()
            else
                NA: Set[String]
        ''',
        '''
        {0}.transcript_id = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.transcript_id).toSet()
            else
                NA: Set[String]
        ''',
        '''
        {0}.csq = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.csq).toSet().head()
            else
                NA: String
        ''',
        '''
        {0}.loftee = 
            if (isDefined({0}.minscore))
                let loftee = {0}.min_tcsq.map(tc => tc.loftee).toSet() in
                if (loftee.contains("HC"))
                    "HC"
                else if (loftee.contains("HCflagged"))
                    "HCflagged"
                else if (loftee.contains("LC"))
                    "LC"
                else 
                    NA: String
            else
                NA: String
        ''',
        '''
        {0}.polyphen = 
            if (isDefined({0}.minscore))
                let loftee = {0}.min_tcsq.map(tc => tc.polyphen).toSet() in
                if (loftee.contains("probably_damaging"))
                    "D"
                else if (loftee.contains("possibly_damaging"))
                    "P"
                else if (loftee.contains("benign"))
                    "B"
                else 
                    NA: String
            else
                NA: String
        ''',
        '''
        {0}.aachange = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.aachange).toSet()
            else
                NA: Set[String]
        '''
    ] 
    vds = vds.annotate_variants_expr(
        [ expr.format(ann) for expr in exprs for ann in [ 'va.ann.canonical', 'va.ann.all' ] ]
    )
    
    logger.info('Annotate variant with all transcripts with lowest minscore values.')
    vds = vds.annotate_variants_expr(
        '''
        va.ann.all.transcript_ids = 
            let gene2transcripts_structlist = 
                va.ann.all.gene_id.toArray()
                    .map(
                        x => {
                            gene_id: x, 
                            transcript_ids: va.ann.all.min_tcsq
                                .filter(tc => tc.gene_id == x)
                                .map(tc => tc.transcript_id)
                                .toSet()
                                .toArray()
                             }
                    ) in 
                index(gene2transcripts_structlist, gene_id)
                    .mapValues(x => x.transcript_ids)
        '''
    )
    return(vds)

def parse_selected_transcript_consequences(vds, transcript_ids, name):
    logger.info(
        'Determine minscores of selected transcripts '
        'and identify all corresponding transcripts with those values.'
    )
    
    vds = vds.annotate_global('global.{}_transcript_ids'.format(name), transcript_ids, TSet(TString()))
    
    vds = (
        vds
            .annotate_variants_expr(
                '''
                va.ann.{0}.minscore = 
                    va.tcsq
                        .filter(tc => global.{0}_transcript_ids.contains(tc.transcript_id))
                        .map(tc => tc.minscore)
                        .min()
                '''.format(name)
            )
            .annotate_variants_expr(
                '''
                va.ann.{0}.min_tcsq = 
                    va.tcsq.filter(
                        tc => global.{0}_transcript_ids.contains(tc.transcript_id) && 
                        tc.minscore == va.ann.{0}.minscore
                    )
                '''.format(name)
            )
    )
    
    logger.info('Annotate the genes and predicted consequences of all minscore transcripts.')
    exprs = [
        '''
        {0}.gene_id = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.gene_id).toSet()
            else
                NA: Set[String]
        ''',
        '''
        {0}.transcript_id = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.transcript_id).toSet()
            else
                NA: Set[String]
        ''',
        '''
        {0}.csq = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.csq).toSet().head()
            else
                NA: String
        ''',
        '''
        {0}.loftee = 
            if (isDefined({0}.minscore))
                let loftee = {0}.min_tcsq.map(tc => tc.loftee).toSet() in
                if (loftee.contains("HC"))
                    "HC"
                else if (loftee.contains("HCflagged"))
                    "HCflagged"
                else if (loftee.contains("LC"))
                    "LC"
                else 
                    NA: String
            else
                NA: String
        ''',
        '''
        {0}.polyphen = 
            if (isDefined({0}.minscore))
                let loftee = {0}.min_tcsq.map(tc => tc.polyphen).toSet() in
                if (loftee.contains("probably_damaging"))
                    "D"
                else if (loftee.contains("possibly_damaging"))
                    "P"
                else if (loftee.contains("benign"))
                    "B"
                else 
                    NA: String
            else
                NA: String
        ''',
        '''
        {0}.aachange = 
            if (isDefined({0}.minscore))
                {0}.min_tcsq.map(tc => tc.aachange).toSet()
            else
                NA: Set[String]
        '''
    ] 
    vds = vds.annotate_variants_expr([ expr.format('va.ann.{}'.format(name)) for expr in exprs  ])
    return(vds)

def parse_basic_transcript_consequences(vds):
    logger.info('Read GENCODE basic transcripts.')
    basic_transcript_ids = set(parse('gs://exome-qc/resources/gencode_v19/2017-10-10_gencode-basic-transcripts.tsv', split = '\t', skip = 1, cols = 1))
    logger.info('Parse consequences of GENCODE basic transcripts.')
    vds = parse_selected_transcript_consequences(vds, basic_transcript_ids, 'basic')
    return(vds)

def get_transcript_consequence_dataframe(vds, columns = [ 'v', 'va.tcsq' ]):
    return(
        vds
            .variants_table()
            .flatten()
            .select(columns)
            .explode('va.tcsq')
            .annotate('v = str(v)')
            .flatten()
            .to_pandas()
    )



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
            .annotate_variants_expr('va = select(va, gnomad, nonpsych_gnomad)')
     )
    vds = vds.annotate_variants_vds(gnomad_vds, expr = 'va = merge(va, vds)')
    return(vds)
    
def annotate_nonpsych_exac_frequencies(vds):
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
    logger.info('Annotate with CADD 1.3 scores.')
    cadd_kt = (
        vds.hc
            .read_table('gs://exome-qc/resources/cadd/cadd1.3_whole_exome_SNVs.kt')
            .key_by('variant')
            .select([ 'variant', 'phred', 'rawscore'])
    )
    vds = vds.annotate_variants_table(cadd_kt, root = 'va.cadd13')
    return(vds)

def annotate_cadd10(vds):
    logger.info('Annotate with CADD 1.0 scores.')
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
                    'pLI',
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

    #logger.info('Annotate with loss-of-function intolerant genes.')
    #lof_intolerant_genes = set(mis_regions.filter('pLI >= 0.9').query('gene_id.collect().toSet()'))
    #vds = vds.annotate_global('global.lof_intolerant_genes', lof_intolerant_genes, TSet(TString()))
    
    logger.info('Annotate with missense intolerant genes.')
    mis_intolerant_genes = set(mis_regions.filter('mis_z >= 3.09').query('gene_id.collect().toSet()'))
    vds = vds.annotate_global('global.mis_intolerant_genes', mis_intolerant_genes, TSet(TString()))
    
    logger.info('Annotate with genes that contain      regions.')
    regional_intolerant_genes = set(mis_regions.filter('obs_exp <= 0.6').query('gene_id.collect().toSet()'))
    vds = vds.annotate_global('global.regional_intolerant_genes', regional_intolerant_genes, TSet(TString()))

    return(vds)

def annotate_nonpsych_lof_intolerant_genes(vds):
    gene_set = (
        vds.hc
            .read_table('gs://exome-qc/resources/exac_release0.3.1/constraint/exac_constraint_nonpsych.kt')
            .filter('pLI >= 0.9')
            .query('gene_id.collect().toSet()')
     )
    
    vds = (
        vds
            .annotate_global('global.lof_intolerant_genes', gene_set, TSet(TString()))
            .annotate_variants_expr(
                '''
                va.genesets.in_lof_intolerant_genes = 
                    if (! va.gene_id.filter(x => global.lof_intolerant_genes.contains(x)).isEmpty())
                        true
                    else
                        false
                '''
            )
    )
    return(vds)

def annotate_geneset(vds, name, geneset = None):
    if isinstance(geneset, set):
        vds = vds.annotate_global('global.{}'.format(name), gene_set, TSet(TString()))
    vds = (
        vds
            .annotate_variants_expr(
                '''
                va.genesets.in_{0} = 
                    if (! va.gene_id.filter(x => global.{0}.contains(x)).isEmpty())
                        true
                    else
                        false
                '''.format(name)
        )
    )
    return(vds)

def get_gencode_keytable():
    return(KeyTable.import_bed('gs://exome-qc/resources/gencode_v19/gencode.v19.cds.merged_by_exonid.merged.bed'))

def get_gencode8_keytable():
    return(KeyTable.import_bed('gs://exome-qc/resources/gencode_v19/gencode.v19.cds.merged_by_exonid.merged_p8.bed'))

def get_mask_keytable():
    return(
        KeyTable.import_bed(
            'gs://sczmeta_exomes/data/coverage/release_v1.5/high_confidence_regions_gencode_annotated_10x80_20170910.bed')
    )

def get_maskext_keytable():
    return(
        KeyTable.import_bed(
            'gs://sczmeta_exomes/data/coverage/release_v1.5/high_confidence_regions_extended_10x80_20170910.bed')
    )

def get_exome_called_keytable():
    return(
        KeyTable.import_bed('gs://sczmeta_exomes/data/regions/exome_calling_regions.merged.v1.bed')
    )

def annotate_paralog(vds):
    logger.info('Import paralog scores.')
    kt = (
        vds.hc
            .read_table('gs://exome-qc/resources/paralog_score/hg19.paralog.allScores.zscore.kt')
            .key_by('variant')
            .select(
                [
                    'variant',
                    'region', 
                    'Gene.refGene',
                    'ExonicFunc.refGene',
                    'PARASUBFAM_ID',
                    'PARASUBSCOREzSCORE'
                ]
            )
    )
    logger.info('Annotate VDS with missense regions.')
    vds = vds.annotate_variants_table(kt, 'va.paralog_score')
    return(vds)

def annotate_dbsnp(vds, version = 'hg19'):
    if version == 'hg19':
        path = 'gs://exome-qc/resources/dbsnp/human_9606_b150_GRCh37p13/dbsnp-150-grch37.vds'
    else:
        path = 'gs://exome-qc/grch38/dbsnp/human_9606_b150_GRCh38p7/dbsnp-150-grch38.vds'

    logger.info('Import dbsnp database (version %s).', version)
    avds = (
        vds.hc
            .read(path)
            .annotate_variants_expr('va.dbsnp = select(va.dbsnp, rsid, build)')
    )
    logger.info('Annotate variant with rsid.')
    vds = vds.annotate_variants_vds(avds, expr = 'va = merge(va, vds)')
    return(vds)

def run_current_pipeline(vds):
    vds = run_vep_pipeline2(vds)
    vds = parse_basic_transcript_consequences(vds)

    # Annotate with frequency databases
    vds = annotate_nonpsych_exac_frequencies(vds)
    vds = annotate_gnomad_frequencies(vds)

    # Annotate with CADD scores
    vds = annotate_cadd10(vds)
    vds = annotate_cadd13(vds)

    # Annotate with MPC, splice, and constrained regions
    vds = annotate_mpc(vds)
    vds = annotate_splice(vds)
    vds = annotate_constraint(vds)
    vds = annotate_paralog(vds)

    # dbSNP
    vds = annotate_dbsnp(vds)

    # Annotate with LoF intolerance
    # vds = annotate_nonpsych_lof_intolerant_genes(vds)
    return(vds)

def get_annotation_table(
    vds, 
    filter_expr = 'va.ann.all.minscore <= 15',
    keep_cols = [
        'v', 
        'va.ann.canonical.csq', 'va.ann.canonical.gene_id', 'va.ann.canonical.minscore',
        'va.ann.canonical.loftee', 'va.ann.canonical.polyphen', 
        'va.gnomad.AC', 'va.nonpsych_gnomad.AC'
    ],
    prettify_headers = ['va.ann.canonical', 'va.ann', 'va'],
    explode_gene = True):
    if filter_expr:
        logger.info('Filter annotation VDS based on expression: %s', filter_expr)
        vds = vds.filter_variants_expr(filter_expr)
    
    logger.info('Select columns and convert to KeyTable.')
    kt = (
        vds
            .variants_table()
            .flatten()
            .select(keep_cols)
    )
    kt = kt.rename(prettify_columns(kt.columns, prettify_headers))
    logger.info('%s variants observed in annotation table.', kt.count())
    
    if explode_gene:
        logger.info('Explode gene column.')
        kt = kt.explode('gene_id')
    return(kt)
    