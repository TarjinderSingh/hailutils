#!/usr/bin/python

import os
import sys
import re
import logging
from pprint import pprint
from hail import *
import subprocess
try:
    from pyrunner.submit import *
except:
    from submit import *

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def mktemp(outdir, suffix):
    return(outdir + tempfile.mktemp(suffix))

def log_summary(vds):
    logger.info('%s samples and %s variants are observed in the data set', *list(vds.count()))
    
def log_variant_count(vds):
    logger.info('%s variants are observed in the data set.', vds.count_variants())
    
def extract_region(vds, locus, start = None, end = None):
    interval_str = ''
    if locus and not start and not end:
        interval_str = '{}'.format(str(locus))
    if locus and start and end:
        interval_str = '{}:{}-{}'.format(locus, str(start), str(end))
    if locus and start and not end:
        interval_str = '{}:{}-{}'.format(locus, str(start), str(start + 1))
    return(vds.filter_intervals(Interval.parse(interval_str)))
    
def extract_region_file(vds, region_file, keep = True, 
                        bed = False, interval_list = False, locus = False, variant = False,
                        header = None):
    ext = re.sub('\.bgz$', '', os.path.splitext(region_file)[1])
    if bed or ext == '.bed':
        kt = KeyTable.import_bed(region_file)
    elif interval_list or ext == '.interval_list':
        kt = KeyTable.import_interval_list(region_file)
    elif locus or variant:
        if header:
            kt = hc.import_table(region_file, impute = True)
        else:
            kt = hc.import_table(region_file, no_header = True, impute = True)
            if locus:
                kt = kt.rename(['Locus'])
            else:
                kt = kt.rename(['Variant'])
        kt = kt.key_by(kt.columns[0])
    return(vds.filter_variants_table(kt, keep = keep))
    
def extract_variants(vds, variant_str_list):
    variant_str_list = [ variant_str_list ] if isinstance(variant_str_list, str) else variant_str_list
    return(vds.filter_variants_list([ Variant.parse(vstr) for vstr in variant_str_list ], keep = True))

def extract_nonref_genotypes(vds, het_only = False, homvar_only = False):
    if het_only:
        expr = 'isHet'
    elif homvar_only:
        expr = 'isHomVar'
    else:
        expr = 'isCalledNonRef'
    return(vds.filter_samples_expr('gs.map(g => g.{}().toInt()).sum() > 0'.format(expr)))

def sites(vds):
    return(vds.drop_samples())

def subset_for_testing(vds, region = '22:15M-20M'):
    return(extract_region(vds, region))
    
def extract_inds(vds, sample_ids, keep = True):
    vds = vds.annotate_global('global.sample_ids', sample_ids, TSet(TString()))
    if keep:
        vds = vds.filter_samples_expr('global.sample_ids.contains(s)')
    else:
        vds = vds.filter_samples_expr('! global.sample_ids.contains(s)')
    return(vds)

def __change_contig(outfile, infile, add = True):
    f = open(outfile, 'w')
    p1 = subprocess.Popen([ 'zless', infile ], stdout = subprocess.PIPE, bufsize = 0)
    opt = '-a' if add else '-r'
    p2 = subprocess.Popen([ 'modify_vcf_contig.py', opt ], stdin = p1.stdout, stdout = subprocess.PIPE, bufsize = 0)
    p1.stdout.close()
    p3 = subprocess.Popen([ 'bgzip', '-c' ], stdin = p2.stdout, stdout = f, bufsize = 0)
    p2.stdout.close()

    output = p3.communicate()[0]
    print(output)

    #cmd = 'modify_vcf_contig.py'.format(infile)
    #cmd += ' -a' if add else ' -r'
    #cmd += ' | bgzip -c > {}'.format(outfile)
    #execute(cmd, silent = True, shell = True)

def change_contig(outfile, infile, add = True):
    make_dir(outfile)
    outroot = change_ext(outfile, '', '.vcf.gz')
    f = gzip.open(infile, 'rb')
    g = open("{0}.vcf.part".format(outroot), 'w')
    if add:
        for line in f:
            line = line.decode()
            if line[0:1] == "#":
                g.write(line)
            else:
                g.write('chr' + line)
    else:
        for line in f:
            line = line.decode()
            if line[0:1] == "#":
                g.write(line)
            else:
                g.write(line[3:])
    f.close()
    g.close()
    rename("{0}.vcf".format(outroot))
    unix_cmd("bgzip -f {0}.vcf".format(outroot))


def export_vcf(
    outfile,
    vds, 
    annotate_info_exprs = [
        'va.info.AC = va.qc.AC',
        'va.info.AF = va.qc.AF',
        'va.info.AN = va.qc.nCalled * 2'
    ],
    filter_variants_expr = 'va.qc.AC >= 1',
    select_info_fields = [
        'AC', 'AF', 'AN',
        'BaseQRankSum', 'ClippingRankSum',
        'DP', 'FS', 'InbreedingCoeff',
        'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR', 
        'VQSLOD', 'culprit'
    ]
):
    
    logger.info('Annotate INFO field: %s', ', '.join(annotate_info_exprs))
    vds = vds.annotate_variants_expr(annotate_info_exprs)
    
    if filter_variants_expr:
        logger.info('Filter variants based on expression: %s', filter_variants_expr)
        vds = vds.filter_variants_expr(filter_variants_expr)
    
    if select_info_fields:
        logger.info('Retain the following field fields: %s', ', '.join(select_info_fields))
    vds = vds.annotate_variants_expr('va.info = select(va.info, {})'.format(', '.join(select_info_fields)))
 
    log_summary(vds)
    logger.info('Writing VDS as VCF: %s', outfile)
    vds.export_vcf(outfile)  

def unflatten_vds(vds):
    all_cols = [ f.name for f in vds.variant_schema.fields ] 
    cols = [ c for c in all_cols if c.startswith('va') ]
    add_cols_expr = [ '{0} = va.`{0}`'.format(c) for c in cols  ] 
    drop_cols = ', '.join([ '`{}`'.format(c) for c in cols ])
    drop_cols_expr = 'va = drop(va, {})'.format(drop_cols)
    return(vds.annotate_variants_expr(add_cols_expr).annotate_variants_expr(drop_cols_expr))

def write_read_vds(vds, tmpdir = 'gs://tsingh'):
    path = mktemp(tmpdir, '.vds')
    logger.info('Write temporary VDS to %s', path)
    vds.write(path)
    logger.info('Read VDS back into memory.')
    return(vds.hc.read(path))

def variants_table(vds, cols = None, flatten = True, variant_string = True):
    kt = vds.variants_table()
    if flatten:
        kt = kt.flatten()
    if variant_string:
        kt = kt.annotate('v = str(v)')
    if cols:
        kt = kt.select(cols)
    return(kt)

def get_ann_type(annotation, schema, root = 'va'):
    return get_ann_field(annotation, schema, root).typ

def get_ann_field(annotation, schema, root='va'):
    anns = flatten_struct(schema, root, leaf_only=False)
    if not annotation in anns:
        logger.error("%s missing from schema.", annotation)
        sys.exit(1)
    return anns[annotation]

def flatten_struct(struct, root ='', leaf_only = True):
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

def get_variant_schema_names(vds):
    return(sorted(flatten_struct(vds.variant_schema, root = 'va').keys()))
