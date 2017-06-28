#!/usr/bin/python

import os
import sys
import re
import logging
from hail import *
from pprint import pprint

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("annot")
logger.setLevel(logging.INFO)

def get_ann_type(annotation, schema, root='va'):
    return get_ann_field(annotation, schema, root).typ

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
