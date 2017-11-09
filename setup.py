#!/usr/bin/env python

from distutils.core import setup

setup(
    name='hailutils',
    version='1.0',
    description='Methods and scripts that compliment Hail.',
    author='Tarjinder Singh',
    author_email='tsingh@broadinstitute.org',
    url='https://github.com/TarjinderSingh/hailutils/',
    scripts=['bin/modify_vcf_contig.py', 'bin/bgzip', 'bin/tabix']
)
