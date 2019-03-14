#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2017 Matthew Stone <mstone5@mgh.harvard.edu>
# Distributed under terms of the MIT license.

"""
Convert CNV records in clustered VCF to RdTest bed format.
"""

import argparse
import sys
from collections import deque
from pysam import VariantFile


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input VCF')
    parser.add_argument('bed', help='Output BED', type=argparse.FileType('w'),
                        default=sys.stdout, nargs='?')
    args = parser.parse_args()

    # Prep bed
    header = '#chrom\tstart\tend\tname\tsamples\tsvtype\n'
    args.bed.write(header)
    entry = '{chrom}\t{start}\t{end}\t{name}\t{samples}\t{svtype}\n'

    vcf = VariantFile(args.vcf)
    for record in vcf:
        # Skip non-CNV
        if record.info['SVTYPE'] not in 'DEL DUP'.split():
            continue

        # Get bed interval and metadata
        chrom = record.chrom
        start = record.pos
        end = record.stop
        name = record.id
        svtype = record.info['SVTYPE']

        # Get list of called samples
        samples = deque()
        null_GTs = [(0, 0), (None, None), (0, ), (None, )]
        for sample in record.samples:
            gt = record.samples[sample]['GT']
            if gt not in null_GTs:
                samples.append(sample)
        if len(samples) == 0:
            continue
        samples = ','.join(sorted(set(samples)))

        args.bed.write(entry.format(**locals()))


if __name__ == '__main__':
    main()
