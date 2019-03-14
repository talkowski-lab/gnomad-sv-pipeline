#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Extract & classify trio allele counts for all complete sets of non-null genotypes
"""

import argparse
import sys
from collections import defaultdict
import pysam


GTs_to_skip = './. None/None 0/None None/0'.split()


def gather_info(vcf, fout, pro, fa, mo, no_header = False, mvr_out = None):
    sex_chroms = 'X Y chrX chrY'.split()
    CNVs = 'DEL DUP'.split()

    labels = 'INCOMPLETE APPARENT_DE_NOVO_HET APPARENT_DE_NOVO_HOM ' + \
             'MENDELIAN_CHILD_HET MENDELIAN_CHILD_HOM ' + \
             'UNTRANSMITTED_HET UNTRANSMITTED_HOM'
    labels = labels.split()
    svtypes = 'ALL DEL DUP INS INV CPX CTX BND'.split()

    #Prep dictionary for counting hets & MVRs by SVTYPE
    counts = {}
    for label in labels:
        for svtype in svtypes:
            for filt in 'PASS FAIL'.split():
                key = '_'.join([svtype, filt, label])
                if key not in counts.keys():
                    counts[key] = 0

    #Write header to output file(s)
    if not no_header:
        keys = list(counts.keys())
        header = '{0}\t{1}\n'.format('#PROBAND', '\t'.join(keys))
        fout.write(header)

        if mvr_out is not None:
            mvr_header = '#chr\tstart\tend\tVID\tSVTYPE\tPROBAND\n'
            mvr_out.write(mvr_header)

    trio_samples = [pro, fa, mo]

    for record in vcf:
        # #Do not include UNRESOLVED variants
        # if 'UNRESOLVED' in record.info.keys() \
        # or 'UNRESOLVED_TYPE' in record.info.keys() \
        # or 'UNRESOLVED' in record.filter:
        #     continue

        # #Only consider PASS variants
        # if 'PASS' not in record.filter:
        #     continue

        #Do not include variants from sex chromosomes
        if record.chrom in sex_chroms:
            continue

        #Do not include multiallelic variants
        if 'MULTIALLELIC' in record.info.keys() \
        or 'MULTIALLELIC' in record.filter \
        or len(record.alts) > 1:
            continue

        #Get GTs for trio
        GTs = [get_GT(record, ID) for ID in trio_samples]

        #Skip sites that are completely reference or null
        if len([g for g in GTs if g == '0/0' or g in GTs_to_skip]) == 3:
            continue

        #Convert to ACs
        ACs = [get_AC(g) for g in GTs]

        #Classify ACs
        label = classify_trio_AC(ACs, record)

        #Add counts to dict, as appropriate
        svtype = record.info['SVTYPE']
        filts = ','.join([f for f in record.filter])
        if filts == 'PASS':
            rfilt = 'PASS'
        else:
            rfilt = 'FAIL'
        counts['ALL_{0}_{1}'.format(rfilt, label)] += 1
        counts['{0}_{1}_{2}'.format(svtype, rfilt, label)] += 1

        #Write MVR coordinates to file, if optioned
        labels_to_skip_coords = 'INCOMPLETE MENDELIAN_CHILD_HET ' + \
                                'MENDELIAN_CHILD_HOM UNTRANSMITTED_HET'
        labels_to_skip_coords.split()
        if mvr_out is not None \
        and label not in labels_to_skip_coords:
            if record.info['SVTYPE'] == 'INS':
                mvr_record_vals = [record.chrom, record.pos, record.pos+1,
                                  record.id, record.info['SVTYPE'], pro, label]
            else:
                mvr_record_vals = [record.chrom, record.pos, record.stop,
                                  record.id, record.info['SVTYPE'], pro, label]
            mvr_newline = '\t'.join(str(i) for i in mvr_record_vals)
            mvr_out.write(mvr_newline + '\n')

    # Write het and MVR counts to file
    newline = '{0}\t{1}\n'.format(pro, '\t'.join(str(i) for i in counts.values()))
    fout.write(newline)


def get_GT(record, ID):
    GT = record.samples[ID]['GT']
    if GT is not None:
        GT_str = '/'.join([str(i) for i in GT])
    else:
        GT_str = 'None/None'
    return GT_str


def get_AC(GT):
    if GT == 'None/None':
        AC = 'NA'
    else:
        AC = sum([int(a) for a in GT.split('/')])
    return str(AC)


#Classify allele count for a single variant
def classify_trio_AC(ACs, record):
    n_missing = len([c for c in ACs if c == 'NA'])
    if n_missing > 0:
        return 'INCOMPLETE'

    ACs = [int(c) for c in ACs]
    n_homref = len([c for c in ACs if c == 0])
    pro_homref = ACs[0] == 0
    n_homref_parents = len([c for c in ACs[1:3] if c == 0])
    n_het = len([c for c in ACs if c == 1])
    pro_het = ACs[0] == 1
    n_het_parents = len([c for c in ACs[1:3] if c == 1])
    n_homalt = len([c for c in ACs if c == 2])
    pro_homalt = ACs[0] == 2
    n_homalt_parents = len([c for c in ACs[1:3] if c == 2])

    if pro_homalt and n_homref_parents > 0:
        return 'APPARENT_DE_NOVO_HOM'
    elif pro_het and n_homref_parents == 2:
        return 'APPARENT_DE_NOVO_HET'
    elif pro_homref and n_homalt_parents > 0:
        return 'UNTRANSMITTED_HOM'
    elif pro_homref and n_het_parents > 0:
        return 'UNTRANSMITTED_HET'
    elif pro_homalt and (n_het_parents + n_homalt_parents) == 2:
        return 'MENDELIAN_CHILD_HOM'
    elif pro_het and n_homref_parents < 2:
        return 'MENDELIAN_CHILD_HET'
    else:
        # import pdb; pdb.set_trace()
        error_message = 'Trio ACs {0} for site {1} do not fit expectations'.format(ACs, record.id)
        raise ValueError(error_message)


#Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('fout', help='Output stats file (supports "stdout").')
    parser.add_argument('pro', help='Proband sample ID.')
    parser.add_argument('fa', help='Father sample ID.')
    parser.add_argument('mo', help='Mother sample ID.')
    parser.add_argument('--coordinates', help='File to write out MVR coordinates.',
                        default = None, dest = 'mvr_out')
    parser.add_argument('--no-header', help='Do not write header line.',
                        action = 'store_true', default = False)

    args = parser.parse_args()

    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    if args.mvr_out is not None:
        mvrfout = open(args.mvr_out, 'w')
    else:
        mvrfout = None

    gather_info(vcf, fout, args.pro, args.fa, args.mo, 
                no_header = args.no_header, mvr_out = mvrfout)

    fout.close()

if __name__ == '__main__':
    main()
