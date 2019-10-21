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
import pybedtools as pbt
from math import floor, ceil



GTs_to_skip = './. None/None 0/None None/0'.split()
sex_chroms = 'X Y chrX chrY'.split()



def get_GT(record, ID):
    """
    Extract a single sample's genotype from VCF record
    """
    GT = record.samples[ID]['GT']
    if GT is not None:
        GT_str = '/'.join([str(i) for i in GT])
    else:
        GT_str = 'None/None'
    return GT_str



def get_AC(GT):
    """
    Convert GT to AC for a single sample
    """
    if GT == 'None/None':
        AC = 'NA'
    else:
        AC = sum([int(a) for a in GT.split('/')])
    return str(AC)


def gather_parent_cnvs(vcf, fa, mo):
    """
    Create BEDTools corresponding to parent CNVs for converage-based inheritance
    """

    cnv_format = '{0}\t{1}\t{2}\t{3}\t{4}\n'

    fa_cnvs = ''
    mo_cnvs = ''

    for record in vcf:

        # Do not include variants from sex chromosomes
        if record.chrom in sex_chroms:
            continue

        # Process biallelic CNVs
        if record.info['SVTYPE'] in 'DEL DUP'.split() \
        and 'MULTIALLELIC' not in record.filter:

            # Father
            fa_ac = get_AC(get_GT(record, fa))
            if fa_ac != 'NA':
                if int(fa_ac) > 0:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                str(record.stop), 
                                                record.info['SVTYPE'], fa_ac)
                    fa_cnvs = fa_cnvs + new_cnv
            
            # Mother
            mo_ac = get_AC(get_GT(record, mo))
            if mo_ac != 'NA':
                if int(mo_ac) > 0:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                str(record.stop), 
                                                record.info['SVTYPE'], mo_ac)
                    mo_cnvs = mo_cnvs + new_cnv

        # Process multiallelic CNVs
        if record.info['SVTYPE'] == 'MCNV' and 'MULTIALLELIC' in record.filter:

            # Father
            fa_ac = get_GT(record, fa).split('/')[1]
            if fa_ac != 'None':
                fa_ac = int(fa_ac)
                if fa_ac < 2:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                    str(record.stop), 'DEL', 
                                                    str(2 - fa_ac))
                    fa_cnvs = fa_cnvs + new_cnv
                elif fa_ac > 2:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                    str(record.stop), 'DUP', 
                                                    str(fa_ac - 2))
                    fa_cnvs = fa_cnvs + new_cnv

            # Mother
            mo_ac = get_GT(record, mo).split('/')[1]
            if mo_ac != 'None':
                mo_ac = int(mo_ac)
                if mo_ac < 2:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                    str(record.stop), 'DEL', 
                                                    str(2 - mo_ac))
                    mo_cnvs = mo_cnvs + new_cnv
                elif mo_ac > 2:
                    new_cnv = cnv_format.format(record.chrom, str(record.pos), 
                                                    str(record.stop), 'DUP', 
                                                    str(mo_ac - 2))
                    mo_cnvs = mo_cnvs + new_cnv

    fa_cnvs = pbt.BedTool(fa_cnvs, from_string=True)
    mo_cnvs = pbt.BedTool(mo_cnvs, from_string=True)

    return fa_cnvs, mo_cnvs


def get_blacklist_hits(vcf, pro, blacklist = None):
    """
    Create BEDTools corresponding to parent CNVs for converage-based inheritance
    """

    breakpoint_format = '{0}\t{1}\t{2}\t{3}\n'

    breakpoints = ''

    for record in vcf:
        if 'MULTIALLELIC' not in record.filter:
            pro_ac = get_GT(record, pro).split('/')[1]
            if pro_ac != 'None' and pro_ac != 'NA':
                if int(pro_ac) > 0:
                    new_bp1 = breakpoint_format.format(record.chrom, str(record.pos),
                                                       str(record.pos + 1), record.id)
                    new_bp2 = breakpoint_format.format(record.chrom, str(record.stop),
                                                       str(record.stop + 1), record.id)
                    breakpoints = breakpoints + new_bp1 + new_bp2

    sv_bt = pbt.BedTool(breakpoints, from_string=True)
    
    if blacklist is not None:
        sv_bt_bl = sv_bt.intersect(blacklist, u=True, wa=True)
        bl_ids = list(set([f[3] for f in sv_bt_bl]))
    else:
        bl_ids = []
    
    return bl_ids


def classify_trio_AC(ACs, record):
    """
    Classify allele count for a single variant
    """
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
        error_message = 'Trio ACs {0} for site {1} do not fit expectations'.format(ACs, record.id)
        raise ValueError(error_message)



def reclassify_cnv_label(record, label, ACs, fa_cnvs, mo_cnvs, min_cov=0.5):
    """
    Attempt to reclassify inheritance labels for CNVs based on parent CNV coverage
    """

    bedline = '{0}\t{1}\t{2}\n'

    child_cnv = pbt.BedTool(bedline.format(record.chrom, 
                                           str(record.pos),
                                           str(record.stop)), 
                            from_string=True)
    cnvtype = record.info['SVTYPE']
    cnvlen = record.info['SVLEN']

    if cnvlen < 1000:
        fa_cov_dat = child_cnv.coverage(fa_cnvs.filter(lambda x: x[3] == cnvtype and x.length < 10000))
        mo_cov_dat = child_cnv.coverage(mo_cnvs.filter(lambda x: x[3] == cnvtype and x.length < 10000))
    else:
        fa_cov_dat = child_cnv.coverage(fa_cnvs.filter(lambda x: x[3] == cnvtype))
        mo_cov_dat = child_cnv.coverage(mo_cnvs.filter(lambda x: x[3] == cnvtype))

    fa_cov = float([f[6] for f in fa_cov_dat][0])
    mo_cov = float([f[6] for f in mo_cov_dat][0])

    if fa_cov >= min_cov or mo_cov >= min_cov:
        label = label.replace('APPARENT_DE_NOVO', 'MENDELIAN_CHILD')
    
    return label


# def _breakpoints_in_blacklist(record, blacklist=None):
#     """
#     Check whether variant breakpoints fall in blacklisted regions
#     """

#     if blacklist is None:
#         return False
#     else:
#         bedline = '{0}\t{1}\t{2}\n'
#         bp1 = pbt.BedTool(bedline.format(record.chrom, 
#                                            str(record.pos),
#                                            str(record.pos + 1)), 
#                             from_string=True)
#         bp2 = pbt.BedTool(bedline.format(record.chrom, 
#                                            str(record.stop),
#                                            str(record.stop + 1)), 
#                             from_string=True)

#         if len(bp1.intersect(blacklist)) > 0 \
#         or len(bp2.intersect(blacklist)) > 0:
#             return True
#         else:
#             return False


def gather_info(vcf, fout, pro, fa, mo, fa_cnvs, mo_cnvs, no_header = False, 
                mvr_out = None, qual_out = None, qual_step=50, 
                blacklisted_ids = None):
    CNVs = 'DEL DUP'.split()

    labels = 'INCOMPLETE APPARENT_DE_NOVO_HET APPARENT_DE_NOVO_HOM ' + \
             'MENDELIAN_CHILD_HET MENDELIAN_CHILD_HOM ' + \
             'UNTRANSMITTED_HET UNTRANSMITTED_HOM'
    labels = labels.split()
    qual_labels = 'APPARENT_DE_NOVO_HET MENDELIAN_CHILD_HET'.split()
    svtypes = 'ALL DEL DUP INS INV CPX CTX BND'.split()
    qual_bins = [i for i in range(0, 1000, qual_step)]

    #Prep dictionary for counting hets & MVRs by SVTYPE
    counts = {}
    for label in labels:
        for svtype in svtypes:
            for filt in 'PASS FAIL'.split():
                key = '_'.join([svtype, filt, label])
                if key not in counts.keys():
                    counts[key] = 0

    #Prep dictionary for counting child hets by QUAL score
    qual_counts = {}
    for label in qual_labels:
        for svtype in svtypes:
            for minQual in qual_bins:
                key = '_'.join([svtype, str(minQual), 
                                str(minQual + qual_step), label])
                if key not in counts.keys():
                    qual_counts[key] = 0

    #Write header to output file(s)
    if not no_header:
        keys = list(counts.keys())
        header = '{0}\t{1}\n'.format('#PROBAND', '\t'.join(keys))
        fout.write(header)

        if mvr_out is not None:
            mvr_header = '#chr\tstart\tend\tVID\tSVTYPE\tPROBAND\tINH\tFILTER\n'
            mvr_out.write(mvr_header)

        if qual_out is not None:
            qual_keys = qual_counts.keys()
            qual_header = '{0}\t{1}\n'.format('#PROBAND', '\t'.join(qual_keys))
            qual_out.write(qual_header)

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

        #Correct inheritance label for apparently de novo heterozygous CNVs
        if record.info['SVTYPE'] in CNVs \
        and label == 'APPARENT_DE_NOVO_HET':
            label = reclassify_cnv_label(record, label, ACs, fa_cnvs, mo_cnvs)

        #Add counts to dict, as appropriate
        svtype = record.info['SVTYPE']
        filts = ','.join([f for f in record.filter])
        if filts == 'PASS':
            rfilt = 'PASS'
        else:
            rfilt = 'FAIL'
        counts['ALL_{0}_{1}'.format(rfilt, label)] += 1
        counts['{0}_{1}_{2}'.format(svtype, rfilt, label)] += 1

        #Add counts to QUAL dict, as appropriate, only for PASS variants with SR support
        if label in qual_labels \
        and filts == 'PASS' \
        and 'SR' in record.info.get('EVIDENCE') \
        and record.id not in blacklisted_ids:
            qual = int(record.qual)
            qual_floor = str(qual_step * floor((qual - 1) / qual_step))
            qual_ceil = str(qual_step * ceil(qual / qual_step))
            qual_counts['ALL_{0}_{1}_{2}'.format(qual_floor, qual_ceil, label)] += 1
            qual_counts['{0}_{1}_{2}_{3}'.format(svtype, qual_floor, qual_ceil, label)] += 1

        #Write MVR coordinates to file, if optioned
        labels_to_skip_coords = 'INCOMPLETE MENDELIAN_CHILD_HET ' + \
                                'MENDELIAN_CHILD_HOM UNTRANSMITTED_HET'
        labels_to_skip_coords.split()
        if mvr_out is not None \
        and label not in labels_to_skip_coords:
            if record.info['SVTYPE'] == 'INS':
                mvr_record_vals = [record.chrom, record.pos, record.pos+1,
                                  record.id, record.info['SVTYPE'], pro, label, 
                                  rfilt]
            else:
                mvr_record_vals = [record.chrom, record.pos, record.stop,
                                  record.id, record.info['SVTYPE'], pro, label, 
                                  rfilt]
            mvr_newline = '\t'.join(str(i) for i in mvr_record_vals)
            mvr_out.write(mvr_newline + '\n')

    # Write het and MVR counts to file
    newline = '{0}\t{1}\n'.format(pro, '\t'.join(str(i) for i in counts.values()))
    fout.write(newline)

    # Write qual counts to file
    newline_qual = '{0}\t{1}\n'.format(pro, '\t'.join(str(i) for i in qual_counts.values()))
    qual_out.write(newline_qual)


#Main block
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf.')
    parser.add_argument('fout', help='Output stats file (supports "stdout").')
    parser.add_argument('pro', help='Proband sample ID.')
    parser.add_argument('fa', help='Father sample ID.')
    parser.add_argument('mo', help='Mother sample ID.')
    parser.add_argument('--coordinates', help='File to write out MVR coordinates.',
                        default = None, dest = 'mvr_out')
    parser.add_argument('--qual-out', help='File to write out de novo analysis ' +
                        ' stratified by QUAL score.',
                        default = None, dest = 'qual_out')
    parser.add_argument('--qual-bin-step', default=50, type=int, 
                        help='Size of QUAL score bins.', dest='qual_bin_size')
    parser.add_argument('--qual-blacklist', default=None, dest='blacklist',
                        help='Blacklist BED file to apply during QUAL score ' +
                        'stratified analysis of de novo rates.')
    parser.add_argument('--no-header', help='Do not write header line.',
                        action = 'store_true', default = False)

    args = parser.parse_args()

    vcf = pysam.VariantFile(args.vcf)

    if args.fout in '- stdout'.split():
        fout = sys.stdout
    else:
        fout = open(args.fout, 'w')

    if args.mvr_out is not None:
        mvrfout = open(args.mvr_out, 'w')
    else:
        mvrfout = None

    if args.qual_out is not None:
        qualfout = open(args.qual_out, 'w')
    else:
        qualfout = None

    fa_cnvs, mo_cnvs = gather_parent_cnvs(vcf, args.fa, args.mo)

    vcf = pysam.VariantFile(args.vcf)

    bl_hit_ids = get_blacklist_hits(vcf, args.pro, args.blacklist)

    vcf = pysam.VariantFile(args.vcf)    

    gather_info(vcf, fout, args.pro, args.fa, args.mo, fa_cnvs, mo_cnvs, 
                no_header = args.no_header, mvr_out = mvrfout, 
                qual_out = qualfout, qual_step = args.qual_bin_size,
                blacklisted_ids = bl_hit_ids)

    fout.close()



if __name__ == '__main__':
    main()


