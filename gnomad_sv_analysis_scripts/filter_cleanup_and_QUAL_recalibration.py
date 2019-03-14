#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2018 Ryan Collins <rlcollins@g.harvard.edu>
# Distributed under terms of the MIT license.

"""
Apply final FILTER cleanup and QUAL score recalibration
"""


import argparse
import sys
import pysam
import csv
from numpy import median


#Define global variables
filts_to_remove = 'HIGH_PCRPLUS_NOCALL_RATE HIGH_PCRMINUS_NOCALL_RATE HIGH_SR_BACKGROUND PESR_GT_OVERDISPERSION'.split(' ')
NULL_GTs = [(0, 0), (None, None), (0, ), (None, ), (None, 2)]
HET_GTs = [(0, 1), (None, 1), (None, 3)]


def get_call_rate(record, samples):
    total_s = [s for s in record.samples if s in samples]
    total = len(total_s)
    nocall_s = [s for s in total_s if record.samples[s]['GT'] in NULL_GTs]
    nocall = len(nocall_s)
    callrate = 1 - ( nocall / total )
    return callrate


def recal_qual_score(record):
    quals = []
    for s in [s for s in record.samples]:
        GT = record.samples[s]['GT']
        if GT in NULL_GTs:
            continue
        elif GT in HET_GTs:
            quals.append(record.samples[s]['GQ'])
        else:
            quals.append(999)

    if len(quals) > 0:
        return int(median(quals))


def cleanup_vcf(vcf, fout, plus_samples, min_callrate=0.9):
    
    minus_samples = [s for s in vcf.header.samples if s not in plus_samples]

    for record in vcf:
        #Mark small (350bp-1kb) deletions with high PCRMINUS nocall rate
        # as PREDICTED_GENOTYPING_ARTIFACT
        if record.info['SVTYPE'] == 'DEL' \
        and record.info['SVLEN'] < 1000 \
        and record.info['SVLEN'] > 350:
            if 'HIGH_PCRMINUS_NOCALL_RATE' in record.filter:
                record.filter.add('PREDICTED_GENOTYPING_ARTIFACT')

        #Move PESR_GT_OVERDISPERSION from FILTER to INFO
        if 'PESR_GT_OVERDISPERSION' in record.filter:
            record.info['PESR_GT_OVERDISPERSION'] = True

        #Remove all HIGH_NOCALL_RATE and HIGH_SR_BACKGROUND tags from FILTER column
        newfilts = [filt for filt in record.filter if filt not in filts_to_remove]
        record.filter.clear()
        for filt in newfilts:
            record.filter.add(filt)
        if len(record.filter) == 0:
            record.filter.add('PASS')

        # #Mark sites with low PCR+ call rate
        # plus_callrate = get_call_rate(record, plus_samples)
        # if plus_callrate < min_callrate:
        #     if 'LOW_PCRPLUS_CALL_RATE' not in record.info.keys():
        #         record.info.keys().append('LOW_PCRPLUS_CALL_RATE')
        #     record.info['LOW_PCRPLUS_CALL_RATE'] = True

        # #Mark sites with low PCR- call rate
        # minus_callrate = get_call_rate(record, minus_samples)
        # if minus_callrate < min_callrate:
        #     if 'LOW_PCRMINUS_CALL_RATE' not in record.info.keys():
        #         record.info.keys().append('LOW_PCRMINUS_CALL_RATE')
        #     record.info['LOW_PCRMINUS_CALL_RATE'] = True

        #Recalibrate QUAL score
        newQUAL = recal_qual_score(record)
        if newQUAL is not None:
            record.qual = newQUAL

        #Only write out non-empty variants to output file
        for s in record.samples:
            if record.samples[s]['GT'] not in NULL_GTs:
                fout.write(record)
                break


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf (supports "stdin").')
    parser.add_argument('PCRPLUS_samples', help='List of PCRPLUS sample IDs.')
    parser.add_argument('fout', help='Output file (supports "stdout").')

    args = parser.parse_args()

    #Open connection to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin) 
    else:
        vcf = pysam.VariantFile(args.vcf)

    #Add new FILTER lines to VCF header
    NEW_FILTERS = ['##FILTER=<ID=PREDICTED_GENOTYPING_ARTIFACT,Description="Site ' + 
                   'is predicted to be a genotyping false-positive based on ' + 
                   'analysis of minimum GQs prior to GQ filtering.">']
    header = vcf.header
    for filt in NEW_FILTERS:
        header.add_line(filt)

    #Remove unused FILTER lines from VCF header
    for filt in filts_to_remove:
        header.filters.remove_header(filt)

    #Add new INFO lines to VCF header
    NEW_INFOS = ['##INFO=<ID=PESR_GT_OVERDISPERSION,Number=0,Type=Flag,Description=' + 
                 '"PESR genotyping data is overdispersed. Flags sites where genotypes' + 
                 ' are likely noisier.">']
    for info in NEW_INFOS:
        header.add_line(info)

    #Read list of PCR+ samples
    f_plus_samples = open(args.PCRPLUS_samples, 'r')
    plus_samples = f_plus_samples.read().splitlines()
    f_plus_samples.close()

    #Open connection to output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    #Cleanup VCF
    cleanup_vcf(vcf, fout, plus_samples)

    fout.close()


if __name__ == '__main__':
    main()
