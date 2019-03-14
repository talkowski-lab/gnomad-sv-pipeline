#!/usr/bin/env python

# Copyright (c) 2018 Harold Wang, Ryan L. Collins, and the Talkowski Lab
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# gnomAD credits: http://gnomad.broadinstitute.org/

"""
Helper script for workflow to calculates B-allele frequencies
per sample from an input VCF file
"""

#Import libraries
import argparse
from collections import deque
import numpy as np
import pandas as pd
import pysam
import boto3
import sys

#Function to load an S3-hosted VCF
def load_s3vcf(bucket, vcf_path, index_filename=None):
    """
    Load an S3-hosted VCF.

    Parameters
    ----------
    bucket : str
        S3 bucket
    vcf_path : str
        S3 key

    Returns
    -------
    vcf : pysam.VariantFile
    """
    s3 = boto3.client('s3')
    url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={'Bucket': bucket, 'Key': vcf_path},
            ExpiresIn=86400)

    return pysam.VariantFile(url, index_filename=index_filename)

#Function to filter VCF records
def filter_records(record):
    """
    Filter VCF records to those informative for BAF genotyping.

    Returns only records which match all of the following criteria:
    1) Biallelic
    2) SNP
    3) FILTER == PASS

    Parameters
    ----------
    records : iterator of pysam.VariantRecords

    Returns
    ------
    record : pysam.VariantRecord
    """

    # for record in records:
    # Restrict to biallelic sites
    if len(record.alleles) > 2:
        return

    # Restrict to variants which PASS
    if record.filter.keys() != ['PASS']:
        return

    # Restrict to SNPs
    ref, alt = record.alleles
    if len(ref) > 1 or len(alt) > 1:
        return

    return record

#Function to calculate BAF per VCF record
def calc_BAF(record, samples=None):
    """

    Parameters
    ----------
    record : pysam.VariantRecord
    samples : list of str, optional
        Subset of samples in record to consider

    Returns
    -------
    bafs : np.ndarray of np.float
        BAF at site for each sample
    """

    def _is_het(sample):
        return record.samples[sample]['GT'] == (0, 1)

    def _calc_BAF(sample):
        if not _is_het(sample):
            return np.nan

        DP = record.samples[sample]['DP']
        AD = record.samples[sample]['AD']

        if DP!=None and DP > 10: # SNP sites with >10 DP are included in BAF profile
            return AD[0] / DP
        else:
            return np.nan

        

    if samples is None:
        samples = record.samples.keys()

    bafs = np.atleast_2d(np.array([_calc_BAF(sample) for sample in samples], dtype=np.float))

    return bafs

#Function to normalize BAF estimations
def normalize_bafs(bafs, max_std=0.2):
    """
    Normalize BAFs and exclude outlying sites
    Normalize so per variant median BAF==0.5. Ignore sites with more than 0.2 standard deviation across samples. 
    
    Parameters
    ----------
    bafs : np.ndarray (n_sites x n_samples)
    max_std : float, optional
        Maximium standard deviation permitted at a site

    Returns
    -------
    normalized_bafs : np.ndarray
    """

    # Convert to n_samples x n_sites
    bafs = bafs.transpose()

    # Remove variants not informative in any sample (all NA)
    bafs = bafs.loc[:, ~bafs.isnull().all()]  # .copy()

    # Remove sites with excessive variance
    # Permit sites with a single sample (SD=NA)
    std = bafs.std()
    bafs = bafs.loc[:, ((std < max_std) | std.isnull())]

    # Center each site's median BAF at 0.5
    bafs = bafs - bafs.median()
    bafs = bafs + 0.5

    return bafs
    
#Main function
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    # parser.add_argument('vcf', help='GATK VCF.')
    # parser.add_argument('--tbi', help='Path to VCF index. (Required for tabix '
                        # 'querying when accessing an S3-hosted VCF.)')
    # parser.add_argument('-w', '--window', type=int, default=None,
                        # help='Window outside variant start and end to query '
                        # 'for SNPs. Defaults to CNV length if not specified.')
    # parser.add_argument('-s', '--samples', type=str, default=None,
                    # help='Samples')
    # parser.add_argument('-i', '--ID', type=str, default='test',help='Samples')
    # parser.add_argument('-t', '--type', type=str, default='test',help='Samples')
    parser.add_argument('-b', '--batch', default='batch.txt')
                    # help='Samples')
    args = parser.parse_args()
    vcf = pysam.VariantFile(sys.stdin)
    while True:
        try:
            record=next(vcf)
            record=filter_records(record)
            if record:
                site=[record.pos]
                site=np.array(site, dtype=np.int)
                samples = list(vcf.header.samples)
                baf=calc_BAF(record)
                # print(baf.shape)
                baf = pd.DataFrame(baf)
                baf.columns = samples
                baf = baf.set_index(site)
                baf = normalize_bafs(baf)
                baf.index.name = 'sample'
                baf = baf.reset_index()
                bf = pd.melt(baf, id_vars=['sample'], var_name='pos', value_name='baf')
                bf = bf.loc[~bf.baf.isnull()]
                called_bafs = bf
                called_bafs['chrom'] = record.chrom
                called_bafs['pos'] = called_bafs.pos.astype(int)
                cols = 'chrom pos baf sample'.split()
                called_bafs = called_bafs[cols]
                if not called_bafs.empty:
                    called_bafs[cols].to_csv(sys.stdout, index=False, mode='a',header=False, sep='\t')
        except StopIteration:
            break

#Main block
if __name__ == '__main__':
    main()
