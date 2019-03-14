#!/usr/bin/env bash

#############################
#    gnomAD SV Discovery    #
#############################

# Copyright (c) 2018 Harold Wang, Ryan L. Collins, and the Talkowski Lab
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>
# gnomAD credits: http://gnomad.broadinstitute.org/

#Wrapper to handle pre-filtering for vcf2baf_helper.py
#Collects BAF data for all samples present in a given VCF input


#####Usage statement
usage(){
cat <<EOF

usage: vcf2baf.sh [-h] [-z] input.vcf.gz output.txt

Wrapper script to generate B-allele frequencies for all heterozygous
sites present in each sample in an input VCF

Positional arguments:
  input.vcf.gz   path to bgzipped input VCF file
  output.txt     path to output BAF file

Optional arguments:
  -h  HELP      Show this help message and exit
  -z  BGZIP     bgzip & tabix index the BAF output file

EOF
}


#####Parse arguments
BGZIP=0
while getopts ":zh" opt; do
  case "$opt" in
    h)
      usage
      exit 0
      ;;
    z)
      BGZIP=1
      ;;
  esac
done
shift $(( ${OPTIND} - 1))
VCF=$1
OUTFILE=$2


#####Check for required input
if [ -z ${VCF} ] || [ -z ${OUTFILE} ]; then
  usage
  exit 0
fi


#Get path to vcf2baf bin
BIN=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )


#####Run BAF generation & write to OUTFILE
zcat ${VCF} | bcftools view -M2 -v snps - | grep -e '^#\|PASS' | \
python ${BIN}/vcf2baf_helper.py > ${OUTFILE}


#####Bgzip & tabix index OUTFILE, if optioned
if [ ${BGZIP} -gt 0 ]; then
  bgzip -f ${OUTFILE}
  tabix -s 1 -b 2 -e 2 -f ${OUTFILE}.gz
fi


