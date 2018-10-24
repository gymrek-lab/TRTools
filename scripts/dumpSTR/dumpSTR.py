#!/usr/bin/env python

"""
Tool for filtering and QC of STR genotypes

Example command:

./dumpSTR.py \
--vcf /storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr22.allfilters.vcf.gz \
--out test \
--min-call-DP 10 \
--max-call-DP 1000 \
--min-call-Q 0.9 \
--max-call-flank-indel 0.15 \
--max-call-stutter 0.15 \
--min-locus-callrate 0.8 \
--min-locus-hwep 0.01 \
--min-locus-het 0 \
--max-locus-het 1 \
--use-length \
--filter-regions /storage/resources/dbase/human/hg19/hg19_segmentalduplications.bed \
--filter-regions-names SEGDUP \
--filter-hrun
"""

# Load external libraries
import argparse
import vcf
from vcf.parser import _Filter

# Load custom libraries
import filters

def BuildFilters(args):
    fdict = []
    if args.min_locus_callrate is not None:
        fdict.append(filters.Filter_MinLocusCallrate(args.min_locus_callrate))
    if args.min_locus_hwep is not None:
        fdict.append(filters.Filter_MinLocusHWEP(args.min_locus_hwep, args.use_length))
    if args.min_locus_het is not None:
        fdict.append(filters.Filter_MinLocusHet(args.min_locus_het, args.use_length))
    if args.max_locus_het is not None:
        fdict.append(filters.Filter_MaxLocusHet(args.max_locus_het, args.use_length))
    if args.filter_hrun is not None:
        fdict.append(filters.Filter_LocusHrun())
    return fdict

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)

    call_group = parser.add_argument_group("Call-level filters") # TODO add GangSTR filters
    call_group.add_argument("--min-call-DP", help="Minimum call coverage", type=int)
    call_group.add_argument("--max-call-DP", help="Maximum call coverage", type=int)
    call_group.add_argument("--min-call-Q", help="Minimum call quality score", type=float)
    call_group.add_argument("--max-call-flank-indel", help="Maximum call flank indel rate", type=float)
    call_group.add_argument("--max-call-stutter", help="Maximum call stutter rate", type=float)

    locus_group = parser.add_argument_group("Locus-level filters")
    locus_group.add_argument("--min-locus-callrate", help="Minimum locus call rate", type=float)
    locus_group.add_argument("--min-locus-hwep", help="Filter loci failing HWE at this p-value threshold", type=float)
    locus_group.add_argument("--min-locus-het", help="Minimum locus heterozygosity", type=float)
    locus_group.add_argument("--max-locus-het", help="Maximum locus heterozygosity", type=float)
    locus_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    locus_group.add_argument("--filter-regions", help="Comma-separated list of BED files of regions to filter", type=str)
    locus_group.add_argument("--filter-regions-names", help="Comma-separated list of filter names for each BED filter file", type=str)
    locus_group.add_argument("--filter-hrun", help="Filter STRs with long homopolymer runs.", action="store_true")

    args = parser.parse_args()

    # Load VCF file
    invcf = vcf.Reader(open(args.vcf, "rb"))

    # Set up filter list
    invcf.filters = {}
    filter_list = BuildFilters(args)
    for f in filter_list:
        short_doc = f.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        invcf.filters[f.filter_name()] = _Filter(f.filter_name(), short_doc)

    # Add new INFO fields - TODO

    # Set up output files
    outvcf = vcf.Writer(open(args.out + ".vcf", "w"), invcf)
    loclog = open(args.out + ".loclog.tab", "w")
    samplog = open(args.out + ".samplog.tab", "w")

if __name__ == "__main__":
    main()
