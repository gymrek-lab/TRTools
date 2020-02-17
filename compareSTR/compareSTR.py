#!/usr/bin/env python3
"""
Tool for comparing two STR callsets
"""

# Load external libraries
import argparse
import os
import numpy as np
import sys
import vcf

# Load local libraries
if __name__ == "compareSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import common
    import mergeutils
    import tr_harmonizer as trh
    import utils
else: # pragma: no cover
    import strtools.utils.common as common  # pragma: no cover
    import strtools.utils.mergeutils as mergeutils  # pragma: no cover
    import strtools.utils.tr_harmonizer as trh # pragma: no cover
    import strtools.utils.utils as utils  # pragma: no cover

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf1", help="First VCF file to compare (must be sorted, bgzipped, and indexed)", type=str, required=True)
    req_group.add_argument("--vcf2", help="First VCF file to compare (must be sorted, bgzipped, and indexed)", type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    req_group.add_argument("--vcftype", help="--vcf1 and --vcf2 must be of the same type. Options=%s"%trh.VCFTYPES.__members__, type=str, default="auto")
    ### Options for filtering input ###
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--region", help="Restrict to this region chrom:start-end", type=str)
    ### Stratify results ###
    stats_group = parser.add_argument_group("Metrics to stratify results")
    stats_group.add_argument("--stratify-fields", help="Comma-separated list of FORMAT fields to stratify by", type=str)
    stats_group.add_argument("--stratify-binsizes", help="Comma-separated list of binsizes to stratify each field on. Must be same length as --stratify-fields.", type=str)
    stats_group.add_argument("--stratify-file", help="Set to 1 to stratify based on --vcf1. Set to 2 to stratify based on --vcf2. Set to 0 to apply stratification to both --vcf1 and --vcf2", default=0, type=int)
    args = parser.parse_args()
    return ags

def main(args):
    ### Check and load VCF files ###
    vcfreaders = mergeutils.LoadReaders([args.vcf1, args.vcf2])
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    ### Check inferred type of each is the same
    vcftype = mergeutils.GetVCFType(vcfreaders, args.vcftype)

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = GetMinRecords(current_records, chroms)
    done = DoneReading(current_records)

    while not done:
        if args.verbose: mergeutils.PrintCurrentRecords(current_records, is_min)
        if mergeutils.CheckMin(is_min): return 1
        # TODO, take info in for each of them
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        is_min = mergeutils.GetMinRecords(current_records, chroms)
        done = mergeutils.DoneReading(current_records)
    return 0 

if __name__ == "__main__":  # pragma: no cover
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)

