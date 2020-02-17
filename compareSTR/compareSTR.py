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

def GetFormatFields(format_fields, format_binsizes, format_fileoption, vcfreaders):
    r"""Get which FORMAT fields to stratify on

    Also perform some checking on user arguments

    Parameters
    ----------
    format_fields : str
        Comma-separated list of FORMAT fields to stratify on
    format_binsizes: str
        Comma-separated list of min:max:binsize, one for each FORMAT field.
    format_fileoption : {0, 1, 2}
        Whether each format field needs to be in both readers (0), reader 1 (1) or reader 2 (2)
    vcfreaders : list of vcf.Reader
        List of readers. Needed to check if required FORMAT fields are present

    Returns
    -------
    formats : list of str
        List of FORMAT fields to stratify on
    binsizes : list of float
        List of binsizes for each FORMAT field
    """
    if format_fields is None or format_binsizes is None:
        return [], []
    formats = format_fields.split(",")
    binsizes = format_binsizes.split(",")
    if len(formats) != len(binsizes):
        raise ValueError("--stratify-formats must be same length as --stratify-binsizes")
    binsizes = [[float(x) for x in item.split(":")] for item in binsizes]
    for fmt in formats:
        check1 = (fmt in vcfreaders[0].formats)
        check2 = (fmt in vcfreaders[1].formats)
        if format_fileoption == 0 and not (check1 and check2):
            raise ValueError("FORMAT field %s must be present in both VCFs if --stratify-file=0"%fmt)
        if format_fileoption == 1 and not check1:
            raise ValueError("FORMAT field %s must be present in --vcf1 if --stratify-file=1"%fmt)
        if format_fileoption == 2 and not check2:
            raise ValueError("FORMAT field %s must be present in --vcf2 if --stratify-file=2"%fmt)
    return formats, binsizes

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
    stats_group.add_argument("--stratify-binsizes", help="Comma-separated list of min:max:binsize to stratify each field on. Must be same length as --stratify-fields.", type=str)
    stats_group.add_argument("--stratify-file", help="Set to 1 to stratify based on --vcf1. Set to 2 to stratify based on --vcf2. Set to 0 to apply stratification to both --vcf1 and --vcf2", default=0, type=int)
    stats_group.add_argument("--period", help="Report results overall and also stratified by repeat unit length (period)", action="store_true")
    ### Optional args ###
    option_group = parser.add_argument_group("Optional arguments")
    option_group.add_argument("--verbose", help="Print helpful debugging info", action="store_true")
    args = parser.parse_args()
    return args

def main(args):
    ### Check and load VCF files ###
    vcfreaders = mergeutils.LoadReaders([args.vcf1, args.vcf2], region=args.region)
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    ### Check inferred type of each is the same
    vcftype = mergeutils.GetVCFType(vcfreaders, args.vcftype)

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = mergeutils.GetMinRecords(current_records, chroms)
    done = mergeutils.DoneReading(current_records)

    ### Load shared samples ###
    samples = mergeutils.GetSharedSamples(vcfreaders)
    print(samples)

    ### Determine FORMAT fields we should look for ###
    format_fields, format_binsizes = GetFormatFields(args.stratify_fields, args.stratify_binsizes, args.stratify_file, vcfreaders)

    ### Keep track of data to summarize at the end ###
    results_dir = {
        "chrom": [],
        "start": [],
        "period": [],
        "sample": [],
        "gtstring1": [],
        "gtstring2": [],
        "gtsum1": [],
        "gtsum2": [],
        "metric-conc": [],
        "metric-alleleconc": []
    }
    for ff in format_fields:
        results_dir[ff+"1"] = []
        results_dir[ff+"2"] = []

    while not done:
        if args.verbose: mergeutils.PrintCurrentRecords(current_records, is_min)
        if mergeutils.CheckMin(is_min): return 1
        # TODO, take info in for each of them and add to results_dir
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        is_min = mergeutils.GetMinRecords(current_records, chroms)
        done = mergeutils.DoneReading(current_records)

    # TODO load results_dir to pandas dataframe
    # TODO make all outputs base on the df
    return 0 

if __name__ == "__main__":  # pragma: no cover
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)

