#!/usr/bin/env python3
"""
Tool for comparing two STR callsets
"""

# Load external libraries
import argparse
import os
import numpy as np
import pandas as pd
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
    option_group.add_argument("--numrecords", help="For debugging, only process this many records", type=int)
    args = parser.parse_args()
    return args

def UpdateComparisonResults(record1, record2, format_fields, samples, results_dir):
    r"""Extract comparable results from a pair of VCF records

    Parameters
    ----------
    record1 : trh.TRRecord
       First record to compare
    record2 : trh.TRRecord
       Second record to compare
    format_fields : list of str
       List of format fields to extract
    samples : list of str
       List of samples to consider
    results_dir : dict
       Results dictionary to update. Need to update:
       chrom, start, period, sample, gtstring1, gtstring2
       gtsum1, gtsum2, metric-conc-seq, metric-conc-len
       for each matching sample
    """
    # Extract shared info
    chrom = record1.vcfrecord.CHROM
    pos = record1.vcfrecord.POS
    for sample in samples:
        call1 = record1.vcfrecord.genotype(sample)
        call2 = record2.vcfrecord.genotype(sample)
        if not (call1.called and call2.called): continue # Skip if not called in both
        gt_string_1 = record1.GetStringGenotype(call1)
        gt_len_1 = record1.GetLengthGenotype(call1)
        gt_string_2 = record2.GetStringGenotype(call2)
        gt_len_2 = record2.GetLengthGenotype(call2)
        # Make sure gts are same ploidy. If not give up
        if len(gt_string_1) != len(gt_string_2) or len(gt_len_1) != len(gt_len_2):
            raise ValueError("Found sample %s of different ploidy at %s:%s"%(sample, chrom, pos))
        # Update results_dir
        results_dir["chrom"].append(chrom)
        results_dir["start"].append(pos)
        results_dir["period"].append(len(record1.motif))
        results_dir["sample"].append(sample)
        results_dir["gtstring1"].append(",".join(gt_string_1))
        results_dir["gtstring2"].append(",".join(gt_string_2))
        results_dir["gtsum1"].append(sum(gt_len_1))
        results_dir["gtsum2"].append(sum(gt_len_2))
        results_dir["metric-conc-seq"].append(int(all([(gt_string_1[i]==gt_string_2[i]) for i in range(len(gt_string_1))])))
        results_dir["metric-conc-len"].append(int(all([(gt_len_1[i]==gt_len_2[i]) for i in range(len(gt_len_1))])))
        for ff in format_fields:
            results_dir[ff+"1"].append("TODO") # TODO
            results_dir[ff+"2"].append("TODO") # TODO

def main(args):
    ### Check and load VCF files ###
    vcfreaders = mergeutils.LoadReaders([args.vcf1, args.vcf2], region=args.region)
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    ### Check inferred type of each is the same
    vcftype = mergeutils.GetVCFType(vcfreaders, args.vcftype)

    ### Set up harmonizers ###
    tr_harmonizers = [trh.TRRecordHarmonizer(args.vcf1, vcftype=vcftype), \
                      trh.TRRecordHarmonizer(args.vcf2, vcftype=vcftype)]

    ### Load shared samples ###
    samples = mergeutils.GetSharedSamples(vcfreaders)

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
        "metric-conc-seq": [],
        "metric-conc-len": [],
    }
    for ff in format_fields:
        results_dir[ff+"1"] = []
        results_dir[ff+"2"] = []

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = mergeutils.GetMinRecords(current_records, chroms)
    done = mergeutils.DoneReading(current_records)
    num_records = 0
    while not done:
        if args.numrecords is not None and num_records > args.numrecords: break
        if args.verbose: mergeutils.PrintCurrentRecords(current_records, is_min)
        if mergeutils.CheckMin(is_min): return 1
        if all([is_min]):
            UpdateComparisonResults(tr_harmonizers[0].HarmonizeRecord(current_records[0]), \
                                    tr_harmonizers[1].HarmonizeRecord(current_records[1]), \
                                    format_fields, samples, results_dir)
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        is_min = mergeutils.GetMinRecords(current_records, chroms)
        done = mergeutils.DoneReading(current_records)
        num_records += 1

    data = pd.DataFrame(results_dir)
    print(data.head())
    # TODO load results_dir to pandas dataframe
    # TODO make all outputs base on the df
    return 0 

if __name__ == "__main__":  # pragma: no cover
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)

