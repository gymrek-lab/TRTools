#!/usr/bin/env python3
"""
Tool for generating various QC plots for TR callsets
"""

# Allow making plots even with no x-forward
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Allow plots to be editable in Adobe Illustrator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Imports
import argparse
import os
import sys
import vcf

# Load local libraries
if __name__ == "qcSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import common
    import tr_harmonizer as trh
else: # pragma: no cover
    import strtools.utils.common as common  # pragma: no cover
    import strtools.utils.tr_harmonizer as trh  # pragma: no cover

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf", help="VCF file to analyze.", type=str, required=True)
    req_group.add_argument("--out", help="Output prefix for files generated", tpe=str, required=True)
    inout_group.add_argument("--vcftype", help="Options=%s"%trh.VCFTYPES.__members__, type=str, default="auto")
    filter_group = parser.add_argument_group("Filtering group")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--period", help="Only consider repeats with this motif length", type=int)
    args = parser.parse_args()
    return args

def main(args):
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    # Set up reader and harmonizer
    invcf = vcf.Reader(filename=args.vcf)
    tr_harmonizer = trh.TRRecordHarmonizer(invcf, vcftype=args.vcftype)

    # Load samples
    if args.samples:
        samplelist = [item.strip() for item in open(args.samples, "r").readlines()]
    else: samplelist = invcf.samples
    
    # Set up data to keep track of
    sample_calls = dict([(sample, 0) for sample in samplelist]) # sample->numcalls
    contigs = invcf.contigs
    if len(contigs) == 0:
        common.MSG("Warning: no contigs found in VCF file.")
    chrom_calls = dict([(chrom, 0) for chrom in contigs]) # chrom->numcalls
    diffs_from_ref = [] # for each allele call, keep track of diff (bp) from ref
    reflens = [] # for each allele call, keep track of reference length (bp)

    for record in invcf:
        chrom = record.CHROM
        trrecord = tr_harmonizer.HarmonizeRecord(record)
        if args.period is not None and len(trrecord.motif) != args.period: continue
        # Extract stats
        rl = len(trrecord.ref_allele)
        allele_counts = trrecord.GetAlleleCounts(uselength=False, samplelist=samplelist)
        called_samples = [item.sample for item in record if item.called]
        # Update data
        num_calls = 0
        for s in called_samples:
            try:
                sample_calls[s] += 1
                num_calls += 1
            except KeyError: pass
        chrom_calls[chrom] = chrom_calls.get(chrom, 0) + num_calls
        for allele in allele_counts.keys():
            allelediff = len(allele)-rl
            count = allele_counts[allele]
            reflens.extend([rl]*count)
            diffs_from_ref.extend([allelediff]*count)

    # TODO output plots: histogram, diff from ref vs. reflen, sample call num, chrom call num

if __name__ == "__main__":  # pragma: no cover
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)

