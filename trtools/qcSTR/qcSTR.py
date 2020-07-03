#!/usr/bin/env python3
# pylint: disable=C0411,C0413
"""
Tool for generating various QC plots for TR callsets
"""

# Allow making plots even with no x-forward
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Allow plots to be editable in Adobe Illustrator
import matplotlib  # pylint: disable=C4013
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Imports
import argparse
import enum
import os
import statistics as stat
import sys
from typing import Dict, List, Union

import numpy as np
import pandas as pd
import sklearn
import vcf

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from trtools import __version__

class _QualityTypes(enum.Enum):
    """Different quality graphs that can be made"""

    per_locus = 'per-locus'
    sample_stratified = 'sample-stratified'
    per_sample = 'per-sample'
    locus_stratified = 'locus-stratified'
    per_call = 'per-call'

    # Don't include the redundant values
    # in how enums are printed out
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)


def OutputDiffRefHistogram(diffs_from_ref, fname):
    r"""Plot histogram of difference in bp from reference allele

    Parameters
    ----------
    diffs_from_ref : list of int
        Difference of each allele call from the ref allele (in units)
    fname : str
        Filename of output plot
    """
    MAXPOSS = 50 # don't let histogram go beyond this
    minval = max(-1*MAXPOSS, min(diffs_from_ref))
    maxval = min(MAXPOSS, max(diffs_from_ref))
    extremeval = max(abs(minval), abs(maxval))
    bins = np.arange(-1*extremeval, extremeval, 1)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(diffs_from_ref, bins=bins, color="black", edgecolor="white", log=True)
    ax.set_xlabel("Difference from ref (rpt. units)", size=15)
    ax.set_ylabel("Number of alleles", size=15)
    fig.savefig(fname)
    plt.close()

def OutputDiffRefBias(diffs_from_ref, reflens, fname, xlim=(0,100), \
                      mingts=100, metric="mean", binsize=5):
    r"""Plot reflen vs. mean difference from ref bias plot

    Parameters
    ----------
    diffs_from_ref : list of int
        Difference of each allele call from the ref allele (in bp)
    reflens : list of int
        List of reference allele lengths for each call (in bp)
    fname : str
        Filename of output plot
    xlim: tuple of int, optional
        Specify the minimum and maximum x-axis range (in bp)
    mingts: int, optional
        Don't plot data points computed based on fewer than
        this many genotypes
    metric: str, optional
        Which metric to plot on the y-axis value. Must be mean or median
    binsize: int, optional
        Size (in bp) of bins on the x-axis.
    """
    data = pd.DataFrame({"diff": diffs_from_ref, "ref": reflens, "count": [1]*len(reflens)})
    data["ref"] = data["ref"].apply(lambda x: int(x/binsize)*binsize)
    if metric == "mean":
        sum_fn = np.mean
    elif metric == "median":
        sum_fn = np.median
    else:
        common.WARNING("Invalid metric ({}) specified. Skipping reference bias plot".format(metric))
        return
    metric = metric.capitalize()
    summ = data.groupby("ref", as_index=False).agg({"diff": sum_fn, "count": len}).sort_values("ref")
    summ = summ[summ["count"]>=mingts] # exclude small counts
    summ = summ[(summ["ref"]>=xlim[0]) & (summ["ref"]<=xlim[1])] # filter by x range
    if summ.shape[0] == 0:
        common.WARNING("No points left to plot in reference bias plot after "
                       "filtering. Skipping")
        return
    common.MSG("Plotting ref bias plot with the following data:")
    common.MSG(summ)
    trcounts = np.cumsum(summ["count"])
    trfreqs = trcounts/np.sum(summ["count"])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(summ["ref"], summ["diff"], marker="o", color="darkblue")
    ax.axhline(y=0, linestyle="dashed", color="gray")
    ax.set_xlabel("Reference length (bp)", size=15)
    ax.set_ylabel("{} diff from ref (bp)".format(metric), size=15)
    ax1 = ax.twinx()
    ax1.plot(summ["ref"], trfreqs, color="darkred")
    ax1.set_ylabel("Cumulative fraction of alleles", size=15)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close()

def OutputSampleCallrate(sample_calls, fname):
    r"""Plot number of calls per sample

    Parameters
    ----------
    sample_calls : dict of str->int
        Number of calls for each sample
    fname : str
        Filename of output plot
    """
    samples = sample_calls.keys()
    data = pd.DataFrame({"sample": samples, "numcalls": [sample_calls[key] for key in samples]})
    #data = data.sort_values("numcalls") # Commented because the order would be incorrect if sorted
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(data.shape[0]), data["numcalls"])
    ax.set_xticks(range(data.shape[0]))
    ax.set_xticklabels(samples, rotation=90)
    ax.set_ylabel("Number of calls", size=15)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close()

def OutputChromCallrate(chrom_calls, fname):
    r"""Plot number of calls per chromosome

    Parameters
    ----------
    chrom_calls : dict of str->int
        Number of calls for each chromosome
    fname : str
        Filename of output plot
    """
    chroms = sorted(chrom for chrom in chrom_calls.keys() 
                    if chrom_calls[chrom] > 0)
    counts = [chrom_calls[chrom] for chrom in chroms]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.bar(range(len(counts)), counts)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(chroms, rotation=90)
    ax.set_ylabel("Number of calls", size=15)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close()


# from https://stackoverflow.com/a/39729964/2966505
def _BetterReverseCDF(data: List[float],
                      ax: matplotlib.axes.Axes):
    # assumes that axes are already set to (max, min)
    data = np.sort(data)[::-1]
    x_axis_max, x_axis_min = ax.get_xlim()
    data = np.hstack(([x_axis_max], data, [x_axis_min]))
    ys = (np.hstack((np.arange(0, len(data) - 1), [1])) /
          np.float(len(data) - 2))
    ax.step(data, ys, where='post')


def _OutputQualityHist(
        data: Union[List[float], Dict[str, List[float]]],
        fname: str,
        dist_name: str):
    spacing = 5e-3
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if isinstance(data, list):
        ax.set_xlim(max(data) + spacing, min(data) - spacing)
        _BetterReverseCDF(data, ax)
    else: # assume dict
        max_val = 0.0
        min_val = 1.0
        for dist in data.values():
            max_val = max(max_val, max(dist))
            min_val = min(min_val, min(dist))
        ax.set_xlim(max_val + spacing, min_val - spacing)
        names = []
        for name, dist in data.items():
            _BetterReverseCDF(dist, ax)
            names.append(name)
        ax.legend(names)
    ax.set_xlabel("Quality", size=15)
    ax.set_ylabel("% of {} with at least this quality".format(dist_name), size=15)
    fig.savefig(fname)
    plt.close()


def OutputQualityPerSample(per_sample_data : List[float], fname: str):
    """Plot quality of calls per sample 

    Parameters
    ----------
    per_sample_data: 
        list of an average quality for each sample, defined as the
        average of qualities of calls across all loci at that sample
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_sample_data, fname, "samples")


def OutputQualityPerLocus(per_locus_data: List[float], fname: str):
    """Plot quality of calls per locus

    Parameters
    ----------
    per_locus_data: 
        list of an average quality for each locus, defined as the
        average of qualities of calls across all samples at that locus
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_locus_data, fname, "loci")


def OutputQualityPerCall(per_call_data : List[float], fname: str):
    """Plot quality of calls as one distribution,
    irrespective of which sample or locus they came from.

    Parameters
    ----------
    per_call_data: 
        List of the qualities of all calls
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_call_data, fname, "calls")


def OutputQualitySampleStrat(
        sample_strat_data: Dict[str, List[float]], 
        fname: str):
    """Plot quality of calls, one line for each sample

    Parameters
    ----------
    sample_strat_data: 
        dict from sample name to a list of the qualities of the calls for that
        sample.
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(sample_strat_data, fname, "calls")


def OutputQualityLocusStrat(
        locus_strat_data: Dict[str, List[float]],
        fname: str):
    """Plot quality of calls, one line for each locus

    Parameters
    ----------
    locus_strat_data: 
        dict from locus ID to a list of the qualities of the calls for that
        locus.
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(locus_strat_data, fname, "calls")

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf", help="VCF file to analyze.", type=str, required=True)
    req_group.add_argument("--out", help="Output prefix for files generated", type=str, required=True)
    inp_group = parser.add_argument_group("Optional input arguments")
    inp_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VcfTypes.__members__], type=str, default="auto")
    inp_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    inp_group.add_argument("--period", help="Only consider repeats with this motif length", type=int)
    quality_group = parser.add_argument_group("Quality plot options")
    quality_group.add_argument(
        "--quality",
        action="append", 
        choices = [option.value for option in
                   _QualityTypes.__members__.values()],
        default = [],
        help = (
            "Which quality plot(s) to produce. May be specified more than "
            " once. See the README for more info"
        )
    )
    quality_group.add_argument(
        "--quality-ignore-no-call",
        action="store_true",
        default=False,
        help=("Exclude no-calls from quality graph distributions instead of "
              "the default, which is to include them as zero quality calls. "
              "Setting this can cause the plotting to crash if it reduces the"
              " number of valid calls (in a strata) to <= 1")
    )
    refbias_group = parser.add_argument_group("Reference bias plot options")
    refbias_group.add_argument(
        "--refbias-metric",
        type=str,
        default="mean",
        help=("Which metric to use for the y-axis on "
              "the reference bias plot. Must be mean or median")
        )
    refbias_group.add_argument(
        "--refbias-mingts",
        type=int,
        default=100,
        help=("Don't plot points in the reference bias plot computed "
              "based on fewer than this many genotypes")
    )
    refbias_group.add_argument(
        "--refbias-xrange-min",
        type=int,
        default=0,
        help=("Minimum x-axis value (bp) to show on the reference bias plot")
        )
    refbias_group.add_argument(
        "--refbias-xrange-max",
        type=int,
        default=100,
        help=("Maximum x-axis value (bp) to show on the reference bias plot")
        )
    refbias_group.add_argument(
        "--refbias-binsize",
        type=int,
        default=5,
        help=("Size (bp) of x-axis bins for the reference bias plot")
    )
    debug_group = parser.add_argument_group("Debug group")
    debug_group.add_argument("--numrecords", help="Only process this many records", type=int)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    if not os.path.exists(args.vcf):
        common.WARNING("The input vcf location %s does not exist"%args.vcf)
        return 1

    containing_dir = os.path.split(os.path.abspath(args.out))[0]
    if not os.path.exists(containing_dir):
        common.WARNING("The directory {} which contains the output location does"
                       " not exist".format(containing_dir))
        return 1

    if os.path.isdir(args.out):
        common.WARNING("The output location {} is a "
                       "directory".format(args.out))
        return 1


    # Set up reader and harmonizer
    invcf = utils.LoadSingleReader(args.vcf, checkgz = False)
    if invcf is None:
        return 1

    if args.vcftype != 'auto':
        harmonizer = trh.TRRecordHarmonizer(invcf, args.vcftype)
    else:
        harmonizer = trh.TRRecordHarmonizer(invcf)

    if len(args.quality) > 0 and not harmonizer.HasQualityScore():
        common.WARNING("Requested a quality plot, but the input vcf doesn't have "
                       "quality scores!")
        return 1

    # Check refbias options
    if args.refbias_metric not in ["mean","median"]:
        common.WARNING("--refbias-metric must be either 'mean' or 'median'")
        return 1
    if args.refbias_binsize < 1:
        common.WARNING("--refbias-binsize must be >=1")
        return 1
    if args.refbias_mingts < 1:
        common.WARNING("--refbias-mingts must be >=1")
        return 1
    if args.refbias_xrange_min >= args.refbias_xrange_max:
        common.WARNING("--refbias-xrange-min ({}) cannot be >= --refbias-xrange-max ({})".format(
            args.refbias_xrange_min, args.refbias_xrange_max))
        return 1

    # Load samples
    if args.samples:
        samplelist = [item.strip()
                      for item
                      in open(args.samples, "r").readlines()
                      if item.strip() in invcf.samples]
    else: samplelist = invcf.samples

    # Figure out which quality plot to produce by default
    default_quality = False
    if len(args.quality) == 0 and harmonizer.HasQualityScore():
        default_quality = True
        if len(samplelist) <= 5:
            args.quality = [_QualityTypes.sample_stratified.value]
        else:
            args.quality = [_QualityTypes.per_locus.value]

    # Set up data to keep track of
    sample_calls = dict([(sample, 0) for sample in samplelist]) # sample->numcalls
    contigs = invcf.contigs
    if len(contigs) == 0:
        common.WARNING("Warning: no contigs found in VCF file.")
    chrom_calls = dict([(chrom, 0) for chrom in contigs]) # chrom->numcalls
    diffs_from_ref = [] # for each allele call, keep track of diff (bp) from ref
    diffs_from_ref_unit = [] # for each allele call, keep track of diff (units) from ref
    reflens = [] # for each allele call, keep track of reference length (bp)
    if _QualityTypes.per_locus.value in args.quality:
        per_locus_data = []
    if _QualityTypes.per_sample.value in args.quality:
        per_sample_data = {}
        for sample in samplelist: 
            per_sample_data[sample] = []
    if _QualityTypes.per_call.value in args.quality:
        per_call_data = []
    if _QualityTypes.sample_stratified.value in args.quality:
        sample_strat_data = {}
        for sample in samplelist: 
            sample_strat_data[sample] = []
    if _QualityTypes.locus_stratified.value in args.quality:
        locus_strat_data = {}

    # read the vcf
    numrecords = 0
    for trrecord in harmonizer:
        if args.numrecords is not None and numrecords >= args.numrecords: break
        if args.period is not None and len(trrecord.motif) != args.period: continue

        record = trrecord.vcfrecord

        # Extract stats
        chrom = record.CHROM
        rl = len(trrecord.ref_allele)
        allele_counts = trrecord.GetAlleleCounts(uselength=False, samplelist=samplelist)

        # Update data
        num_calls = 0
        if _QualityTypes.per_locus.value in args.quality:
            per_locus_data.append([])
        if _QualityTypes.locus_stratified.value in args.quality:
            locus_strat_data[trrecord.record_id] = []

        # loop over sample data
        for call in record:
            s = call.sample
            if s not in samplelist:
                continue
            if call.called:
                sample_calls[s] += 1
                num_calls += 1

            if len(args.quality) == 0:
                continue

            # set non-calls to zero quality
            if call.called:
                quality_score = trrecord.GetQualityScore(call)
            elif args.quality_ignore_no_call:
                continue
            else:
                quality_score = 0

            if _QualityTypes.per_sample.value in args.quality:
                per_sample_data[s].append(quality_score)
            if _QualityTypes.sample_stratified.value in args.quality:
                sample_strat_data[s].append(quality_score)
            if _QualityTypes.per_locus.value in args.quality:
                per_locus_data[-1].append(quality_score)
            if _QualityTypes.locus_stratified.value in args.quality:
                locus_strat_data[trrecord.record_id].append(quality_score)
            if _QualityTypes.per_call.value in args.quality:
                per_call_data.append(quality_score)

        chrom_calls[chrom] = chrom_calls.get(chrom, 0) + num_calls
        for allele in allele_counts.keys():
            allelediff = len(allele)-rl
            count = allele_counts[allele]
            reflens.extend([rl]*count)
            diffs_from_ref.extend([allelediff]*count)
            diffs_from_ref_unit.extend([allelediff/len(trrecord.motif)]*count)

        numrecords += 1

    print("Producing " + args.out + "-diffref-bias.pdf ... ", end='')
    OutputDiffRefBias(diffs_from_ref, reflens, args.out + "-diffref-bias.pdf", \
                      xlim=(args.refbias_xrange_min, args.refbias_xrange_max), \
                      mingts=args.refbias_mingts, metric=args.refbias_metric, \
                      binsize=args.refbias_binsize)
    if len(samplelist) > 1:
        print("Done.\nProducing " + args.out + "-sample-callnum.pdf ... ",
              end='')
        OutputSampleCallrate(sample_calls, args.out+"-sample-callnum.pdf")
        print("Done.")
    else:
        print("Done.\nOnly one sample, so skipping " + args.out + "-sample-callnum.pdf ...")
    if 1 < len(list(chrom for chrom, value in chrom_calls.items() 
                    if value > 0)):
        print("Producing " + args.out + "-chrom-callnum.pdf ... ", end='')
        OutputChromCallrate(chrom_calls, args.out+"-chrom-callnum.pdf")
        print("Done.\n", end='')
    else:
        print("Only one chromosome, so skipping " + args.out + "-chrom-callnum.pdf ...")
    print("Producing " + args.out + "-diffref-histogram.pdf ... ", end='')
    OutputDiffRefHistogram(diffs_from_ref_unit, args.out + "-diffref-histogram.pdf")
    print("Done.")

    if default_quality:
        def quality_output_loc(quality_value):
            return args.out+"-quality.pdf"
    else:
        def quality_output_loc(quality_value):
            return args.out+"-quality-{}.pdf".format(quality_value)

    if _QualityTypes.per_sample.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.per_sample.value) +
              " ... ", end='')
        new_per_sample_data = []
        for sample_data in per_sample_data.values():
            new_per_sample_data.append(stat.mean(sample_data))
        OutputQualityPerSample(new_per_sample_data,
                               quality_output_loc(_QualityTypes.per_sample.value))

    if _QualityTypes.sample_stratified.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.sample_stratified.value) +
              " ... ", end='')
        OutputQualitySampleStrat(sample_strat_data,
                                 quality_output_loc(_QualityTypes.sample_stratified.value))

    if _QualityTypes.per_locus.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.per_locus.value) +
              " ... ", end='')
        new_per_locus_data = []
        for locus_data in per_locus_data:
            new_per_locus_data.append(stat.mean(locus_data))
        OutputQualityPerLocus(new_per_locus_data,
                              quality_output_loc(_QualityTypes.per_locus.value))

    if _QualityTypes.locus_stratified.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.locus_stratified.value) +
              " ... ", end='')
        OutputQualityLocusStrat(locus_strat_data,
                                quality_output_loc(_QualityTypes.locus_stratified.value))

    if _QualityTypes.per_call.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.per_call.value) +
              " ... ", end='')
        OutputQualityPerCall(per_call_data,
                             quality_output_loc(_QualityTypes.per_call.value))

    if len(args.quality) == 0:
        print("This vcf does not have quality scores, so skipping all "
              "quality plots.")

    print("Done.")
    return 0

def run(): # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
