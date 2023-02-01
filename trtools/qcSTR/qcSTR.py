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
import sys
from typing import Dict, List, Optional, Union

import cyvcf2
import pandas as pd
import numpy as np
import sklearn

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

def OutputDiffRefBias(diffs_from_ref, reflens, fname, xlim=(0, 100), \
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

def OutputSampleCallrate(sample_calls: np.ndarray,
                         samples: List[str],
                         fname: str):
    r"""Plot number of calls per sample

    Parameters
    ----------
    sample_calls :
        1D array, number of calls for each sample
    samples:
        List of names of samples, same len as sample_calls
    fname :
        Filename of output plot
    """
    if len(sample_calls.shape) > 1:
        raise ValueError("sample_calls should be 1D")
    if len(samples) != sample_calls.shape[0]:
        raise ValueError("samples should have the same length as"
                         " sample_calls")
    data = pd.DataFrame({"sample": samples, "numcalls": sample_calls})
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


# original idea https://stackoverflow.com/a/39729964/2966505
def _BetterCDF(data: np.ndarray,
               ax: matplotlib.axes.Axes):
    # assumes that axes are already set to (min, max)
    data = np.sort(data)
    x_axis_min, x_axis_max = ax.get_xlim()
    n_points = len(data)
    has_quality_1_point = data[-1] == 1
    if has_quality_1_point:
        # don't print a drop off if the last data point(s)
        # have quality 1
        n_ones = sum(data == data[-1])
        data = np.hstack((
            [x_axis_min],
            data[0:(len(data) - n_ones)],
            [x_axis_max]
        ))
        ys = np.hstack((
            [1],
            np.arange(n_points - 1, n_ones - 1, -1) / n_points,
            [n_ones / n_points]
        ))
    else:
        data = np.hstack((
            [x_axis_min],
            data,
            [x_axis_max]
        ))
        ys = np.hstack((
            [1],
            np.arange(n_points - 1, -1, -1) / n_points,
            [0]
        ))
    #ax.step(data, ys)#, where='post')
    ax.step(data, ys, where='post')


def _OutputQualityHist(
        data: np.ndarray,
        fname: str,
        dist_name: str,
        strat_names: Optional[List[str]] = None):
    #data is 1D if there is no stratification,
    #2D if there is, with stratification over rows
    #and strat_names being set to the names of each row
    spacing = 5e-3
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(np.nanmin(data) - spacing, np.nanmax(data) + spacing)
    if len(data.shape) == 1:
        _BetterCDF(data, ax)
    else: # assume dict
        names = []
        for stratum, name in enumerate(strat_names):
            _BetterCDF(data[stratum, ~np.isnan(data[stratum, :])], ax)
            names.append(name)
        ax.legend(names)
    ax.set_xlabel("Quality", size=15)
    ax.set_ylabel("% of {} with at least this quality".format(dist_name), size=15)
    fig.savefig(fname)
    plt.close()


def OutputQualityPerSample(per_sample_data: np.ndarray, fname: str):
    """Plot quality of calls per sample

    Parameters
    ----------
    per_sample_data:
        1D array of the average quality for each sample, defined as the
        average of qualities of calls across all loci at that sample
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_sample_data, fname, "samples")


def OutputQualityPerLocus(per_locus_data: np.ndarray, fname: str):
    """Plot quality of calls per locus

    Parameters
    ----------
    per_locus_data:
        1D array of an average quality for each locus, defined as the
        average of qualities of calls across all samples at that locus
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_locus_data, fname, "loci")


def OutputQualityPerCall(per_call_data: np.ndarray, fname: str):
    """Plot quality of calls as one distribution,
    irrespective of which sample or locus they came from.

    Parameters
    ----------
    per_call_data:
        1D array of the qualities of all calls
    fname :
        Location to save the output plot
    """
    _OutputQualityHist(per_call_data, fname, "calls")


def OutputQualitySampleStrat(
        per_call_data: np.ndarray,
        samples: List[str],
        fname: str):
    """Plot quality of calls, one line for each sample

    Parameters
    ----------
    per_call_data:
        2D array of qualities of calls where each row is a locus and each col
        is a sample.
    samples:
        List of the names of samples
    fname :
        Location to save the output plot
    """
    if len(per_call_data.shape) != 2:
        raise ValueError("per_call_data should be 2D")
    if len(samples) != per_call_data.shape[1]:
        raise ValueError("samples should have the same length as"
                         " the number of cols in per_call_data")
    _OutputQualityHist(per_call_data.T, fname, "calls", strat_names=samples)


def OutputQualityLocusStrat(
        per_call_data: np.ndarray,
        loci: List[str],
        fname: str):
    """Plot quality of calls, one line for each locus

    Parameters
    ----------
    per_call_data:
        2D array of qualities of calls where each row is a locus and each col
        is a sample.
    loci:
        List of the IDs of loci
    fname :
        Location to save the output plot
    """
    if len(per_call_data.shape) != 2:
        raise ValueError("per_call_data should be 2D")
    if len(loci) != per_call_data.shape[0]:
        raise ValueError("loci should have the same length as"
                         " the number of rows in per_call_data")
    _OutputQualityHist(per_call_data, fname, "calls", strat_names=loci)

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf", help="VCF file to analyze.", type=str, required=True)
    req_group.add_argument("--out", help="Output prefix for files generated", type=str, required=True)
    inp_group = parser.add_argument_group("Optional input arguments")

    vcftype_options = [str(item) for item in trh.VcfTypes.__members__]
    vcftype_options.append("auto")
    inp_group.add_argument(
        "--vcftype",
        type=str,
        help="Which type of VCF to restrict the input to, or 'auto' for no"
             " restrction",
        default="auto",
        choices = vcftype_options
    )

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
        help=("Exclude no-calls and calls without quality scores from quality"
              " graph distributions instead of the default, which is to "
              "include them as zero quality calls. "
              "Setting this can cause the plotting to crash if it reduces the"
              " number of valid calls (in a strata) to <= 1")
    )
    # TODO add n loci argument
    refbias_group = parser.add_argument_group("Reference bias plot options")
    refbias_group.add_argument(
        "--refbias-metric",
        type=str,
        default="mean",
        help=("Which metric to use for the y-axis on "
              "the reference bias plot."),
        choices=['mean', 'median']
    )
    refbias_group.add_argument(
        "--refbias-mingts",
        type=int,
        default=100,
        help=("Don't compute points for the reference bias plot "
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

    if not os.path.exists(os.path.dirname(os.path.abspath(args.out))):
        common.WARNING("Error: The directory which contains the output location {} does"
                       " not exist".format(args.out))
        return 1

    if os.path.isdir(args.out) and args.out.endswith(os.sep):
        common.WARNING("Error: The output location {} is a "
                       "directory".format(args.out))
        return 1

    # Set up reader and harmonizer
    invcf = utils.LoadSingleReader(args.vcf, checkgz=False)
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
    if args.refbias_binsize < 1:
        common.WARNING("--refbias-binsize must be >=1")
        return 1
    if args.refbias_mingts < 0:  # allow for 0 mingts as a synonym for 1
        common.WARNING("--refbias-mingts must be >=1")
        return 1
    if args.refbias_xrange_min >= args.refbias_xrange_max:
        common.WARNING("--refbias-xrange-min ({}) cannot be >= --refbias-xrange-max ({})".format(
            args.refbias_xrange_min, args.refbias_xrange_max))
        return 1

    # Load samples
    if args.samples:
        sample_list = [item.strip()
                       for item
                       in open(args.samples, "r").readlines()]
        sample_index = np.isin(np.array(invcf.samples), sample_list)
        sample_list = list(np.array(invcf.samples)[sample_index])
    else:
        sample_list = invcf.samples
        sample_index = np.ones(len(sample_list), dtype=bool)

    # Figure out which quality plot to produce by default
    default_quality = False
    if len(args.quality) == 0 and harmonizer.HasQualityScore():
        default_quality = True
        if len(sample_list) <= 5:
            args.quality = [_QualityTypes.sample_stratified.value]
        else:
            args.quality = [_QualityTypes.per_locus.value]

    # Set up data to keep track of
    sample_calls = np.zeros(len(sample_list))
    chrom_calls = {} # chrom->numcalls
    diffs_from_ref_bp = [] # for each allele call, keep track of diff (bp) from ref
    diffs_from_ref_unit = [] # for each allele call, keep track of diff
                             # (repeat units) from ref
    reflens_bp = [] # for each allele call, keep track of reference length (bp)
    if _QualityTypes.per_locus.value in args.quality:
        per_locus_data = []
    if _QualityTypes.per_sample.value in args.quality:
        per_sample_total_qual = np.zeros(len(sample_list))
    if (_QualityTypes.per_call.value in args.quality or
            _QualityTypes.sample_stratified.value in args.quality or
            _QualityTypes.locus_stratified.value in args.quality):
        per_call_data = []
    if _QualityTypes.locus_stratified.value in args.quality:
        locus_ids = []


    # read the vcf
    numrecords = 0
    for trrecord in harmonizer:
        if args.numrecords is not None and numrecords >= args.numrecords: break
        if args.period is not None and len(trrecord.motif) != args.period: continue

        chrom = trrecord.chrom
        if chrom not in chrom_calls:
            chrom_calls[chrom] = 0
        allele_counts = trrecord.GetAlleleCounts(uselength=True,
                                                 sample_index=sample_index)

        idx_gts = trrecord.GetGenotypeIndicies()[sample_index, :-1]
        nocall = np.full((1, idx_gts.shape[1]), -1)
        calls = ~np.all(idx_gts == nocall, axis=1)
        sample_calls += calls
        chrom_calls[chrom] += np.sum(calls)

        if len(args.quality) != 0:
            quality_scores = trrecord.GetQualityScores()[sample_index, :]
            quality_scores[~calls] = np.nan
            if not args.quality_ignore_no_call:
                quality_scores[np.isnan(quality_scores)] = 0
            else:
                quality_idxs = ~np.isnan(quality_scores)

        if _QualityTypes.per_sample.value in args.quality:
            if not args.quality_ignore_no_call:
                per_sample_total_qual += quality_scores.reshape(-1)
            else:
                per_sample_total_qual[quality_idxs.reshape(-1)] += \
                    quality_scores[quality_idxs].reshape(-1)
        if _QualityTypes.per_locus.value in args.quality:
            if not args.quality_ignore_no_call:
                per_locus_data.append(np.mean(quality_scores))
            else:
                per_locus_data.append(np.mean(quality_scores[quality_idxs]))
        if (_QualityTypes.sample_stratified.value in args.quality or
                _QualityTypes.locus_stratified.value in args.quality or
                _QualityTypes.per_call.value in args.quality):
            per_call_data.append(quality_scores)
        if _QualityTypes.locus_stratified.value in args.quality:
            locus_ids.append(trrecord.record_id)

        for allele in allele_counts.keys():
            allelediff_unit = allele - trrecord.ref_allele_length
            count = allele_counts[allele]
            reflens_bp.extend([trrecord.ref_allele_length*len(trrecord.motif)]*count)
            diffs_from_ref_unit.extend([allelediff_unit]*count)
            diffs_from_ref_bp.extend([allelediff_unit*len(trrecord.motif)]*count)

        numrecords += 1

    # now rows are loci, cols are samples
    if (_QualityTypes.sample_stratified.value in args.quality or
            _QualityTypes.locus_stratified.value in args.quality or
            _QualityTypes.per_call.value in args.quality):
        per_call_data = np.concatenate(per_call_data, axis=1).T
        if not args.quality_ignore_no_call:
            per_call_data[np.isnan(per_call_data)] = 0

    print("Producing " + args.out + "-diffref-bias.pdf ... ", end='',
          flush=True)
    OutputDiffRefBias(diffs_from_ref_bp, reflens_bp, args.out + "-diffref-bias.pdf", \
                      xlim=(args.refbias_xrange_min, args.refbias_xrange_max), \
                      mingts=args.refbias_mingts, metric=args.refbias_metric, \
                      binsize=args.refbias_binsize)
    if len(sample_list) > 1:
        print("Done.\nProducing " + args.out + "-sample-callnum.pdf ... ",
              end='', flush=True)
        OutputSampleCallrate(sample_calls, sample_list, args.out+"-sample-callnum.pdf")
        print("Done.")
    else:
        print("Done.\nOnly one sample, so skipping " + args.out + "-sample-callnum.pdf ...")
    if 1 < len(list(chrom for chrom, value in chrom_calls.items()
                    if value > 0)):
        print("Producing " + args.out + "-chrom-callnum.pdf ... ", end='',
              flush=True)
        OutputChromCallrate(chrom_calls, args.out+"-chrom-callnum.pdf")
        print("Done.\n", end='')
    else:
        print("Only one chromosome, so skipping " + args.out + "-chrom-callnum.pdf ...")
    print("Producing " + args.out + "-diffref-histogram.pdf ... ", end='',
          flush=True)
    OutputDiffRefHistogram(diffs_from_ref_unit, args.out + "-diffref-histogram.pdf")
    print("Done.")

    if default_quality:
        def quality_output_loc(quality_value):
            return args.out+"-quality.pdf"
    else:
        def quality_output_loc(quality_value):
            return args.out+"-quality-{}.pdf".format(quality_value)

    prior_qual_plot = False
    if _QualityTypes.per_sample.value in args.quality:
        print("Producing " +
              quality_output_loc(_QualityTypes.per_sample.value) +
              " ... ", end='', flush=True)
        # turn totals into means
        if not args.quality_ignore_no_call:
            per_sample_total_qual /= numrecords
        else:
            per_sample_total_qual /= sample_calls
        OutputQualityPerSample(per_sample_total_qual,
                               quality_output_loc(_QualityTypes.per_sample.value))
        prior_qual_plot = True

    if _QualityTypes.sample_stratified.value in args.quality:
        if prior_qual_plot:
            print("Done.")
        print("Producing " +
              quality_output_loc(_QualityTypes.sample_stratified.value) +
              " ... ", end='', flush=True)
        OutputQualitySampleStrat(
            per_call_data,
            sample_list,
            quality_output_loc(_QualityTypes.sample_stratified.value))
        prior_qual_plot = True

    if _QualityTypes.per_locus.value in args.quality:
        if prior_qual_plot:
            print("Done.")
        print("Producing " +
              quality_output_loc(_QualityTypes.per_locus.value) +
              " ... ", end='', flush=True)
        OutputQualityPerLocus(np.array(per_locus_data),
                              quality_output_loc(_QualityTypes.per_locus.value))
        prior_qual_plot = True

    if _QualityTypes.locus_stratified.value in args.quality:
        if prior_qual_plot:
            print("Done.")
        print("Producing " +
              quality_output_loc(_QualityTypes.locus_stratified.value) +
              " ... ", end='', flush=True)
        OutputQualityLocusStrat(
            per_call_data,
            locus_ids,
            quality_output_loc(_QualityTypes.locus_stratified.value))
        prior_qual_plot = True

    if _QualityTypes.per_call.value in args.quality:
        if prior_qual_plot:
            print("Done.")
        print("Producing " +
              quality_output_loc(_QualityTypes.per_call.value) +
              " ... ", end='', flush=True)
        OutputQualityPerCall(
            per_call_data[~np.isnan(per_call_data)].reshape(-1),
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
