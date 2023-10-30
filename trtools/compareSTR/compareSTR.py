# pylint: disable=C0411,C0413
"""
Tool for comparing genotypes from two TR VCFs
"""

# Allow making plots even with no x-forward
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Allow plots to be editable in Adobe Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import argparse
import os

# Load external libraries
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scipy.stats
import sys

import trtools.utils.common as common
import trtools.utils.mergeutils as mergeutils
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from trtools import __version__

from typing import List, Any, Callable, Optional


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
    format_bins : List[List[float]]
        List of bin start/stop coords for each FORMAT field
    """
    if format_fields is None or format_binsizes is None:
        return [], []

    def get_formats(vcf):
        formats = []
        for header in vcf.header_iter():
            if header['HeaderType'] == 'FORMAT':
                formats.append(header['ID'])
        return formats

    formats1 = get_formats(vcfreaders[0])
    formats2 = get_formats(vcfreaders[1])

    formats = format_fields.split(",")
    # TODO confirm that these are all Int of Float type and Number == 1
    binsizes = format_binsizes.split(",")
    if len(formats) != len(binsizes):
        raise ValueError("--stratify-formats must be same length as --stratify-binsizes")
    binsizes = [[float(x) for x in item.split(":")] for item in binsizes]
    # TODO confirm len(binsize) == 3 for each binsize
    bins = []
    for start, stop, step in binsizes:
        bins.append(np.arange(start, stop, step).tolist())
        bins[-1].append(stop)

    for fmt in formats:
        check1 = fmt in formats1
        check2 = fmt in formats2
        if format_fileoption == 0 and not (check1 and check2):
            raise ValueError("FORMAT field %s must be present in both VCFs if --stratify-file=0" % fmt)
        if format_fileoption == 1 and not check1:
            raise ValueError("FORMAT field %s must be present in --vcf1 if --stratify-file=1" % fmt)
        if format_fileoption == 2 and not check2:
            raise ValueError("FORMAT field %s must be present in --vcf2 if --stratify-file=2" % fmt)

    return formats, bins


def OutputLocusMetrics(locus_results, outprefix, noplot):
    r"""Output per-locus metrics

    Outputs text file and plot of per-locus metrics
    outprefix + "-locuscompare.tab"
    outprefix + "-locuscompare.pdf"

    Parameters
    ----------
    locus_results: Dict[str, Any]
        The info needed to write the output file
    outprefix : str
        Prefix to name output file
    noplot : bool
        If True, don't output plots
    """
    with open(outprefix + '-locuscompare.tab', 'w') as tabfile:
        tabfile.write('chrom\tstart\tmetric-conc-seq\tmetric-conc-len\tnumcalls\n')
        for chrom, start, metric_conc_seq, metric_conc_len, numcalls in zip(
                locus_results['chrom'],
                locus_results['start'],
                locus_results['metric-conc-seq'],
                locus_results['metric-conc-len'],
                locus_results['numcalls']
        ):
            tabfile.write('{}\t{}\t{}\t{}\t{}\n'.format(
                chrom, start, metric_conc_seq, metric_conc_len, numcalls
            ))

    # Create per-locus plot
    if noplot: return

    fig = plt.figure()
    ax = fig.add_subplot(111)

    nloci = len(locus_results['chrom'])
    if nloci <= 20:
        sort_idx = np.argsort(locus_results['metric-conc-len'])[::-1]
        for key in {'chrom', 'start', 'metric-conc-len'}:
            locus_results[key] = np.array(locus_results[key])[sort_idx]
        ax.scatter(np.arange(nloci), locus_results['metric-conc-len'], color="darkblue")
        ax.set_xticks(np.arange(nloci))
        ax.set_xticklabels(
            ["{}:{}".format(chrom, start) for chrom, start in zip(
                locus_results['chrom'], locus_results['start']
            )], size=12, rotation=90
        )
    else:
        sorted_results = np.sort(locus_results['metric-conc-len'])[::-1]
        ax.scatter(np.arange(nloci), sorted_results, color="darkblue")
        ax.set_xlabel("Successive TR Loci", size=15)
    ax.set_ylabel("Length Concordance", size=15)
    plt.tight_layout()
    fig.savefig(outprefix + "-locuscompare.pdf")
    plt.close()


def OutputSampleMetrics(sample_results, sample_names, outprefix, noplot):
    r"""Output per-sample metrics

    Outputs text file and plot of per-sample metrics
    outprefix + "-samplecompare.tab"
    outprefix + "-samplecompare.pdf"

    Parameters
    ----------
    sample_results : Dict[str, any]
        The info needed to write the output file
    sample_names : List[str]
    outprefix : str
        Prefix to name output file
    noplot : bool
        If True, don't output plots
    """
    sample_results['conc-seq-count'] = \
        sample_results['conc-seq-count'] / sample_results['numcalls']
    sample_results['conc-len-count'] = \
        sample_results['conc-len-count'] / sample_results['numcalls']
    with open(outprefix + '-samplecompare.tab', 'w') as tabfile:
        tabfile.write('sample\tmetric-conc-seq\tmetric-conc-len\tnumcalls\n')
        for idx, sample in enumerate(sample_names):
            tabfile.write('{}\t{}\t{}\t{}\n'.format(
                sample,
                sample_results['conc-seq-count'][idx],
                sample_results['conc-len-count'][idx],
                sample_results['numcalls'][idx]
            ))

    # Create per-locus plot
    if noplot: return
    nsamples = len(sample_names)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if nsamples <= 20:
        sort_idx = np.argsort(sample_results['conc-len-count'])[::-1]
        ax.scatter(np.arange(nsamples),
                   sample_results['conc-len-count'][sort_idx],
                   color="darkblue")
        ax.set_xticks(np.arange(nsamples))
        ax.set_xticklabels(np.array(sample_names)[sort_idx], size=12, rotation=90)
    else:
        sorted_results = np.sort(sample_results['conc-len-count'])[::-1]
        ax.scatter(np.arange(nsamples), sorted_results, color="darkblue")
        ax.set_xlabel("Successive samples", size=15)
    ax.set_ylabel("Length Concordance", size=15)
    plt.tight_layout()
    fig.savefig(outprefix + "-samplecompare.pdf")
    plt.close()


def OutputOverallMetrics(overall_results, format_fields, format_bins, outprefix):
    r"""Output overall accuracy metrics

    Output metrics overall, by period, and by FORMAT bins
    Output results to outprefix+"-overall.tab"

    Parameters
    ----------
    overall_results : Dict[str, Any]
        Info needed to write the tabfile
    format_fields : List[str]
        List of FORMAT fields to stratify by
    format_bins : List[List[float]]
        List of bin start/stop coords for each FORMAT field
    outprefix : str
        Prefix to name output file
    """
    # get periods in sort order
    periods = set(overall_results.keys())
    periods.remove('ALL')
    periods = list(periods)
    periods.sort()
    periods.insert(0, 'ALL')

    def write_format_bin(tabfile, format_bin_results, per, fmt_idx,
                         format_bin_string):
        numcalls = format_bin_results['numcalls']
        if numcalls == 0:
            return

        tabfile.write(str(per))
        tabfile.write('\t')
        for idx in range(len(format_fields)):
            if idx == fmt_idx:
                tabfile.write(format_bin_string)
                tabfile.write('\t')
            else:
                tabfile.write('NA\t')
        tabfile.write('{}\t{}\t{}\t{}\n'.format(
            format_bin_results['conc_seq_count'] / numcalls,
            format_bin_results['conc_len_count'] / numcalls,
            CalcR2(format_bin_results),
            numcalls
        ))

    with open(outprefix + "-overall.tab", "w") as tabfile:
        tabfile.write('period\t')
        for fmt in format_fields:
            tabfile.write(fmt)
            tabfile.write('\t')
        tabfile.write("concordance-seq\tconcordance-len\tr2\tnumcalls\n")

        for per in periods:
            # write the entry that is not stratified across formats
            write_format_bin(tabfile, overall_results[per]['ALL'],
                             per, None, None)
            # stratify across formats
            for fmt_idx, (fmt, bins) in enumerate(zip(format_fields, format_bins)):
                for bin_idx in range(len(bins) - 2):
                    bin_string = "[{}, {})".format(bins[bin_idx],
                                                   bins[bin_idx + 1])
                    write_format_bin(tabfile,
                                     overall_results[per][fmt][bins[bin_idx]],
                                     per,
                                     fmt_idx,
                                     bin_string)
                bin_string = "[{}, {}]".format(bins[-2],
                                               bins[-1])
                write_format_bin(tabfile,
                                 overall_results[per][fmt][bins[-2]],
                                 per,
                                 fmt_idx,
                                 bin_string)


def GetBubbleLegend(coordinate_counts):
    r"""Get three good bubble legend sizes to use

    They should be nice round numbers spanning the orders of magnitude of the dataset

    Parameters
    ----------
    coordinate_counts :
        set of counts for coordinates in the graph

    Returns
    -------
    legend_values : list of int
        List of three or fewer representative sample sizes to use for bubble legend
    """
    if len(coordinate_counts) <= 3: return list(coordinate_counts)  # if only three values, return three of them
    # Determine if we do log10 or linear scale
    minval = min(coordinate_counts)
    maxval = max(coordinate_counts)
    if maxval / minval > 10:
        # Do log10 scale
        # Find max power of 10
        max10 = int(np.log10(maxval))
        # Find min power of 10
        min10 = int(np.log10(minval))
        # Find power of 10 in between
        mid10 = int((max10 + min10) / 2)
        return sorted(list(set([10 ** min10, 10 ** mid10, 10 ** max10])))
    else:
        # Do linear scale
        mid = int((minval + maxval) / 2)
        return sorted(list(set([minval, mid, maxval])))


def OutputBubblePlot(bubble_results, outprefix, minval=None, maxval=None):
    r"""Output bubble plot of gtsum1 vs. gtsum2

    Parameters
    ----------
    bubble_results :
        counts of sum1 vs sum2
    outprefix : str
        Prefix to name output file
    """
    # get periods in sort order
    periods = set(bubble_results.keys())
    periods.remove('ALL')
    periods = list(periods)
    periods.sort()
    periods.insert(0, 'ALL')

    for per in periods:
        per_results = bubble_results[per]
        x_vals = [x for x, y in per_results.keys()]
        y_vals = [y for x, y in per_results.keys()]
        scale = 10000 / np.mean(list(per_results.values()))
        if minval is None:
            minval = min(min(x_vals), min(y_vals))
        if maxval is None:
            maxval = max(max(x_vals), max(y_vals))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Plot (0,0) separately so everything else is in front of it
        if (0, 0) in per_results:
            ax.scatter(0, 0,
                       s=np.sqrt(per_results[(0, 0)] * scale),
                       color="darkblue",
                       alpha=0.5)
        for coord, count in per_results.items():
            if coord == (0, 0):
                continue
            ax.scatter(coord[0], coord[1],
                       s=np.sqrt(count * scale),
                       color="darkblue",
                       alpha=0.5)
        ax.set_xlabel("sum # repeats - file 1\n(diff from ref)", size=15)
        ax.set_ylabel("sum # repeats - file 2\n(diff from ref)", size=15)
        ax.plot([minval, maxval], [minval, maxval], linestyle="dashed",
                color="gray", alpha=0.75)
        ax.set_xlim(left=minval, right=maxval)
        ax.set_ylim(bottom=minval, top=maxval)
        ax.axhline(y=0, linestyle="dashed", color="gray", alpha=0.75)
        ax.axvline(x=0, linestyle="dashed", color="gray", alpha=0.75)
        # plot dummy points for legend
        legend_values = GetBubbleLegend(set(per_results.values()))
        xval = (maxval - minval) / 10 + minval
        for i, val in enumerate(legend_values):
            step = (maxval - minval) / 15
            yval = step * (i + 3)
            ax.scatter([xval], [yval], color="darkblue", s=np.sqrt(val * scale))
            ax.annotate(val, xy=(xval + step, yval))
        fig.savefig(outprefix + "-bubble-period%s.pdf" % per,
                    bbox_inches='tight')
        plt.close()


def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf1", help="First VCF file to compare (must be sorted, bgzipped, and indexed)", type=str,
                           required=True)
    req_group.add_argument("--vcf2", help="Second VCF file to compare (must be sorted, bgzipped, and indexed)",
                           type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    ### Options for filtering input ###
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--region", help="Restrict to this region chrom:start-end", type=str)
    ### Stratify results ###
    stats_group = parser.add_argument_group("Metrics to stratify results")
    stats_group.add_argument("--stratify-fields", help="Comma-separated list of FORMAT fields to stratify by", type=str)
    stats_group.add_argument("--stratify-binsizes",
                             help="Comma-separated list of min:max:binsize to stratify each field on. Must be same length as --stratify-fields.",
                             type=str)
    stats_group.add_argument("--stratify-file",
                             help="Set to 1 to stratify based on --vcf1. Set to 2 to stratify based on --vcf2. Set to 0 to apply stratification to both --vcf1 and --vcf2",
                             default=0, type=int)
    stats_group.add_argument("--period",
                             help="Report results overall and also stratified by repeat unit length (period)",
                             action="store_true")
    ### Plotting args ###
    plot_group = parser.add_argument_group("Plotting options")
    plot_group.add_argument("--bubble-min", help="Minimum x/y axis value to display on bubble plots", type=int)
    plot_group.add_argument("--bubble-max", help="Maximum x/y axis value to display on bubble plots", type=int)
    ### Optional args ###
    option_group = parser.add_argument_group("Optional arguments")
    option_group.add_argument("--verbose", help="Print helpful debugging info", action="store_true")
    option_group.add_argument("--numrecords", help="For debugging, only process this many records", type=int)
    option_group.add_argument("--noplot", help="Don't output any plots. Only produce text output", action="store_true")
    option_group.add_argument("--vcftype1",
                              help="Type of --vcf1. Options=%s" % [str(item) for item in trh.VcfTypes.__members__],
                              type=str, default="auto")
    option_group.add_argument("--vcftype2",
                              help="Type of --vcf2. Options=%s" % [str(item) for item in trh.VcfTypes.__members__],
                              type=str, default="auto")
    option_group.add_argument("--ignore-phasing", help="Treat all calls as if they are unphased", action="store_true")
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version='{version}'.format(version=__version__))
    args = parser.parse_args()
    return args


def NewOverallFormatBin():
    """
    Return an empty bin for the overall dictionary.

    Returns
    -------
    Dict[str, Union[int, float]] :
        Contains the fields:
        conc_len_count
        conc_seq_cont
        numcalls
        total_len_1
        total_len_2
        total_len_11
        total_len_12
        total_len_22
    """
    return {
        'conc_seq_count': 0,
        'conc_len_count': 0,
        'numcalls': 0,
        'total_len_1': 0,
        'total_len_2': 0,
        'total_len_11': 0,
        'total_len_12': 0,
        'total_len_22': 0
    }


def CalcR2(format_bin_results):
    """
    Calculate the squared (pearson) correlation coefficient
    for the values in this bin.

    Calculation is done using the formulas:
        n = numcalls
        var(X) = sum(X_i**2)/n - [sum(X_i)/n]**2
        covar(X,Y) = sum(X_i*Y_i)/n - sum(X_i)/n * sum(Y_i)/n
        r^2 = covar(X,Y)**2/(var(X) * var(Y))

    Parameters
    ----------
    format_bin_results : Dict[str, int]
        See the method NewOverallForamtBin
    
    Returns
    -------
    float:
        r^2, or np.nan if one of the two vcfs has
        no variance in this format bin
    """
    f = format_bin_results
    n = f['numcalls']
    var1 = f['total_len_11'] / n - (f['total_len_1'] / n) ** 2
    var2 = f['total_len_22'] / n - (f['total_len_2'] / n) ** 2
    if var1 == 0 or var2 == 0:
        return np.nan
    covar = f['total_len_12'] / n - f['total_len_1'] * f['total_len_2'] / n ** 2
    return covar ** 2 / (var1 * var2)


def NewOverallPeriod(format_fields, format_bins):
    """
    Return an empty dictionary containing
    bins for each format stratification and for
    'ALL' (no format stratification).

    Returns
    -------
    The empty dictionary.
    """
    period_dict = {
        'ALL': NewOverallFormatBin()
    }
    for fmt, bins in zip(format_fields, format_bins):
        period_dict[fmt] = {}
        for _bin in bins[:-1]:
            period_dict[fmt][_bin] = NewOverallFormatBin()
    return period_dict


def UpdateComparisonResults(record1, record2, sample_idxs,
                            ignore_phasing,
                            stratify_by_period,
                            format_fields, format_bins, stratify_file,
                            overall_results, locus_results, sample_results,
                            bubble_results):
    r"""Extract comparable results from a pair of VCF records

    Parameters
    ----------
    record1 : trh.TRRecord
       First record to compare
    record2 : trh.TRRecord
       Second record to compare
    sample_idxs : list of np.array
        Two arrays, one for each vcf
        Each array is a list of indicies so that
        vcf1.samples[index_array1] == vcf2.samples[index_array2]
        and that this is the set of shared samples
    stratify_by_period : bool
        If True, also stratify results by period
    format_fields : list of str
       List of format fields to extract
    format_bins : List[List[float]]
        List of bin start/stop coords for each FORMAT field
    stratify_file : {0, 1, 2}
        Specify whether to apply FORMAT stratification to both files (0), or only (1) or (2)
    overall_results : dict
        Period and format nested dictionary to update.
    locus_results : dict
       Locus-stratified results dictionary to update.
    sample_results : dict
       Sample-stratified results dictionary to update.
    bubble_results : dict
        dictionary of counts to update
    """
    # Extract shared info
    chrom = record1.chrom
    pos = record1.pos
    period = len(record1.motif)
    reflen = len(record1.ref_allele) / period

    both_called = np.logical_and(
        record1.GetCalledSamples()[sample_idxs[0]],
        record2.GetCalledSamples()[sample_idxs[1]]
    )
    numcalls = np.sum(both_called)
    if numcalls == 0:
        return

    locus_results["chrom"].append(chrom)
    locus_results["start"].append(pos)
    locus_results["numcalls"].append(numcalls)
    sample_results['numcalls'] += both_called

    # build this so indexing later in the method is more intuitive
    called_sample_idxs = []
    for sample_idx in sample_idxs:
        called_sample_idxs.append(sample_idx[both_called])

    ploidies1 = record1.GetSamplePloidies()[called_sample_idxs[0]]
    ploidies2 = record2.GetSamplePloidies()[called_sample_idxs[1]]
    # Make sure gts are same ploidy. If not give up
    if not np.all(ploidies1 == ploidies2):
        raise ValueError("Found sample(s) of different ploidy at %s:%s" % (chrom, pos))

    gts_string_1 = record1.GetStringGenotypes()[called_sample_idxs[0], :]
    gts_string_2 = record2.GetStringGenotypes()[called_sample_idxs[1], :]

    # Make sure same phasedness between calls. If not, give up
    if ignore_phasing:
        all_unphased = True
    else:
        unphased = (gts_string_1[:, -1] == '0') & (gts_string_2[:, -1] == '0')
        all_unphased = np.all(unphased)
        if not (all_unphased or np.all(~unphased)):
            raise ValueError("Found sample(s) with different phasedness at %s:%s" % (chrom, pos))

    gts_string_1 = gts_string_1[:, :-1]
    gts_string_2 = gts_string_2[:, :-1]
    if all_unphased:
        gts_string_1 = np.sort(gts_string_1, axis=1)
        gts_string_2 = np.sort(gts_string_2, axis=1)
    conc_seq = np.all(gts_string_1 == gts_string_2, axis=1)

    locus_results["metric-conc-seq"].append(np.sum(conc_seq) / numcalls)
    sample_results['conc-seq-count'][both_called] += conc_seq

    gts_length_1 = record1.GetLengthGenotypes()[called_sample_idxs[0], :-1]
    gts_length_2 = record2.GetLengthGenotypes()[called_sample_idxs[1], :-1]
    if all_unphased:
        gts_length_1 = np.sort(gts_length_1, axis=1)
        gts_length_2 = np.sort(gts_length_2, axis=1)
    conc_len = np.all(gts_length_1 == gts_length_2, axis=1)

    locus_results["metric-conc-len"].append(np.sum(conc_len) / numcalls)
    sample_results['conc-len-count'][both_called] += conc_len

    sum_length_1 = np.sum(gts_length_1 - reflen, axis=1)
    sum_length_2 = np.sum(gts_length_2 - reflen, axis=1)

    outer_keys = ['ALL']
    if stratify_by_period:
        outer_keys.append(period)
        if period not in overall_results:
            overall_results[period] = NewOverallPeriod(format_fields, format_bins)
            if bubble_results:
                bubble_results[period] = {}

    # handle bubble results
    if bubble_results:
        length_sums = np.stack((sum_length_1, sum_length_2)).T
        coords, counts = np.unique(length_sums, axis=0, return_counts=True)
        for coord, count in zip((tuple(row) for row in coords), counts):
            if coord not in bubble_results['ALL']:
                bubble_results['ALL'][coord] = 0
            if stratify_by_period and coord not in bubble_results[period]:
                bubble_results[period][coord] = 0

            bubble_results['ALL'][coord] += count
            if stratify_by_period:
                bubble_results[period][coord] += count

    # handle overall results
    for key in outer_keys:
        overall_results[key]['ALL']['numcalls'] += numcalls
        overall_results[key]['ALL']['conc_seq_count'] += np.sum(conc_seq)
        overall_results[key]['ALL']['conc_len_count'] += np.sum(conc_len)
        overall_results[key]['ALL']['total_len_1'] += np.sum(sum_length_1)
        overall_results[key]['ALL']['total_len_2'] += np.sum(sum_length_2)
        overall_results[key]['ALL']['total_len_11'] += np.sum(sum_length_1 ** 2)
        overall_results[key]['ALL']['total_len_12'] += np.sum(sum_length_1 * sum_length_2)
        overall_results[key]['ALL']['total_len_22'] += np.sum(sum_length_2 ** 2)

    for fmt, bins in zip(format_fields, format_bins):
        fmt1 = record1.format[fmt][sample_idxs[0], 0]
        fmt2 = record2.format[fmt][sample_idxs[1], 0]
        masks = []
        for idx in range(len(bins) - 2):
            if stratify_file == 0:
                mask = ((fmt1 >= bins[idx]) & (fmt1 < bins[idx + 1]) &
                        (fmt2 >= bins[idx]) & (fmt2 < bins[idx + 1]))
            elif stratify_file == 1:
                mask = (fmt1 >= bins[idx]) & (fmt1 < bins[idx + 1])
            elif stratify_file == 2:
                mask = (fmt2 >= bins[idx]) & (fmt2 < bins[idx + 1])
            masks.append(mask[both_called])

        # last bin has inclusive end
        if stratify_file == 0:
            mask = ((fmt1 >= bins[-2]) & (fmt1 <= bins[-1]) &
                    (fmt2 >= bins[-2]) & (fmt2 <= bins[-1]))
        elif stratify_file == 1:
            mask = (fmt1 >= bins[-2]) & (fmt1 <= bins[-1])
        elif stratify_file == 2:
            mask = (fmt2 >= bins[-2]) & (fmt2 <= bins[-1])
        masks.append(mask[both_called])

        for _bin, mask in zip(bins[:-1], masks):
            ncalls = np.sum(mask)
            if ncalls == 0:
                continue
            conc_seq_count = np.sum(conc_seq[mask])
            conc_len_count = np.sum(conc_len[mask])
            total_len_1 = np.sum(sum_length_1[mask])
            total_len_2 = np.sum(sum_length_2[mask])
            total_len_11 = np.sum(sum_length_1[mask] ** 2)
            total_len_12 = np.sum(sum_length_1[mask] * sum_length_2[mask])
            total_len_22 = np.sum(sum_length_2[mask] ** 2)
            for key in outer_keys:
                overall_results[key][fmt][_bin]['numcalls'] += \
                    ncalls
                overall_results[key][fmt][_bin]['conc_seq_count'] += \
                    conc_seq_count
                overall_results[key][fmt][_bin]['conc_len_count'] += \
                    conc_len_count
                overall_results[key][fmt][_bin]['total_len_1'] += \
                    total_len_1
                overall_results[key][fmt][_bin]['total_len_2'] += \
                    total_len_2
                overall_results[key][fmt][_bin]['total_len_11'] += \
                    total_len_11
                overall_results[key][fmt][_bin]['total_len_12'] += \
                    total_len_12
                overall_results[key][fmt][_bin]['total_len_22'] += \
                    total_len_22


def check_region(contigs1, contigs2, region_str):
    def check_contig(contig):
        if contig not in contigs1 or contig not in contigs2:
            common.WARNING("contig {} was not found in both input vcfs".format(contig))
            return 1
        return 0

    if ':' not in region_str:
        return check_contig(region_str)

    parts = region_str.split(':')
    if not len(parts) == 2:
        common.WARNING("--region should have the format contig:range")
        return 1

    contig, _range = parts
    if check_contig(contig) == 1:
        return 1

    def bad_range():
        common.WARNING("The range portion of --region should have one of "
                       "the forms: 42, -42, 42- or 13-42")
        return 1

    try:
        if '-' not in _range:
            int(_range)
            return 0

        parts = _range.split('-')
        if not len(parts) == 2:
            return bad_range()
        start, end = parts
        if start != '':
            int(start)
        if end != '':
            int(end)
        if end == '' and start == '':
            return bad_range()
        if end != '' and start != '' and int(end) <= int(start):
            common.WARNING("Cannot have range portion of --region "
                           "start-end where end <= start")
            return 1
    except ValueError:  # an int() cast failed
        return bad_range()

    return 0


def handle_overlaps(records: List[Optional[trh.TRRecord]], chrom_indices: List[int], min_chrom_index: int) -> bool:
    """
        This function determines whether (two) records in list are comparable
        Currently only works with record lists which are two records long

        Parameters
        ----------
        records: List[Optional[trh.TRRecord]]
            List of TRRecords whose comparability is to be determined. If any of them is None,
            they are not comparable
        chrom_indices: List[int]
            List of indices of chromosomes of current records
        min_chrom_index: int
            Smallest index in chrom_indices. All records should have the same chrom_index,
            otherwise they are not comparable

        Returns
        -------
        comparable: bool
            Result, that says whether records are comparable
        """
    assert len(records) == 2
    # for now, this is just a constant, but this value might become configurable in future releases
    min_overlap = 1.0

    if any(record is None for record in records):
        return False

    left, right = records[0], records[1]
    if chrom_indices[0] != chrom_indices[1] or \
            chrom_indices[0] != min_chrom_index or \
            chrom_indices[1] != min_chrom_index:
        return False

    left_start, left_end = left.pos, left.end_pos
    right_start, right_end = right.pos, right.end_pos

    overlap = min(left_end, right_end) - max(left_start, right_start) + 1
    # This calculation contains max() + 1 to compensate
    # that both start and end coordinates that are used in previous calculation are inclusive
    comparable = \
        overlap / max(left.ref_allele_length * len(left.motif), right.ref_allele_length * len(right.motif)) \
        >= min_overlap

    if overlap >= 1 and not comparable:
        common.WARNING(f"Records {left.record_id} and {right.record_id} overlap:\n"
                       f"{left.record_id}: {left_start, left_end}\n"
                       f"{right.record_id}: {right_start, right_end},\n"
                       f"but are NOT comparable!")

    return comparable


def main(args):
    if not os.path.exists(os.path.dirname(os.path.abspath(args.out))):
        common.WARNING("Error: The directory which contains the output location {} does"
                       " not exist".format(args.out))
        return 1

    if os.path.isdir(args.out) and args.out.endswith(os.sep):
        common.WARNING("Error: The output location {} is a "
                       "directory".format(args.out))
        return 1

    ### Check and load VCF files ###
    vcfreaders = utils.LoadReaders([args.vcf1, args.vcf2], checkgz=True)
    if vcfreaders is None or len(vcfreaders) != 2:
        return 1
    chroms = utils.GetContigs(vcfreaders[0])

    ### Load shared samples ###
    samples = mergeutils.GetSharedSamples(vcfreaders)
    if len(samples) == 0:
        common.WARNING("No shared smaples found between the vcfs")
        return 1
    if args.samples:
        usesamples = set([item.strip() for item in open(args.samples, "r").readlines()])
        samples = list(set(samples).intersection(usesamples))
    if len(samples) == 0:
        common.WARNING("No shared samples found between the vcfs and the "
                       "--samples file")
        return 1
    samples.sort()
    sample_idxs = []
    for vcf in vcfreaders:
        sort = np.argsort(vcf.samples)
        rank = np.searchsorted(vcf.samples, samples, sorter=sort)
        sample_idxs.append(sort[rank])
    # now we have vcfreaders[i].samples[sample_idxs[i]] == samples

    ### Determine FORMAT fields we should look for ###
    if args.stratify_file is not None and args.stratify_file not in [0, 1, 2]:
        common.MSG("--stratify-file must be 0,1, or 2")
        return 1
    format_fields, format_bins = GetFormatFields(args.stratify_fields, args.stratify_binsizes, args.stratify_file,
                                                 vcfreaders)

    ### Keep track of data to summarize at the end ###
    locus_results = {
        "chrom": [],
        "start": [],
        "numcalls": [],
        "metric-conc-seq": [],
        "metric-conc-len": [],
    }
    sample_results = {
        "numcalls": np.zeros((len(samples)), dtype=int),
        "conc-seq-count": np.zeros((len(samples)), dtype=int),
        "conc-len-count": np.zeros((len(samples)), dtype=int)
    }
    # nested dicts period -> format -> val
    # record running totals so results do not need to be stored in memory
    overall_results = {
        'ALL': NewOverallPeriod(format_fields, format_bins)
    }

    if not args.noplot:
        bubble_results = {'ALL': {}}
    else:
        bubble_results = None

    try:
        vcftype1 = trh.InferVCFType(vcfreaders[0], args.vcftype1)
    except TypeError as te:
        common.WARNING("Error with type of vcf1: " + str(te))
        return 1

    try:
        vcftype2 = trh.InferVCFType(vcfreaders[1], args.vcftype2)
    except TypeError as te:
        common.WARNING("Error with type of vcf2: " + str(te))
        return 1

    if not args.region:
        vcfregions = vcfreaders
    else:
        contigs1 = utils.GetContigs(vcfreaders[0])
        contigs2 = utils.GetContigs(vcfreaders[0])
        if check_region(contigs1, contigs2, args.region) == 1:
            return 1
        vcfregions = [vcfreaders[0](args.region), vcfreaders[1](args.region)]

    ### Walk through sorted readers, merging records as we go ###
    current_records = mergeutils.InitReaders(vcfreaders)
    done = mergeutils.DoneReading(current_records)
    vcf_types = [vcftype1, vcftype2]
    num_records = 0
    compared_records = 0

    while not done:
        if any([item is None for item in current_records]): break
        if args.numrecords is not None and num_records >= args.numrecords: break

        harmonized_records = [trh.HarmonizeRecord(vcf_types[i], current_records[i]) for i in
                              range(len(current_records))]
        # increments contains information about which record should be
        # skipped in next iteration
        increment, comparable = mergeutils.GetIncrementAndComparability(harmonized_records, chroms,
                                                                        handle_overlaps)

        if args.verbose: mergeutils.DebugPrintRecordLocations(current_records, increment)
        if mergeutils.CheckMin(increment): return 1
        if comparable:
            UpdateComparisonResults(*harmonized_records,
                                    sample_idxs,
                                    args.ignore_phasing, args.period,
                                    format_fields, format_bins,
                                    args.stratify_file,
                                    overall_results, locus_results,
                                    sample_results, bubble_results)
            compared_records += 1

        current_records = mergeutils.GetNextRecords(vcfregions, current_records, increment)
        done = mergeutils.DoneReading(current_records)
        num_records += 1

    if compared_records == 0:
        common.WARNING("No comparable records were found, exiting!")
        return 1

    ### Overall metrics ###
    OutputOverallMetrics(overall_results, format_fields, format_bins, args.out)
    if not args.noplot: OutputBubblePlot(bubble_results, args.out, minval=args.bubble_min, maxval=args.bubble_max)

    ### Per-locus metrics ###
    OutputLocusMetrics(locus_results, args.out, args.noplot)

    ### Per-sample metrics ###
    OutputSampleMetrics(sample_results, samples, args.out, args.noplot)

    return 0


def run():  # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)


if __name__ == "__main__":  # pragma: no cover
    run()
