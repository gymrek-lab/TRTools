#!/usr/bin/env python3
"""
Tool for computing stats on STR VCF files
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
import numpy as np
import os
import sys
import vcf

MAXPLOTS = 10 # don't plot more than this many allele freqs

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from trtools import __version__


def PlotAlleleFreqs(trrecord, outprefix, samplelists=None, sampleprefixes=None):
    r"""Plot allele frequencies for a locus

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    outprefix : str
          Prefix for output file
    samplelists: list of list of str, optional
          List of lists of the samples that we include when compute the statistic
    sampleprefixes : list of str, optional
          Prefixes for each sample list to use in legend
    """
    if samplelists is None or samplelists == []:
        samplelists = [None]
        sampleprefixes = ["sample"]
    allele_freqs_list = []
    allele_set = set()
    for sl in samplelists:
        afreqs = trrecord.GetAlleleFreqs(uselength=True, samplelist=sl)
        allele_freqs_list.append(afreqs)
        allele_set = allele_set.union(afreqs.keys())
    min_allele = min(allele_set)-2
    max_allele = max(allele_set)+2
    bins = np.arange(min_allele, max_allele, 1)

    fname = outprefix + "-%s-%s.pdf"%(trrecord.vcfrecord.CHROM, trrecord.vcfrecord.POS)
    w = 1.0/(len(samplelists)+0.3)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(samplelists)):
        ax.bar([item+i*w for item in bins], [allele_freqs_list[i].get(item, 0) for item in bins],
               label=sampleprefixes[i], width=w*1.1)
    ax.legend()
    ax.set_xlabel("TR allele (num. %s rpts)"%trrecord.motif, size=15)
    ax.set_ylabel("Frequency", size=15)
    ax.set_xticklabels([int(item) for item in ax.get_xticks()], size=12)
    ax.set_yticklabels(["%.2f"%item for item in ax.get_yticks()], size=12)
    fig.tight_layout()
    fig.savefig(fname)
    plt.close()

def GetHeader(header, sample_prefixes):
    r"""Return header items for a column

    Parameters
    ----------
    header : str
       Header item
    sample_prefixes : list of str
       List of sample prefixes. empty if no sample groups used

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelists: list of list of str
          List of lists of the samples that we include when compute the statistic

    Returns
    -------
    thresh: list of float
          List of Maximum allele length observed in each sample group
    """
    if len(samplelists) == 0:
        return [trrecord.GetMaxAllele()]
    else:
        return [trrecord.GetMaxAllele(samplelist=sl) for sl in samplelists]

def GetAFreq(trrecord, samplelists=[], count=False, uselength=True):
    r"""Return allele frequency for a TR

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of lists of the samples that we include when compute the statistic
    count: bool
          If True, return allele counts rather than allele frequencies
    uselength: bool
          Whether we should collapse alleles by length

    Returns
    -------
    allele_freqs_strs: list of str
          Format: allele1:freq1,allele2:freq2,etc. for each sample group
    """
    if len(samplelists) == 0:
        samplelists.append(None)
    allele_freqs_strs = []
    for sl in samplelists:
        if count:
            allele_counts = trrecord.GetAlleleCounts(uselength=uselength, samplelist=sl)
            allele_freqs_strs.append(",".join(["%s:%i"%(a, allele_counts.get(a, 0)) for a in sorted(allele_counts.keys())]))
        else:
            allele_freqs = trrecord.GetAlleleFreqs(uselength=uselength, samplelist=sl)
            allele_freqs_strs.append(",".join(["%s:%.3f"%(a, allele_freqs.get(a, 0)) for a in sorted(allele_freqs.keys())]))
    return allele_freqs_strs

def GetHWEP(trrecord, samplelists=[], uselength=True):
    r"""Compute Hardy Weinberg p-value

    Tests whether the number of observed heterozygous vs.
    homozygous individuals is different than expected
    under Hardy Weinberg Equilibrium given the observed
    allele frequencies, based on a binomial test.

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of list of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns
    -------
    p-value: list of float
          The two-sided p-value returned by a binomial test (scipy.stats.binom_test)
          If the allele frequencies dictionary is invalid, return np.nan
          If the genotype alleles not included in frequencies dictionary, return np.nan
          One value returned for each samplelist
    """
    if len(samplelists)==0: samplelists.append(None)
    pvals = []
    for sl in samplelists:
        allele_freqs = trrecord.GetAlleleFreqs(samplelist=sl, uselength=uselength)
        genotype_counts = trrecord.GetGenotypeCounts(samplelist=sl, uselength=uselength)
        pvals.append(utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts))
    return pvals

def GetHet(trrecord, samplelists=[], uselength=True):
    r"""Compute heterozygosity of a locus

    Heterozygosity is defined as the probability
    that two randomly drawn allele are different.

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of list of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns
    -------
    heterozygosity: list of float
          The heterozygosity of the locus. One value for each sample list.
          If the allele frequencies dictionary is invalid, return np.nan
    """
    if len(samplelists) == 0: samplelists.append(None)
    hetvals = []
    for sl in samplelists:
        allele_freqs = trrecord.GetAlleleFreqs(samplelist=sl, uselength=uselength)
        hetvals.append(utils.GetHeterozygosity(allele_freqs))
    return hetvals

def GetMean(trrecord, samplelists=[], uselength=True):
    r"""Compute the mean allele length

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of list of the samples that we include when compute the statistic

    Returns
    -------
    mean: list of float
          The mean allele length. One value for each sample list
          If the allele frequencies dictionary is invalid, return np.nan
    """
    if len(samplelists) == 0: samplelists.append(None)
    return [utils.GetMean(trrecord.GetAlleleFreqs(samplelist=sl, uselength=True)) for sl in samplelists]

def GetMode(trrecord, samplelists=[]):
    r"""Compute the mode of the allele frequencies

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of the samples that we include when compute the statistic

    Returns
    -------
    mode: list of float
	  The mode of the allele frequencies. One value for each sample list
          If the allele frequencies dictionary is invalid, return np.nan
    """
    if len(samplelists) == 0: samplelists.append(None)
    return [utils.GetMode(trrecord.GetAlleleFreqs(samplelist=sl, uselength=True)) for sl in samplelists]

def GetVariance(trrecord, samplelists=[]):
    r"""Compute the variance of the allele lengths

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelists: list of list of str
          List of list of the samples that we include when compute the statistic

    Returns
    -------
    variance: list of float
          The variance of the allele lengths. One value for each sample list
          If the allele frequencies dictionary is invalid, return np.nan
    """
    if len(samplelists) == 0: samplelists.append(None)
    return [utils.GetVariance(trrecord.GetAlleleFreqs(samplelist=sl, uselength=True)) for sl in samplelists]

def GetNumSamples(trrecord, samplelists=[]):
    r"""Compute the number of samples

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of list of str
          List of list of the samples that we include when compute the statistic

    Returns
    -------
    numSamples: list of int
          The number of samples. One value for each sample list
          If the allele frequencies dictionary is invalid, return np.nan
    """
    if len(samplelists) == 0: samplelists.append(None)
    return [sum(trrecord.GetGenotypeCounts(samplelist=sl).values()) for sl in samplelists]

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Output file prefix. Use stdout to print file to standard output.", type=str, required=True)
    inout_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VcfTypes.__members__], type=str, default="auto")
    filter_group = parser.add_argument_group("Filtering group")
    filter_group.add_argument("--samples", help="File containing list of samples to include. Or a comma-separated list of files to compute stats separate for each group of samples", type=str)
    filter_group.add_argument("--sample-prefixes", help="Prefixes to name output for each samples group. By default uses 1,2,3 etc.", type=str)
    filter_group.add_argument("--region", help="Restrict to this region chrom:start-end", type=str)
    stat_group_name = "Stats group"
    stat_group = parser.add_argument_group(stat_group_name)
    stat_group.add_argument("--thresh", help="Output threshold field (max allele size, used for GangSTR strinfo).", action="store_true")
    stat_group.add_argument("--afreq", help="Output allele frequencies", action="store_true")
    stat_group.add_argument("--acount", help="Output allele counts", action="store_true")
    stat_group.add_argument("--hwep", help="Output HWE p-values per loci.", action="store_true")
    stat_group.add_argument("--het", help="Output heterozygosity of each locus.", action="store_true")
    stat_group.add_argument("--mean", help="Output mean of allele frequencies.", action="store_true")
    stat_group.add_argument("--mode", help="Output mode of allele frequencies.", action="store_true")
    stat_group.add_argument("--var", help="Output variance of allele frequencies.", action="store_true")
    stat_group.add_argument("--numcalled", help="Output number of samples called.", action="store_true")
    stat_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    plot_group = parser.add_argument_group("Plotting group")
    plot_group.add_argument("--plot-afreq", help="Output allele frequency plot. Will only do for a maximum of 10 TRs.", action="store_true")
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    # If no stat selected, print an error message and terminate
    stat_dict = {}
    for grp in parser._action_groups:
        if grp.title == stat_group_name:
            stat_dict = {a.dest:getattr(args,a.dest,None) for a in grp._group_actions}

    if not any(stat_dict.values()):
        common.WARNING("Error: Please use at least one of the flags in the Stats group. See statSTR --help for options.")
        return None
    return args

def main(args):
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    # Load samples
    sample_lists = []
    sample_prefixes = []
    if args.samples:
        sfiles = args.samples.split(",")
        if args.sample_prefixes:
            sample_prefixes = args.sample_prefixes.split(",")
        else:
            sample_prefixes = [str(item) for item in range(1, len(sfiles)+1)]
        if len(sfiles) != len(sample_prefixes):
            common.MSG("--sample-prefixes must be same length as --samples")
            return 1
        for sf in sfiles:
            sample_lists.append([item.strip() for item in open(sf, "r").readlines()])

    invcf = utils.LoadSingleReader(args.vcf, checkgz = False)
    if invcf is None:
        return 1
    if args.vcftype != 'auto':
        vcftype = trh.VcfTypes[args.vcftype]
    else:
        vcftype = trh.InferVCFType(invcf)

    header = ["chrom","start","end"]
    if args.thresh: header.extend(GetHeader("thresh", sample_prefixes))
    if args.afreq: header.extend(GetHeader("afreq", sample_prefixes))
    if args.acount: header.extend(GetHeader("acount", sample_prefixes))
    if args.hwep: header.extend(GetHeader("hwep", sample_prefixes))
    if args.het: header.extend(GetHeader("het", sample_prefixes))
    if args.mean: header.extend(GetHeader("mean", sample_prefixes))
    if args.mode: header.extend(GetHeader("mode", sample_prefixes))
    if args.var: header.extend(GetHeader("var", sample_prefixes))
    if args.numcalled: header.extend(GetHeader("numcalled", sample_prefixes))
    if args.out == "stdout":
        if args.plot_afreq:
            common.MSG("Cannot use --out stdout when generating plots")
            return 1
        outf = sys.stdout
    else:
        outf = open(args.out + ".tab", "w")
    outf.write("\t".join(header)+"\n")

    if args.region:
        if not os.path.isfile(args.vcf+".tbi"):
            common.MSG("Make sure %s is bgzipped and indexed"%args.vcf)
            return 1
        regions = invcf.fetch(args.region)
    else: regions = invcf
    num_plotted = 0
    for record in regions:
        trrecord = trh.HarmonizeRecord(vcftype, record)
        if args.plot_afreq and num_plotted <= MAXPLOTS:
            PlotAlleleFreqs(trrecord, args.out, samplelists=sample_lists, sampleprefixes=sample_prefixes)
            num_plotted += 1
        items = [record.CHROM, record.POS, record.INFO["END"]]
        if args.thresh:
            items.extend(GetThresh(trrecord, samplelists=sample_lists))
        if args.afreq:
            items.extend(GetAFreq(trrecord, samplelists=sample_lists, uselength=args.use_length))
        if args.acount:
            items.extend(GetAFreq(trrecord, samplelists=sample_lists, uselength=args.use_length, count=True))
        if args.hwep:
            items.extend(GetHWEP(trrecord, samplelists=sample_lists, uselength=args.use_length))
        if args.het:
            items.extend(GetHet(trrecord, samplelists=sample_lists, uselength=args.use_length))
        if args.mean:
            items.extend(GetMean(trrecord, samplelists=sample_lists))
        if args.mode:
            items.extend(GetMode(trrecord, samplelists=sample_lists))
        if args.var:
            items.extend(GetVariance(trrecord, samplelists=sample_lists))
        if args.numcalled:
            items.extend(GetNumSamples(trrecord, samplelists=sample_lists))
        outf.write("\t".join([str(item) for item in items])+"\n")
    outf.close()
    return 0

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
