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
import numpy as np
import os
import pandas as pd
import sys
import vcf

# Load local libraries
if __name__ == "qcSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trtools", "utils"))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trtools"))
    import common
    import tr_harmonizer as trh
    import utils
    import version
else: # pragma: no cover
    import trtools.utils.utils as utils
    import trtools.utils.common as common  # pragma: no cover
    import trtools.utils.tr_harmonizer as trh  # pragma: no cover
    import trtools.version as version

__version__ = version.__version__


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

def OutputDiffRefBias(diffs_from_ref, reflens, fname):
    r"""Plot reflen vs. mean difference from ref bias plot

    Parameters
    ----------
    diffs_from_ref : list of int
        Difference of each allele call from the ref allele (in bp)
    reflens : list of int
        List of reference allele lengths for each call (in bp)
    fname : str
        Filename of output plot
    """
    data = pd.DataFrame({"diff": diffs_from_ref, "ref": reflens, "count": [1]*len(reflens)})
    data["ref"] = data["ref"].apply(lambda x: int(x/5)*5) # bin by 5bp
    summ = data.groupby("ref", as_index=False).agg({"diff": np.mean, "count": len}).sort_values("ref") # median or mean?
    summ = summ[summ["count"]>=25] # exclude small counts
    trcounts = np.cumsum(summ["count"])
    trfreqs = trcounts/np.sum(summ["count"])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(summ["ref"], summ["diff"], marker="o", color="darkblue")
    ax.axhline(y=0, linestyle="dashed", color="gray")
    ax.set_xlabel("Reference length (bp)", size=15)
    ax.set_ylabel("Median diff from ref (bp)", size=15)
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
    chroms = sorted(chrom_calls.keys())
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

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf", help="VCF file to analyze.", type=str, required=True)
    req_group.add_argument("--out", help="Output prefix for files generated", type=str, required=True)
    req_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VCFTYPES.__members__], type=str, default="auto")
    filter_group = parser.add_argument_group("Filtering group")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--period", help="Only consider repeats with this motif length", type=int)
    debug_group = parser.add_argument_group("Debug group")
    debug_group.add_argument("--numrecords", help="Only process this many records", type=int)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    # Set up reader and harmonizer
    invcf = utils.LoadSingleReader(args.vcf, checkgz = False)
    if invcf is None:
        return 1
    if args.vcftype != 'auto':
        vcftype = trh.VCFTYPES[args.vcftype]
    else:
        vcftype = trh.InferVCFType(invcf)

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
    diffs_from_ref_unit = [] # for each allele call, keep track of diff (units) from ref
    reflens = [] # for each allele call, keep track of reference length (bp)

    numrecords = 0
    for record in invcf:
        if args.numrecords is not None and numrecords >= args.numrecords: break
        chrom = record.CHROM
        trrecord = trh.HarmonizeRecord(vcftype, record)
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
            diffs_from_ref_unit.extend([allelediff/len(trrecord.motif)]*count)
        numrecords += 1

    OutputDiffRefHistogram(diffs_from_ref_unit, args.out + "-diffref-histogram.pdf")
    OutputDiffRefBias(diffs_from_ref, reflens, args.out + "-diffref-bias.pdf")
    OutputSampleCallrate(sample_calls, args.out+"-sample-callnum.pdf")
    OutputChromCallrate(chrom_calls, args.out+"-chrom-callnum.pdf")
    return 0

def run(): # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
