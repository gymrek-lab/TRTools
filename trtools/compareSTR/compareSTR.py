#!/usr/bin/env python3
"""
Tool for comparing two STR callsets
"""

# Allow making plots even with no x-forward
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# Allow plots to be editable in Adobe Illustrator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Load external libraries
import argparse
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
import pandas as pd
import scipy.stats
import sys
import vcf

# Load local libraries
if __name__ == "compareSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trtools", "utils"))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trtools"))
    import common
    import mergeutils
    import utils
    import tr_harmonizer as trh
    import version
else: # pragma: no cover
    import trtools.utils.common as common  # pragma: no cover
    import trtools.utils.mergeutils as mergeutils  # pragma: no cover
    import trtools.utils.utils as utils  # pragma: no cover
    import trtools.utils.tr_harmonizer as trh # pragma: no cover
    import trtools.version as version

__version__ = version.__version__


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

def OutputLocusMetrics(data, outprefix, noplot):
    r"""Output per-locus metrics

    Outputs text file and plot of per-locus metrics
    outprefix + "-locuscompare.tab"
    outprefix + "-locuscompare.pdf"

    Parameters
    ----------
    data : pd.Dataframe 
        Call comparison results
    outprefix : str
        Prefix to name output file
    noplot : bool
        If True, don't output plots
    """
    # collapse data by locus and output
    perloc = data.groupby(["chrom","start"], as_index=False).agg({"metric-conc-seq": np.mean, "metric-conc-len": np.mean, "sample": len})
    perloc = perloc.sort_values("metric-conc-len", ascending=False)
    perloc.to_csv(outprefix+"-locuscompare.tab", sep="\t", index=False)

    # Create per-locus plot
    if noplot: return
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([item for item in range(perloc.shape[0])], perloc["metric-conc-len"], color="darkblue")
    ax.set_ylabel("Concordance", size=15)
    if perloc.shape[0]<=20:
        ax.set_xticks([item for item in range(perloc.shape[0])])
        ax.set_xticklabels(perloc.apply(lambda x: "%s:%s"%(x["chrom"],x["start"]), 1), size=12, rotation=90)
    else:
        ax.set_xlabel("TR Locus", size=15)
    plt.tight_layout()
    fig.savefig(outprefix+"-locuscompare.pdf")
    plt.close()

def OutputSampleMetrics(data, outprefix, noplot):
    r"""Output per-sample metrics

    Outputs text file and plot of per-sample metrics
    outprefix + "-samplecompare.tab"
    outprefix + "-samplecompare.pdf"

    Parameters
    ----------
    data : pd.Dataframe 
        Call comparison results
    outprefix : str
        Prefix to name output file
    noplot : bool
        If True, don't output plots
    """
    # collapse data by locus and output
    persamp = data.groupby(["sample"], as_index=False).agg({"metric-conc-seq": np.mean, "metric-conc-len": np.mean, "start": len})
    persamp.columns = ["sample","metric-conc-seq","metric-conc-len","numcalls"]
    persamp = persamp.sort_values("metric-conc-len", ascending=False)
    persamp.to_csv(outprefix+"-samplecompare.tab", sep="\t", index=False)

    # Create per-locus plot
    if noplot: return
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter([item for item in range(persamp.shape[0])], persamp["metric-conc-len"], color="darkblue")
    ax.set_ylabel("Concordance", size=15)
    if persamp.shape[0]<=20:
        ax.set_xticks([item for item in range(persamp.shape[0])])
        ax.set_xticklabels(persamp["sample"], size=12, rotation=90)
    else:
        ax.set_xlabel("Sample", size=15)
    plt.tight_layout()
    fig.savefig(outprefix+"-samplecompare.pdf")
    plt.close()

def OutputOverallMetrics(data, format_fields, format_binsizes, stratify_file, period, outprefix):
    r"""Output overall accuracy metrics

    Output metrics overall, by period, and by FORMAT bins
    Output results to outprefix+"-overall.tab"

    Parameters
    ----------
    data : pd.Dataframe 
        Call comparison results
    format_fields : list of str
        List of FORMAT fields to stratify by
    format_binsizes : list of float,float,float
        List of min,max,binsize to stratify formats
    stratify_file : {0, 1, 2}
        Specify whether to apply FORMAT stratification to both files (0), or only (1) or (2)
    period : bool
        If True, also stratify results by period
    outprefix : str
        Prefix to name output file
    """
    f = open(outprefix+"-overall.tab", "w")
    header = ["period"]
    for ff in format_fields: header.append(ff)
    header.extend(["concordance-seq","concordance-len","r2","numcalls"])
    f.write("\t".join(header)+"\n")

    useperiods = ["ALL"]
    if period: useperiods.extend(list(set(data["period"])))
    for per in useperiods:
        if per == "ALL":
            usedata = data
        else: usedata = data[data["period"]==per]
        if usedata.shape[0] < 2: continue
        # Overall
        items = [per] + ["NA"]*len(format_fields)+[np.mean(usedata["metric-conc-seq"]), np.mean(usedata["metric-conc-len"]), \
                                                scipy.stats.pearsonr(usedata["gtsum1"], usedata["gtsum2"])[0],
                                                usedata.shape[0]]
        f.write("\t".join([str(item) for item in items])+"\n")
        
        # By format fields (each separately)
        for i in range(len(format_fields)):
            ff = format_fields[i]
            minval, maxval, binsize = format_binsizes[i]
            bins = np.arange(minval, maxval+binsize, binsize)
            for j in range(len(bins)-1):
                lb = bins[j]
                ub = bins[j+1]
                if stratify_file == 0:
                    ffdata = usedata[(usedata[ff+"1"]>=lb) & (usedata[ff+"1"]<ub) & \
                                     (usedata[ff+"2"]>=lb) & (usedata[ff+"2"]<ub)]
                elif stratify_file == 1:
                    ffdata = usedata[(usedata[ff+"1"]>=lb) & (usedata[ff+"1"]<ub)]
                else:
                    ffdata = usedata[(usedata[ff+"2"]>=lb) & (usedata[ff+"2"]<ub)]
                if ffdata.shape[0] < 2: continue
                ff_vals = ["NA"]*len(format_fields)
                ff_vals[i] = "%s-%s"%(lb,ub)
                items = [per]+ff_vals+[np.mean(ffdata["metric-conc-seq"]), np.mean(ffdata["metric-conc-len"]), \
                                       scipy.stats.pearsonr(ffdata["gtsum1"], ffdata["gtsum2"])[0],
                                       ffdata.shape[0]]
                f.write("\t".join([str(item) for item in items])+"\n")
    f.close()

def GetBubbleLegend(sample_counts):
    r"""Get three good bubble legend sizes to use
    
    They should be nice round numbers spanning the orders of magnitude of the dataset

    Parameters
    ----------
    sample_counts : array-like of int
        Sample counts for each bubble size

    Returns
    -------
    legend_values : list of int
        List of three or fewer representative sample sizes to use for bubble legend
    """
    values = set(sample_counts)
    if len(values) <= 3: return list(values) # if only three values, return three of them
    # Determine if we do log10 or linear scale
    minval = min(values)
    maxval = max(values)
    if maxval/minval > 10:
        # Do log10 scale
        # Find max power of 10
        max10 = int(np.log10(maxval))
        # Find min power of 10
        min10 = int(np.log10(minval))
        # Find power of 10 in between
        mid10 = int((max10+min10)/2)
        return sorted(list(set([10**min10, 10**mid10, 10**max10])))
    else:
        # Do linear scale
        mid = int((minval+maxval)/2)
        return sorted(list(set([minval, mid, maxval])))

def OutputBubblePlot(data, period, outprefix, minval=None, maxval=None):
    r"""Output bubble plot of gtsum1 vs. gtsum2

    Parameters
    ----------
    data : pd.Dataframe 
        Call comparison results
    period : bool
        If True, also stratify results by period
    outprefix : str
        Prefix to name output file
    """
    useperiods = ["ALL"]
    if period: useperiods.extend(list(set(data["period"])))
    for per in useperiods:
        if per == "ALL":
            usedata = data
        else: usedata = data[data["period"]==per]
        bubble_counts = usedata.groupby(["gtsum1","gtsum2"], as_index=False).agg({"sample": len})
        scale = 10000/np.mean(bubble_counts["sample"])
        if minval is None:
            minval = min(min(usedata["gtsum1"]), min(usedata["gtsum2"]))
        if maxval is None:
            maxval = max(max(usedata["gtsum1"]), max(usedata["gtsum2"]))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        # Plot (0,0) separately
        b00 = bubble_counts[(bubble_counts["gtsum1"]==0) & (bubble_counts["gtsum2"]==0)]
        brest = bubble_counts[~((bubble_counts["gtsum1"]==0) & (bubble_counts["gtsum2"]==0))]
        scatter = ax.scatter(b00["gtsum1"], b00["gtsum2"], s=np.sqrt(b00["sample"]*scale), color="darkblue", alpha=0.5)
        scatter = ax.scatter(brest["gtsum1"], brest["gtsum2"], s=np.sqrt(brest["sample"]*scale), color="darkblue", alpha=0.5)
        ax.set_xlabel("GT sum - file 1", size=15)
        ax.set_ylabel("GT sum - file 2", size=15)
        ax.plot([minval, maxval], [minval, maxval], linestyle="dashed", color="gray")
        ax.set_xlim(left=minval, right=maxval)
        ax.set_ylim(bottom=minval, top=maxval)
        ax.axhline(y=0, linestyle="dashed", color="gray")
        ax.axvline(x=0, linestyle="dashed", color="gray")
        # plot dummy points for legend
        legend_values = GetBubbleLegend(bubble_counts["sample"])
        handles = []
        xval = (maxval-minval)/10+minval
        for i in range(len(legend_values)):
            val = legend_values[i]
            step=(maxval-minval)/15
            yval = step*(i+3)
            ax.scatter([xval], [yval], color="darkblue", s=np.sqrt(val*scale))
            ax.annotate(val, xy=(xval+step,yval))
        fig.savefig(outprefix + "-bubble-period%s.pdf"%per)
        plt.close()

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcf1", help="First VCF file to compare (must be sorted, bgzipped, and indexed)", type=str, required=True)
    req_group.add_argument("--vcf2", help="Second VCF file to compare (must be sorted, bgzipped, and indexed)", type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    req_group.add_argument("--vcftype1", help="Type of --vcf1. Options=%s"%[str(item) for item in trh.VCFTYPES.__members__], type=str, default="auto")
    req_group.add_argument("--vcftype2", help="Type of --vcf2. Options=%s"%[str(item) for item in trh.VCFTYPES.__members__], type=str, default="auto")
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
    ### Plotting args ###
    plot_group = parser.add_argument_group("Plotting options")
    plot_group.add_argument("--bubble-min", help="Minimum x/y axis value to display on bubble plots", type=int)
    plot_group.add_argument("--bubble-max", help="Maximum x/y axis value to display on bubble plots", type=int)
    ### Optional args ###
    option_group = parser.add_argument_group("Optional arguments")
    option_group.add_argument("--verbose", help="Print helpful debugging info", action="store_true")
    option_group.add_argument("--numrecords", help="For debugging, only process this many records", type=int)
    option_group.add_argument("--noplot", help="Don't output any plots. Only produce text output", action="store_true")
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
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
    period = len(record1.motif)
    reflen = len(record1.ref_allele)/period
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
        results_dir["period"].append(period)
        results_dir["sample"].append(sample)
        results_dir["gtstring1"].append(",".join(gt_string_1))
        results_dir["gtstring2"].append(",".join(gt_string_2))
        results_dir["gtsum1"].append((sum(gt_len_1)-reflen*2))
        results_dir["gtsum2"].append((sum(gt_len_2)-reflen*2))
        results_dir["metric-conc-seq"].append(int(all([(gt_string_1[i]==gt_string_2[i]) for i in range(len(gt_string_1))])))
        results_dir["metric-conc-len"].append(int(all([(gt_len_1[i]==gt_len_2[i]) for i in range(len(gt_len_1))])))
        for ff in format_fields:
            try:
                val1 = call1[ff]
            except: val1 = np.nan
            try:
                val2 = call2[ff]
            except: val2 = np.nan
            results_dir[ff+"1"].append(val1)
            results_dir[ff+"2"].append(val2)

def main(args):
    ### Check and load VCF files ###
    vcfreaders = utils.LoadReaders([args.vcf1, args.vcf2], checkgz=True, region=args.region)
    if vcfreaders is None or len(vcfreaders) != 2:
        return 1
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)
    
    ### Load shared samples ###
    samples = mergeutils.GetSharedSamples(vcfreaders)
    if len(samples) == 0:
        common.WARNING("No shared smaples found between vcf readers")
        return 1
    if args.samples:
        usesamples = set([item.strip() for item in open(args.samples, "r").readlines()])
        samples = list(set(samples).intersection(usesamples))
    if len(samples) == 0:
        common.WARNING("No shared samples found between files")
        return 1
    
    ### Determine FORMAT fields we should look for ###
    if args.stratify_file is not None and args.stratify_file not in [0,1,2]:
        common.MSG("--stratify-file must be 0,1, or 2")
        return 1
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

    vcftype1 = trh.GetVCFType(vcfreaders[0], args.vcftype1)
    vcftype2 = trh.GetVCFType(vcfreaders[1], args.vcftype2)

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = mergeutils.GetMinRecords(current_records, chroms)

    done = mergeutils.DoneReading(current_records)
    num_records = 0
    while not done:
        if any([item is None for item in current_records]): break
        if args.numrecords is not None and num_records >= args.numrecords: break
        if args.verbose: mergeutils.PrintCurrentRecords(current_records, is_min)
        if mergeutils.CheckMin(is_min): return 1
        if all([is_min]):
            if (current_records[0].CHROM == current_records[1].CHROM and \
                current_records[0].POS == current_records[1].POS):
                UpdateComparisonResults(trh.HarmonizeRecord(vcftype1, current_records[0]), \
                                        trh.HarmonizeRecord(vcftype2, current_records[1]), \
                                        format_fields, samples, results_dir)
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        is_min = mergeutils.GetMinRecords(current_records, chroms)
        done = mergeutils.DoneReading(current_records)
        num_records += 1

    ### Load all results to a dataframe and output full results ###
    data = pd.DataFrame(results_dir)
    data.to_csv(args.out + "-callcompare.tab", sep="\t", index=False)

    ### Overall metrics ###
    OutputOverallMetrics(data, format_fields, format_binsizes, args.stratify_file, args.period, args.out)
    if not args.noplot: OutputBubblePlot(data, args.period, args.out, minval=args.bubble_min, maxval=args.bubble_max)

    ### Per-locus metrics ###
    OutputLocusMetrics(data, args.out, args.noplot)

    ### Per-sample metrics ###
    OutputSampleMetrics(data, args.out, args.noplot)

    return 0 

def run(): # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()

