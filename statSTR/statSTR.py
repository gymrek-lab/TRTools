#!/usr/bin/env python3

"""
Tool for computing stats on STR VCF files
"""
# TODO add arguments for other stats here: mean, mode, variance, nheterozygosity

# Imports
import argparse
import os
import sys
import vcf

# Load local libraries
if __name__ == "statSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import common
    import tr_harmonizer as trh
    import utils
else:
    import strtools.utils.common as common
    import strtools.utils.tr_harmonizer as trh
    import strtools.utils.utils as utils

def GetThresh(trrecord, samplelist=[]):
    return trrecord.GetMaxAllele(self, samplelist=samplelist)

def GetAFreq(trrecord, samplelist=[], count=False, uselength=True):
    r"""Return allele frequency 

    Parameters
    ----------
    trrecord: trh.TRRecord object 
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic 
    count: 
    uselength: bool 
          Whether we should collapse alleles by length

    Returns
    -------

    """
    if count:
        allele_counts = trrecord.GetAlleleCcounts(samplelist=samplelist, uselength=uselength)
        return ",".join(["%s:%i"%(a, acounts.get(a, 0)) for a in sorted(allele_counts.keys())])
    else:
        allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
        return ",".join(["%s:%.3f"%(a, acounts.get(a, 0)*1.0/total) for a in sorted(allele_freqs.keys())])

def GetHWEP(trrecord, samplelist=[], uselength=True):
    r"""Compute Hardy Weinberg p-value 

    Tests whether the number of observed heterozygous vs.
    homozygous individuals is different than expected
    under Hardy Weinberg Equilibrium given the observed
    allele frequencies, based on a binomial test.

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length 

    Returns 
    -------
    p-value: float 
          The two-sided p-value returned by a binomial test (scipy.stats.binom_test) 
          If the allele frequencies dictionary is invalid, return np.nan 
          If the genotype alleles not included in frequencies dictionary, return np.nan

    """
 
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    genotype_counts = trrecord.GetGenotypeCounts(samplelist=samplelist, uselength=uselength)
    return utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)

def GetHet(trrecord, samplelist=[], uselength=True):
    r"""Compute heterozygosity of a locus 

    Heterozygosity is defined as the probability 
    that two randomly drawn allele are different. 

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns 
    -------
    heterozygosity: float
          The heterozygosity of the locus
          If the allele frequencies dictionary is invalid, return np.nan 

    """
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    return utils.GetHeterozygosity(allele_freqs)

def GetMean(trrecord, samplelist=[], uselength=True):
    r"""Compute the mean of the allele frequencies 

    Parameters 
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns 
    -------
    mean: float 
          The mean of the allele frequencies 
          If the allele frequencies dictionary is invalid, return np.nan 

    """
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    return utils.getMean(allele_freqs)

def GetMode(trrecord, samplelist=[], uselength=True):
    r"""Compute the mode of the allele frequencies 

    Parameters
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns 
    -------
    mode: float 
	  The mode of the allele frequencies 
          If the allele frequencies dictionary is invalid, return np.nan 

    """
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    return utils.getMode(allele_freqs)

def GetVariance(trrecord, samplelist=[], uselength=True):
    r"""Compute the variance of the allele frequencies 

    Parameters 
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length 

    Returns 
    -------
    variance: float 
          The variance of the allele frequencies 
          If the allele frequencies dictionary is invalid, return np.nan 

    """
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    return utils.getVariance(allele_freqs)

def GetNumSamples(trrecord, samplelist=[], uselength=True):
    r"""Compute the number of samples 

    Parameters 
    ----------
    trrecord: trh.TRRecord object
          The record that we are computing the statistic for
    samplelist: list of str
          List of the samples that we include when compute the statistic
    uselength: bool
          Whether we should collapse alleles by length

    Returns
    -------
    numSamples: int 
          The number of samples 
          If the allele frequencies dictionary is invalid, return np.nan 

    """ 
    allele_freqs = trrecord.GetAlleleFreqs(samplelist=samplelist, uselength=uselength)
    return utils.getNumSamples(allele_freqs)

def getargs(): 
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Name of output file. Use stdout for standard output.", type=str, required=True)
    inout_group.add_argument("--vcftype", help="Options=%s"%trh.VCFTYPES.__members__, type=str, default="auto")
    filter_group = parser.add_argument_group("Filtering group")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--region", help="Restrict to this region chrom:start-end", type=str)
    stat_group = parser.add_argument_group("Stats group")
    stat_group.add_argument("--thresh", help="Output threshold field (max allele size, used for GangSTR strinfo).", action="store_true")
    stat_group.add_argument("--afreq", help="Output allele frequencies", action="store_true")
    stat_group.add_argument("--acount", help="Output allele counts", action="store_true")
    stat_group.add_argument("--hwep", help="Output HWE p-values per loci.", action="store_true")
    stat_group.add_argument("--het", help="Output observed heterozygote counts used for HWE per loci.", action="store_true")
    stat_group.add_argument("--mean", help="Output mean of allele frequencies.", action="store_true")
    stat_group.add_argument("--mode", help="Output mode of allele frequencies.", action="store_true")
    stat_group.add_argument("--var", help="Output variance of allele frequencies.", action="store_true")
    stat_group.add_argument("--numSamples", help="Output number of samples.", action="store_true")
    stat_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    args = parser.parse_args()
    return args 

def main(args=None):
    if args is None:
        args = getargs()
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    # Load samples
    if args.samples:
        samplelist = [item.strip() for item in open(args.samples, "r").readlines()]
    else: samplelist = []

    invcf = vcf.Reader(filename=args.vcf)
    tr_harmonizer = trh.TRRecordHarmonizer(invcf, vcftype=args.vcftype)

    header = ["chrom","start","end"]
    if args.thresh: header.append("thresh")
    if args.afreq: header.append("afreq")
    if args.acount: header.append("acount")
    if args.hwep: header.append("hwep")
    if args.het: header.append("obs_het")
    if args.mean: header.append("mean")
    if args.mode: header.append("mode")
    if args.var: header.append("var")
    if args.numSamples: header.append("numSamples")
    if args.out == "stdout":
        outf = sys.stdout
    else:
        outf = open(args.out, "w")
    outf.write("\t".join(header)+"\n")

    if args.region: regions = invcf.fetch(args.region)
    else: regions = invcf
    for record in regions:
        trrecord = tr_harmonizer.HarmonizeRecord(record)
        items = [record.CHROM, record.POS, record.INFO["END"]]
        if args.thresh:
            items.append(GetThresh(trrecord, samplelist=samplelist))
        if args.afreq:
            items.append(GetAFreq(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.acount:
            items.append(GetAFreq(trrecord, samplelist=samplelist, uselength=args.uselength, count=True))
        if args.hwep:
            items.append(GetHWEP(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.het:
            items.append(GetHet(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.mean: 
            items.append(GetMean(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.mode:
            items.append(GetMode(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.var:
            items.append(GetVariance(trrecord, samplelist=samplelist, uselength=args.uselength))
        if args.numSamples:
            items.append(GetNumSamples(trrecord, samplelist=samplelist, uselength=args.uselength))
        outf.write("\t".join([str(item) for item in items])+"\n")
    outf.close()
    return 0 

if __name__ == "__main__":
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)
