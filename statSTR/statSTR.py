#!/usr/bin/env python3

"""
Tool for computing stats on STR VCF files
"""
# Example: statSTR --vcf /storage/mgymrek/gangstr-analysis/round2/trio/trio_gangstr_filtered_level1.vcf.gz --out test.tab --thresh
# Example: statSTR --vcf
# TODO add arguments for other stats here: mean, heterozygosity

# Imports
import argparse
import os
import sys
import vcf

# Load local libraries
import dumpSTR.filters as filters
import strtools.utils.common as common
import strtools.utils.utils as utils

"""
Compute the maximum allele length seen
"""
def GetThresh(record, samplelist=[]):
    values = []
    for sample in record:
        if len(samplelist) > 0:
            if sample.sample not in samplelist: continue
        if sample.called:
            values.extend(sample["REPCN"])
    if len(values) == 0: return -1
    return max(values)

def GetAFreq(record, samplelist=[], count=False):
    acounts = {}
    reflen = len(record.REF)
    for sample in record:
        if len(samplelist) > 0:
            if sample.sample not in samplelist: continue
        if sample.called:
            alleles = sample.gt_bases.split(sample.gt_phase_char())
            alens = [len(item)-reflen for item in alleles]
            for a in alens: acounts[a] = acounts.get(a, 0) + 1
    total = sum(acounts.values())
    if len(acounts.keys()) == 0: return "."
    if total == 0: return "."
    if count:
        return ",".join(["%s:%i"%(a, acounts.get(a, 0)) for a in sorted(acounts.keys())])
    afreqs = ",".join(["%s:%.3f"%(a, acounts.get(a, 0)*1.0/total) for a in sorted(acounts.keys())])
    return afreqs

def GetHWEP(record, samplelist=[], use_length=False, het_output=False):
    '''For each STR loci, it outputs the Hardy-Weinberg equilibrium exact test p-value.
    If het=True, then Extended output is HWE P and obs_het.
    If het=False, then Standard output is HWE P.'''
    return utils.GetSTRHWE(record, samples=samplelist, uselength=use_length, het_output=het_output)



def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Name of output file. Use stdout for standard output.", type=str, required=True)
    filter_group = parser.add_argument_group("Filtering group")
    filter_group.add_argument("--samples", help="File containing list of samples to include", type=str)
    filter_group.add_argument("--region", help="Restrict to this region chrom:start-end", type=str)
    stat_group = parser.add_argument_group("Stats group")
    stat_group.add_argument("--thresh", help="Output threshold field (for GangSTR strinfo). Threshold is set to the max observed allele length", action="store_true")
    stat_group.add_argument("--afreq", help="Output allele frequencies", action="store_true")
    stat_group.add_argument("--acount", help="Output allele counts", action="store_true")
    stat_group.add_argument("--hwep", help="Output HWE p-values per loci.", action="store_true")
    stat_group.add_argument("--het", help="Output observed heterozygote counts used for HWE per loci.", action="store_true")
    stat_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")

    args = parser.parse_args()

    # Load samples
    if args.samples:
        samplelist = [item.strip() for item in open(args.samples, "r").readlines()]
    else: samplelist = []

    invcf = vcf.Reader(filename=args.vcf)
    header = ["chrom","start","end"]
    if args.thresh: header.append("thresh")
    if args.afreq: header.append("afreq")
    if args.acount: header.append("acount")
    if args.hwep: header.append("hwep")
    if args.het: header.append("obs_het")
    if args.out == "stdout":
        outf = sys.stdout
    else:
        outf = open(args.out, "w")
    outf.write("\t".join(header)+"\n")

    if args.region: regions = invcf.fetch(args.region)
    else: regions = invcf
    for record in regions:
        items = [record.CHROM, record.POS, record.INFO["END"]]
        if args.thresh:
            items.append(GetThresh(record, samplelist=samplelist))
        if args.afreq:
            items.append(GetAFreq(record, samplelist=samplelist))
        if args.acount:
            items.append(GetAFreq(record, samplelist=samplelist, count=True))
        if args.hwep & ~args.het:
            items.append(GetHWEP(record, samplelist=samplelist, use_length=args.use_length, het_output=args.het))
        if args.hwep & args.het:
            items.append(GetHWEP(record, samplelist=samplelist, use_length=args.use_length, het_output=args.het)[0])
        if args.het:
            items.append(GetHWEP(record, samplelist=samplelist, use_length=args.use_length, het_output=args.het)[1])
        outf.write("\t".join([str(item) for item in items])+"\n")

    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()
