#!/usr/bin/env python3

"""
Tool for computing stats on STR VCF files
"""
# Example: statSTR --vcf /storage/mgymrek/gangstr-analysis/round2/trio/trio_gangstr_filtered_level1.vcf.gz --out test.tab --thresh
# TODO add arguments for other stats here: mean, heterozygosity, allele freqs

# Imports
import argparse
import os
import sys
import vcf

"""
Compute the maximum allele length seen
"""
def GetThresh(record):
    values = []
    for sample in record:
        if sample.called:
            values.extend(sample["REPCN"])
    if len(values) == 0: return -1
    return max(values)

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Name of output file", type=str, required=True)
    stat_group = parser.add_argument_group("Stats group")
    stat_group.add_argument("--thresh", help="Output threshold field (for GangSTR strinfo). Threshold is set to the max observed allele length", action="store_true")

    args = parser.parse_args()
    invcf = vcf.Reader(filename=args.vcf)
    header = ["chrom","start","end"]
    if args.thresh: header.append("thresh")
    outf = open(args.out, "w")
    outf.write("\t".join(header)+"\n")

    for record in invcf:
        items = [record.CHROM, record.POS, record.INFO["END"]]
        if args.thresh:
            items.append(GetThresh(record))
        outf.write("\t".join([str(item) for item in items])+"\n")

    outf.close()
    sys.exit(0)

if __name__ == "__main__":
    main()
