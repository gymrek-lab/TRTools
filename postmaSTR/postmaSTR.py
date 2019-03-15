#!/usr/bin/env python3

"""
Tool for post-processing and priotization of STR genotypes

"""

# Imports


# Load local libraries
import strtools.utils.common as common

# Load external libraries
import argparse
import vcf
from vcf.parser import _Filter
from vcf.parser import _Format
from vcf.parser import _Info

def CheckFilters(invcf, args):
    """
    Perform checks on user input for filters
    Input:
    - invcf (vcf.Reader)
    - args (argparse namespace)
    Exit program if checks fail
    """
    if args.affec_expansion_prob_het is not None:
        if args.affec_expansion_prob_het < 0 or args.affec_expansion_prob_het > 1:
            common.ERROR("--affec-expansion-prob-het must be between 0 and 1")
    if args.unaff_expansion_prob_het is not None:
        if args.unaff_expansion_prob_het < 0 or args.unaff_expansion_prob_het > 1:
            common.ERROR("--unaff-expansion-prob-het must be between 0 and 1")
    if args.affec_expansion_prob_hom is not None:
        if args.affec_expansion_prob_hom < 0 or args.affec_expansion_prob_hom > 1:
            common.ERROR("--affec_expansion-prob-hom must be between 0 and 1")
    if args.unaff_expansion_prob_hom is not None:
        if args.unaff_expansion_prob_hom < 0 or args.unaff_expansion_prob_hom > 1:
            common.ERROR("--unaff_expansion-prob-hom must be between 0 and 1")
    if args.affec_expansion_prob_total is not None:
        if args.affec_expansion_prob_total < 0 or args.affec_expansion_prob_total > 1:
            common.ERROR("--affec_expansion-prob-total must be between 0 and 1")
    if args.unaff_expansion_prob_total is not None:
        if args.unaff_expansion_prob_total < 0 or args.unaff_expansion_prob_total > 1:
            common.ERROR("--unaff_expansion-prob-total must be between 0 and 1")

def BuildPostMaSTRCallFilters(args):
    """
    Build list of locus-level filters to include for affected and unaffected samples.
    Input:
    - args (namespace from parser.parse_args)
    
    Output:
    - cdict ({'affec': list<filters.Filter>, 'unaff': list<filters.Filter>}): list of call-level filters
    """
    cdict = {'affec':[], 'unaff':[]}
    if args.affec_expansion_prob_het is not None:
        cdict['affec'].append(filters.ProbHet(args.affec_expansion_prob_het))
    if args.affec_expansion_prob_hom is not None:
        cdict['affec'].append(filters.ProbHom(args.affec_expansion_prob_hom))
    if args.affec_expansion_prob_total is not None:
        cdict['affec'].append(filters.ProbTotal(args.affec_expansion_prob_total))
    if args.unaff_expansion_prob_het is not None:
        cdict['unaff'].append(filters.ProbHet(args.unaff_expansion_prob_het))
    if args.unaff_expansion_prob_hom is not None:
        cdict['unaff'].append(filters.ProbHom(args.unaff_expansion_prob_hom))
    if args.unaff_expansion_prob_total is not None:
        cdict['unaff'].append(filters.ProbTotal(args.unaff_expansion_prob_total))
    return cdict


def ParseFam(filename):
    isAffected = {}
    with open(filename, 'r') as f:
        i = 0
        for line in f:
            i = i + 1
            recs = line.strip().split('\t')
            if len(recs) < 6:
                common.ERROR("Insufficient number of columns in line " + str(i) + " of fam file: " + filename)
            sid = recs[1]
            phe = recs[5]
            isAffected[sid] = phe
    return isAffected
                

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)    
    inout_group.add_argument("--fam", help="Input fam file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)



    #### Affected sample filters
    affec_group = parser.add_argument_group("Call-level filters specific to affected samples")
    affec_group.add_argument("--affec-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    affec_group.add_argument("--affec-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    affec_group.add_argument("--affec-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)
    
    #### Unaffected sample filters
    unaff_group = parser.add_argument_group("Call-level filters specific to unaffected samples")
    unaff_group.add_argument("--unaff-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    unaff_group.add_argument("--unaff-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    unaff_group.add_argument("--unaff-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)

    debug_group = parser.add_argument_group("Debugging parameters")
    debug_group.add_argument("--num-records", help="Only process this many records", type=int)
    debug_group.add_argument("--die-on-warning", help="Quit if a record can't be parsed", action="store_true")
    debug_group.add_argument("--verbose", help="Print out extra info", action="store_true")

    args = parser.parse_args()
    # Load VCF file
    invcf = vcf.Reader(filename=args.vcf)
    # Load FAM file
    isAffected = ParseFam(filename=args.fam)
    
    # Set up filter list
    CheckFilters(invcf, args)
    
    # No need to set up locus level filters.

    # Two sets of parameters set for call filters (for affected and unaffected)
    call_filters = BuildPostMaSTRCallFilters(args)

if __name__ == "__main__":
    main()
