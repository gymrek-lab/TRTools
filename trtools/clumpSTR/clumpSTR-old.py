"""
Tool for performing LD-clumping on SNPs, STRs, or both
"""

import argparse
import sys

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

from haptools import data

from trtools import __version__

class SummaryStats:
    """
    Keep track of summary statistics

    TODO: add detailed class and methods documentation
    """

    def __init__(self):
        pass # TODO

    def Load(self, statsfile, vartype="SNP", pthresh=1.0):
        """
        Load summary statistics
        Ignore variants with pval < pthresh
        Not yet implemented
        """
        pass # TODO

    def GetNextIndexVariant(self, index_pval_thresh):
        """
        Get the next index variant, which is the 
        variant with the best p-value

        If no more variants below the clump-p1 threshold,
        return None

        Not yet implemented
        """
        return None # TODO

    def QueryWindow(self, indexvar, window_kb):
        """
        Find all candidate variants in the specified
        window around the index variant

        Not yet implemented
        """
        return [] # TODO

    def RemoveClump(self, indexvar, clumpvars):
        """
        Remove the variants from a clump 
        from further consideration
        """
        pass # TODO

def ComputeLD(var1, var2):
    """
    Compute the LD between two variants
    """
    return 0 # TODO

def WriteClump(indexvar, clumped_vars):
    """
    Write a clump to the output file
    Not yet implemented
    """
    pass # TODO

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    ### Input genotypes and summ stats ###
    input_group = parser.add_argument_group("Inputs")
    input_group.add_argument("--summstats-snps", help="SNP summary statistics", type=str)
    input_group.add_argument("--summstats-strs", help="STR summary statistics", type=str)
    input_group.add_argument("--gts-snps", help="SNP genotypes", type=str)
    input_group.add_argument("--gts-strs", help="STR genotypes (VCF format)", type=str)
    ### Output arguments ###
    # TODO: add --out for output prefix
    ### Summary statistics parsing arguments ###
    # TODO: add --clump-snp-field, --clump-field (see plink)
    # TODO: probably also want to add --clump-chr-field --clump-pos-field
    ### Clumping arguments ###
    clump_group = parser.add_argument_group("Clumping Parameters")
    clump_group.add_argument("--clump-p1", help="index variant p-value threshold", type=float, default=0.0001)
    clump_group.add_argument("--clump-p2", help="SP2 column p-value threshold", type=float, default=0.01)
    clump_group.add_argument("--clump-r2", help="r2 threshold", type=float, default=0.5)
    clump_group.add_argument("--clump-kb", help="clump kb radius", type=float, default=250)
    ### Optional args ###
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version='{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    #### Check all input arguments ####
    # TODO make sure sum stats and genotype files exist if specified
    # TODO need one of summstats-snps, summstats-strs
    # TODO if summstats-snps, need gts-snps
    # TODO if summstats-strs, need gts-strs

    #### Set up summary statistics object ####
    summstats = SummaryStats()
    if args.summstats_snps is not None:
        summstats.Load(args.summstats_snps, vartype="SNP", pthresh=args.clump_p2)
    if args.summstats_strs is not None:
        summstats.Load(args.summstats_strs, vartype="STR", pthresh=args.clump_p2)

    #### Perform clumping ####
    indexvar = summstats.GetNextIndexVariant(args.clump_p1)
    while indexvar is not None:
        candidates = summstats.QueryWindow(indexvar, args.clump_kb)
        clumpvars = []
        for c in candidates:
            r2 = ComputeLD(c, indexvar)
            if r2 > args.clump_r2:
                clumpvars.append(c)
        WriteClump(indexvar, clumpvars)
        summstats.RemoveClump(indexvar, clumpvars)
        indexvar = summstats.GetNextIndexVariant(args.clump_p1)

def run():  # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)


if __name__ == "__main__":  # pragma: no cover
    run()
