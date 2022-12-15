"""
Tool for performing LD-clumping on SNPs, STRs, or both
"""

import argparse
import sys

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

from haptools import data

from trtools import __version__

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
    ### Clumping arguments ###
    clump_group = parser.add_argument_group("Clumping Parameters")
    clump_group.add_argument("--clump-p1", help="index variant p-value threshold", type=float, default=0.0001)
    clump_group.add_argument("--clump-p2", help="SP2 column p-value threshold", type=float, default=0.01)
    clump_group.add_argument("--clump-r2", help="r2 treshold", type=float, default=0.5)
    clump_group.add_argument("--clump-kb", help="clump kb radius", type=float, default=250)
    ### Optional args ###
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version='{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    pass # TODO - do all of the things

def run():  # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)


if __name__ == "__main__":  # pragma: no cover
    run()
