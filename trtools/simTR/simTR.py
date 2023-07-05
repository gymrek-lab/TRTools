#!/usr/bin/env python3
"""
Performs ART simulation including the
characteristic stutter error to STRs
"""

import argparse
import os
import shutil
import sys
from pyfaidx import Fasta
import trtools.utils.utils as utils

def main(args):
	pass # TODO

def getargs():
	parser = argparse.ArgumentParser(
		__doc__,
		formatter_class=utils.ArgumentDefaultsHelpFormatter
	)
	inout_group = parser.add_argument_group("Input/output")
	inout_group.add_argument("--ref", help="Path to reference genome", type=str, required=True)
	inout_group.add_argument("--coords", help="Path to file containing coordinates for the target TR", \
			type=str, required=True)
	inout_group.add_argument("--output_dir", help="Name of output directory", type=str, default="test_dir")
	stutter_group = parser.add_argument_group("Stutter simulation parameters")
	stutter_group.add_argument("--u", help="Probability of adding additional copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--d", help="Probability of deleting copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--rho", help="Size of stutter-induced changes", type=float, default=0.9)
	stutter_group.add_argument("--p_thresh", help="Ignore stutter alleles expected to have lower than this frequency", \
		type=float, default=0.01)
	seq_group = parser.add_argument_group("Sequencing parameters")
	seq_group.add_argument("--coverage", help="Target coverage level", type=int, default=1000)
	seq_group.add_argument("--read_length", help="Length of each read (bp)", type=int, default=100)
	seq_group.add_argument("--insert", help="Mean fragment length", type=int, default=350)
	seq_group.add_argument("--sd", help="Std. deviation of fragmen tlength", type=int, default=50)
	seq_group.add_argument("--window", help="Size of window around target TR to sequence (bp)", type=int, default=1000)
	other_group = parser.add_argument_group("Other options")
	other_group.add_argument("--art", help="Path to ART simulator package", type=str, required=True)
	ver_group = parser.add_argument_group("Version")
	ver_group.add_argument("--version", action="version", \
			version='{version}'.format(version=__version__))
	args = parser.parse_args()
	return args

def run():  # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()