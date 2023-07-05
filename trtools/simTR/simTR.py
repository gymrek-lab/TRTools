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
from trtools import __version__

def main(args):
	if not os.path.exists(args.coords):
		common.WARNING("Error: {} does not exist".format(args.coords))
		return 1
	if not os.path.exists(args.ref):
		common.WARNING("Error: {} does not exist".format(args.ref))
		return 1
	if args.art is not None:
		if not os.path.exists(args.art):
			common.WARNING("Error: ART path {} does not exist".format(args.art))
			return 1
	else:
		if shutil.which("art_illumina") is None:
			common.WARNING("Error: Could not find art_illumina executable")
			return 1
	if args.u < 0 or args.u > 1:
		common.WARNING("Error: --u u ({}) is not between 0 and 1".format(args.u))
		return 1
	if args.d < 0 or args.d > 1:
		common.WARNING("Error: --d ({}) is not between 0 and 1".format(args.d))
		return 1
	if args.rho < 0 or args.rho > 1:
		common.WARNING("Error: --rho ({}}) is not between 0 and 1".format(args.rho))
		return 1
	if args.p_thresh < 0 or args.p_thresh > 1:
		common.WARNING("Error: --p_thresh ({}) is not between 0 and 1".format(args.p_thresh))
		return 1
	if args.coverage < 0:
		common.WARNING("Error: --coverage ({}}) cannot be less than 0".format(args.coverage))
		return 1
	if args.read_length < 0:
		common.WARNING("Error: --read_length ({}) cannot be less than 0".format(args.read_length))
		return 1
	if args.insert < 0:
		common.WARNING("Error: --insert ({}) cannot be less than 0".format(args.insert))
		return 1
	if args.sd < 0:
		common.WARNING("Error: --sd ({}) cannot be less than 0".format(args.sd))
		return 1
	if args.window < 0:
		common.WARNING("Error: --window ({}) cannot be less than 0".format(args.window))
		return 1
	if not os.path.exists(os.path.dirname(os.path.abspath(args.output_dir))):
		common.WARNING("Error: The directory which contains the output location {} does"
			" not exist".format(args.output_dir))
		return 1

	# TODO

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