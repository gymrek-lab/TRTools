"""
Generate the SISTR index
"""

import sys

from . import sistr_utils as sutils

def main(args):
	#### Load configuration ####
	config = sutils.LoadSISTRConfig(args)
	if args.verbose:
		sutils.PrintConfigInfo(config)

	# TODO

	if config is None:
		return 1
	return 0