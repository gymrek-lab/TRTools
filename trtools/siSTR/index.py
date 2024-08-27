"""
Generate the SISTR index
"""

import sys

from . import sistr_utils as sutils

def main(args):
	config = sutils.LoadSISTRConfig(args)
	print(config)

	# TODO
	
	if config is None:
		return 1
	return 0