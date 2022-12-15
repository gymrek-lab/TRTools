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