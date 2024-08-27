"""
Tool for performing selection inference at STRs
"""

import argparse
import enum
import sys

import trtools.utils.common as common
import trtools.utils.utils as utils
from trtools import __version__

from . import index as index
from . import sistr_utils as sutils

class SISTRCommands(enum.Enum):
    """Possible SISTR commands to run."""
    index = "index"
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("command", help="Options=%s"%[str(item) for item in SISTRCommands.__members__])
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)
    index_group = parser.add_argument_group("Indexing options")
    index_group.add_argument(
        "--config",
        help="JSON file with indexing options. Index command line arguments "
             "override individual parameter values. See example: trtools/siSTR/config.json",
        type=str
    )
    index_group.add_argument(
        "--periods",
        help="Comma-separated list of rpt. unit lengths (in bp) to include. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["periods"]])),
        type=str,
    )
    index_group.add_argument(
        "--opt-allele-ranges",
        help="Comma-separated list of optimal allele ranges to include for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item[0])+"-"+str(item[1]) for item in sutils.DEFAULTS["opt_allele_ranges"]])),
        type=str
    )
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    try:
        command = SISTRCommands[args.command]
    except KeyError:
        common.WARNING("Command {cmd} invalid".format(cmd=args.command))
        sys.exit(1)

    if command == SISTRCommands.index:
        return index.main(args)

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()

