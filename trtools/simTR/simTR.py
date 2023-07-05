#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
from pyfaidx import Fasta

def main(args):
	pass # TODO
	
def getargs():
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    # TODO

def run():  # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()