import argparse
import os

import pytest

from ..simTR import *

# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.ref = None
    args.repeat_unit = None
    args.outprefix = None
    args.tmpdir = None
    args.u = None
    args.d = None
    args.rho = None
    args.p_thresh = None
    args.coverage = None
    args.read_length = None
    args.insert = None
    args.sd = None
    args.window = None
    args.art = None
    args.single = False
    args.debug = False
    return args

# TODO add tests