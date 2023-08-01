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

def test_ParseCoordinates():
    chrom, start, end = ParseCoordinates("chr20:40-100")
    assert chrom == "chr20"
    assert start == 40
    assert end == 100

def test_GetMaxDelta():
    pthresh = 0.1
    delta = GetMaxDelta(0.2, 0.5, pthresh)
    assert delta >= pthresh

# TODO add tests