import argparse
import os

import pytest

from ..prancSTR import *

# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.out = str(tmpdir / "test")
    args.region = None
    args.only_passing = False
    args.debug = False
    args.vcftype = "hipstr"
    return args

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

# Test the right file or directory
def test_RightFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0