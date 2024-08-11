import argparse
import gzip
import os

import pytest

from ..annotaTR import *

# Set up base argparser
@pytest.fixture(name='args')
def args_fixture(tmpdir):
    return args(tmpdir)

def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.vcftype = "auto"
    args.out = str(tmpdir / "test")
    args.outtype = ["vcf"]
    args.dosages = None
    args.ref_panel = None
    return args

@pytest.fixture
def testBeagleDir(vcfdir):
    return vcfdir + "/beagle"

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1
    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.ref_panel = fname
    retcode = main(args)
    assert retcode==1

def test_WrongOutDir(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.out = os.path.join("fakedir","test")
    retcode = main(args)
    assert retcode==1
    args.out = vcfdir + "/"
    retcode = main(args)
    assert retcode==1
   
def test_OutTypes(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.outtype = ["vcf", "pgen"]
    args.dosages = "bestguess_norm"
    retcode = main(args)
    assert retcode==0
    args.outtype = ["pgen", "vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["pgen"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["vcf", "vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["dummy"]
    retcode = main(args)
    assert retcode==1
    
   