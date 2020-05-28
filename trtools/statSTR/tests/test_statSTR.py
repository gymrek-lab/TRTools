import argparse
import os, sys
import pytest

from ..statSTR import *


# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = str(tmpdir /  "test")
    args.vcftype = "auto"
    args.samples = None 
    args.sample_prefixes = None
    args.plot_afreq = False
    args.region = None
    args.thresh = False
    args.afreq = False
    args.acount = False
    args.hwep = False
    args.het = False
    args.use_length = False
    args.mean = False
    args.mode = False 
    args.var = False 
    args.numcalled = False 
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
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

# Test all the statSTR options
def test_Stats(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.thresh = True
    args.het = True
    args.hwep = True
    args.numcalled = True
    args.mode = True
    args.mean = True
    args.var = True
    args.acount = True
    args.afreq = True
    assert main(args)==0
    args.uselength = True
    assert main(args)==0
    args.samples = os.path.join(vcfdir, "test_gangstr_samples.txt") 
    args.out = "stdout"
    args.region = "chr1:3045469-3045470"
    assert main(args)==1

def test_Stats2(args, vcfdir):
    args.vcf = os.path.join(vcfdir, "mergeSTR_vcfs", "test_file_gangstr1.vcf.gz")
    args.region = "chr1:3045469-3045470"
    assert main(args)==0

def test_PlotAfreq(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.plot_afreq = True
    assert main(args)==0
