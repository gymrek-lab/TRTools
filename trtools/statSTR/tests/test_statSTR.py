import argparse
import filecmp
import os

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
    args.entropy = False
    args.precision = 2
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
    fname = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.vcf = fname
    args.thresh = True
    args.het = True
    args.entropy = True
    args.hwep = True
    args.numcalled = True
    args.mode = True
    args.mean = True
    args.var = True
    args.acount = True
    args.afreq = True
    assert main(args) == 0
    args.uselength = True
    assert main(args) == 0
    args.region = "chr1:3045469-3045470"
    assert main(args) == 0
    args.samples = os.path.join(vcfdir, "fewer_samples.txt")
    assert main(args) == 0

def test_require_tabix_for_regions(args, vcfdir, capsys):
    # args.samples = os.path.join(vcfdir, "test_gangstr_samples.txt")
    args.region = "chr1:3045469-3045470"
    args.vcf = os.path.join(vcfdir, "test_gangstr.vcf")
    args.thresh = True
    assert main(args) == 1
    assert 'bgzipped' in capsys.readouterr().err

def test_Stats2(args, vcfdir):
    args.vcf = os.path.join(vcfdir, "mergeSTR_vcfs", "test_file_gangstr1.vcf.gz")
    args.region = "chr1:3045469-3045470"
    assert main(args)==0

def test_PlotAfreq(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.plot_afreq = True
    assert main(args)==0

def test_stats_output(args, vcfdir, statsdir):
    """
    Run statSTR on a file which statSTR has been run on in the past
    and confirm the results haven't changed
    """
    fname = os.path.join(vcfdir, "many_samples.vcf.gz")
    args.vcf = fname
    args.thresh = True
    args.afreq = True
    args.acount = True
    args.hwep = True
    args.het = True
    args.entropy = True
    args.mean = True
    args.mode = True
    args.var = True
    args.numcalled = True
    # exclude an allele which doesn't have
    # reproducible stats even up to two decimal places
    args.region = "1:1-3683401"
    assert main(args) == 0
    assert filecmp.cmp(
        args.out + ".tab",
        os.path.join(statsdir, "many_samples_all.tab")
    )

