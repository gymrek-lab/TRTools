import argparse
import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','statSTR'))
from statSTR import *

COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
DUMPDIR = os.path.join(COMMDIR, "dump")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")

# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = os.path.join(DUMPDIR, "test")
    args.vcftype = "auto"
    args.samples = None 
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
def test_WrongFile():
    args = base_argparse()
    fname = os.path.join(VCFDIR, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

# Test the right file or directory
def test_RightFile():
    args = base_argparse()
    fname = os.path.join(VCFDIR, "test_gangstr.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

# Test all the statSTR options
def test_Stats():
    args = base_argparse()
    fname = os.path.join(VCFDIR, "test_gangstr.vcf")
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
    args.samples = os.path.join(VCFDIR, "test_gangstr_samples.txt") 
    args.out = "stdout"
    args.region = "chr1:3045469-3045470"
    assert main(args)==1
    args = base_argparse()
    args.vcf = os.path.join(VCFDIR, "mergeSTR_vcfs", "test_file_gangstr1.vcf.gz")
    args.region = "chr1:3045469-3045470"
    assert main(args)==0
