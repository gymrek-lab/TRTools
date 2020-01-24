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


