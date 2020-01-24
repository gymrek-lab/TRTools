import argparse
import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','mergeSTR'))
from mergeSTR import * 

TESTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files")
COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
DUMPDIR = os.path.join(COMMDIR, "dump")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")
# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = os.path.join(DUMPDIR, "test")
    args.update_sample_from_file = False 
    args.merge_ggl = False
    args.quiet = False
    args.verbose = False
    return args

# Test no such file or directory
def test_WrongFile():
    args = base_argparse()
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname1 = "xxxx"
    fname2 = "yyyy"
    if os.path.exists(fname1):
        os.remove(fname1)
    if os.path.exists(fname2):
        os.remove(fname2)
    args.vcfs = fname1 + "," + fname2
    retcode = main(args)
    assert retcode==1

# Test right files or directory
def test_RightFile():
    args = base_argparse()
    fname1 = os.path.join(TESTDIR, "test_file.vcf.gz")
    fname2 = os.path.join(TESTDIR, "test_file2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    retcode = main(args)
    assert retcode==0
