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
    args.quiet = False
    args.verbose = False
    args.vcftype = "auto"
    return args

# Test no such file or directory
def test_WrongFile():
    args = base_argparse()
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname1 = os.path.join(VCFDIR, "test_non_existent1.vcf")
    fname2 = os.path.join(VCFDIR, "test_non_existent2.vcf")
    if os.path.exists(fname1):
        os.remove(fname1)
    if os.path.exists(fname2):
        os.remove(fname2)
    args.vcfs = fname1 + "," + fname2
    retcode = main(args)
    assert retcode==1

# Test right files or directory - GangSTR
def test_GangSTRRightFile():
    args = base_argparse()
    fname1 = os.path.join(VCFDIR, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(VCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcftype = "gangstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

def test_GangSTRDuplicate():
    args = base_argparse()
    fname1 = os.path.join(VCFDIR, "test_file_gangstr1.vcf.gz")
    args.vcfs = fname1 + "," + fname1
    assert main(args)==1

# TODO - test EH, advntr, hipSTR, popstr
# TODO - test unzipped, unindexed VCFs return 1
# TODO - test VCFs with different ref genome contigs return 1
