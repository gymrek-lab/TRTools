import argparse
import os
import pytest
from .mergeSTR import * 

# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = "test"
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
    fname1 = "./test_files/test_file_short1.vcf.gz"
    fname2 = "./test_files/test_file_short2.vcf.gz"
    args.vcfs = fname1 + "," + fname2
    retcode = main(args)
    assert retcode==0
