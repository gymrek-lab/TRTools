import argparse
import os
import pytest
from .statSTR import *

# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = "test"
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
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname1 = "xxxx"
    if os.path.exists(fname1):
        os.remove(fname1)
    args.vcf = fname1
    retcode = main(args)
    assert retcode==1

# Test the right file or directory
def test_RightFile():
    args = base_argparse()
    fname1 = "/home/npusarla/workspace/str/STRTools/statSTR/test_files/test_file.vcf"
    args.vcf = fname1
    retcode = main(args)
    assert retcode==0


