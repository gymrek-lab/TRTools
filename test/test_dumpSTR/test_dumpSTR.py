import argparse
import os, sys
import pytest
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','dumpSTR'))
from dumpSTR import *

TESTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files")

# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcf = None
    args.out = "test"
    args.min_locus_callrate = None
    args.min_locus_hwep = None
    args.min_locus_het = None
    args.max_locus_het = None
    args.use_length = False
    args.filter_regions = None
    args.filter_regions_names = None
    args.filter_hrun = False
    args.drop_filtered = False
    args.min_call_DP = None
    args.max_call_DP = None
    args.min_call_Q = None
    args.max_call_flank_indel = None
    args.max_call_stutter = None
    args.min_supp_reads = 0
    args.expansion_prob_het = None
    args.expansion_prob_hom = None
    args.expansion_prob_total = None
    args.filter_span_only = False
    args.filter_spanbound_only = False 
    args.filter_badCI = None
    args.require_support = 0
    args.readlen = None
    args.num_records = None
    args.die_on_warning = False
    args.verbose = False
    return args

# Test no such file or directory
def test_WrongFile():
    args = base_argparse()
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname = "xxxx"
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1


#Test if test works for the right filename 
def test_RightFile():
    args = base_argparse()
    # Try an actual file 
    fname = os.path.join(TESTDIR, "test_file.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

def test_Filters(): 
    args = base_argparse() 
    fname = os.path.join(TESTDIR, "test_file.vcf")
    args.vcf = fname
    invcf = vcf.Reader(filename=args.vcf)
    retcode = main(args)
    assert retcode==0
