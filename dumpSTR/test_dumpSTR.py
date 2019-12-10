import argparse
import os
import pytest
from .dumpSTR import *

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
    args.min_supp_reads = None
    args.expansion_prob_het = None
    args.expansion_prob_hom = None
    args.expansion_prob_total = None
    args.filter_span_only = None
    args.filter_spanbound_only = None
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
