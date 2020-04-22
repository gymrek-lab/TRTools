import argparse
import os, sys
import numpy as np
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..'))
from qcSTR import * 

TESTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files")
COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
DUMPDIR = os.path.join(COMMDIR, "dump")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")
QCDIR = os.path.join(VCFDIR, "qc_vcfs")
# Set up base argparser
def base_argparse():
    args = argparse.ArgumentParser()
    args.vcf = None
    args.out = os.path.join(DUMPDIR, "test_qc")
    args.vcftype = "auto"
    args.samples = None
    args.numrecords = None
    args.period = None
    return args

def test_OutputDiffRefHistogram():
    diffs_from_ref = [0,0,0,0,1,0,-1,-2,-4,-5]
    fname = os.path.join(DUMPDIR, "test_hist.pdf")
    try:
        OutputDiffRefHistogram(diffs_from_ref, fname)
    except:
        pytest.fail("Unexpected error")

def test_OutputDiffRefBias():
    diffs_from_ref = [0,0,0,0,1,0,-1,-2,-4,-5]
    reflens = [1,2,3,4,5,6,7,8,9,10]
    fname = os.path.join(DUMPDIR, "test_bias.pdf") # will be empty because of low count <25
    try:
        OutputDiffRefBias(diffs_from_ref, reflens, fname)
    except:
        pytest.fail("Unexpected error")

def test_OutPutSampleCallrate():
    sample_calls = {'s1': 120, 's2': 10}
    fname = os.path.join(DUMPDIR, "test_qc1.pdf")
    try:
        OutputSampleCallrate(sample_calls, fname)
    except:
        pytest.fail("Unexpected error")

def test_OutPutChromCallrate():
    chrom_calls = {'chr1': 100, 'chr2': 200}
    fname = os.path.join(DUMPDIR, "test_qc2.pdf")
    try:
        OutputChromCallrate(chrom_calls, fname)
    except:
        pytest.fail("Unexpected error")

def test_main():
    # correct vcf
    args = base_argparse()
    args.vcf = os.path.join(QCDIR, "test_popstr.vcf")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # vcf file with no contig
    args = base_argparse()
    args.vcf = os.path.join(QCDIR, "test_popstr_nocontig.vcf")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0
    
    # Set sample list 
    args.samples = os.path.join(QCDIR, "test_samplelist.txt")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # Trying to test line 186 but seems to not work. TODO update with something else
    # Set sample list with sample that doesn't exist in VCF (should skip it)
    args.samples = os.path.join(QCDIR, "test_samplelist.txt")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # Non existent vcf
    args = base_argparse()
    args.vcf = os.path.join(QCDIR, "test_non_exist.vcf")
    retcode = main(args)
    assert retcode == 1    
