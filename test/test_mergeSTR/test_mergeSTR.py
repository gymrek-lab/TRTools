import argparse
import os, sys
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','mergeSTR'))
from mergeSTR import * 

TESTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_files")
COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
DUMPDIR = os.path.join(COMMDIR, "dump")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")
MRGVCFDIR = os.path.join(VCFDIR, "mergeSTR_vcfs")

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
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname1 = os.path.join(MRGVCFDIR, "test_non_existent1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_non_existent2.vcf.gz")
    if os.path.exists(fname1):
        os.remove(fname1)
    if os.path.exists(fname2):
        os.remove(fname2)
    print (os.path.isfile(fname1))
    with pytest.raises(ValueError) as info:
        LoadReaders([fname1, fname2])
    assert "Could not find VCF file" in str(info.value)

# Test right files or directory - GangSTR
def test_GangSTRRightFile():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcftype = "gangstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

# TODO fails bc no contig line in VCFs
# Test right files or directory - advntr
#def test_AdVNTRRightFile():
#    args = base_argparse()
#    fname1 = os.path.join(MRGVCFDIR, "test_file_advntr1.vcf.gz")
#    fname2 = os.path.join(MRGVCFDIR, "test_file_advntr2.vcf.gz")
#    args.vcftype = "advntr"
#    args.vcfs = fname1 + "," + fname2
#    assert main(args)==0
#    args.vcftype = "auto"
#    assert main(args)==0
#    args.update_sample_from_file = True
#    assert main(args)==0
#    args.verbose = True
#    assert main(args)==0

# Test right files or directory - hipstr
def test_hipSTRRightFile():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_hipstr2.vcf.gz")
    args.vcftype = "hipstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

# Test right files or directory - ExpansionHunter
# TODO fails bc no contig line in VCFs
#def test_ExpansionHunterRightFile():
#    args = base_argparse()
#    fname1 = os.path.join(MRGVCFDIR, "test_file_eh1.vcf.gz")
#    fname2 = os.path.join(MRGVCFDIR, "test_file_eh2.vcf.gz")
#    args.vcftype = "eh"
#    args.vcfs = fname1 + "," + fname2
#    assert main(args)==0
#    args.vcftype = "auto"
#    assert main(args)==0
#    args.update_sample_from_file = True
#    assert main(args)==0
#    args.verbose = True
#    assert main(args)==0

def test_GangSTRDuplicate():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr1.vcf.gz")
    args.vcfs = fname1 + "," + fname1
    assert main(args)==1

# Test right files or directory - popstr
def test_PopSTRRightFile():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_popstr1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_popstr2.vcf.gz")
    args.vcftype = "popstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

    
# test unzipped, unindexed VCFs return 1
def test_UnzippedUnindexedFile():
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_unzipped1.vcf")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_unzipped2.vcf")
    with pytest.raises(ValueError) as info:
        LoadReaders([fname1, fname2])
    assert "is bgzipped and indexed" in str(info.value)

    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_unindexed1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_unindexed2.vcf.gz")
    with pytest.raises(ValueError) as info:
        LoadReaders([fname1, fname2])
    assert "Could not find VCF index" in str(info.value)

# test VCFs with different ref genome contigs return 1
def test_WrongContigFile():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_wrongcontig1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_wrongcontig2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "is not in list" in str(info.value)

    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_wrongcontig3.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_wrongcontig4.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "Different contigs found across VCF files." in str(info.value)


def test_MissingFieldWarnings():
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_missinginfo1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.warns(RuntimeWarning) as info:
        main(args)
    assert "Expected info field STUTTERP not found" in str(info[0].message.args[0])

    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_missingformat1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.warns(RuntimeWarning) as info:
        main(args)
    assert "Expected format field DP not found" in str(info[0].message.args[0])
