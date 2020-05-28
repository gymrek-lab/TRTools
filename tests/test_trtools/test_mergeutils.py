import os, sys
import numpy as np
import pytest
import vcf
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..'))

import trtools.utils.mergeutils as mergeutils

COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")
MRGVCFDIR = os.path.join(VCFDIR, "mergeSTR_vcfs")

# Set up dummy class
class DummyRecord:
    def __init__(self, chrom, pos, ref, alts=[], info = {}):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALTS = alts
        self.INFO = info

def test_PrintCurrentRecords(capsys):
    # Set up dummy class
    class DummyRecordNoChrom:
        def __init__(self, chrom, pos, ref, alts=[], info = {}):
            self.POS = pos
            self.REF = ref
            self.ALTS = alts
            self.INFO = info
    # Set up dummy records
    dummy_records = [] 
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyRecordNoChrom('chr1', 100, 'CAGCAG', info={'END': 120}))
    
    mergeutils.PrintCurrentRecords(dummy_records, [True, True])
    captured = capsys.readouterr()
    assert "Missing CHROM and POS in record" in captured.err

def test_CheckMin():
    assert mergeutils.CheckMin([True, False]) == False
    with pytest.raises(ValueError) as info:
        mergeutils.CheckMin([False, False])
    assert "Unexpected error. Stuck in infinite loop and exiting." in str(info.value)

def test_CheckVCFType():
    snps_path = os.path.join(VCFDIR, "snps.vcf")
    gangstr_path = os.path.join(VCFDIR, "test_gangstr.vcf")
    hipstr_path = os.path.join(VCFDIR, "test_hipstr.vcf")
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    hipstr_vcf = vcf.Reader(filename=hipstr_path)
    snps_vcf = vcf.Reader(filename=snps_path)
    # TODO add tests to infer vcf type
    assert "gangstr" == mergeutils.InferAndCheckVCFType([gangstr_vcf], "gangstr")

    with pytest.raises(ValueError) as info:
        print(mergeutils.InferAndCheckVCFType([gangstr_vcf, hipstr_vcf], "auto"))
    assert "VCF files are of mixed types." in str(info.value)
    
    with pytest.raises(TypeError) as info:
        print(mergeutils.InferAndCheckVCFType([gangstr_vcf, snps_vcf], "auto"))
    assert "Could not identify the type of this vcf" in str(info.value)

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
        mergeutils.LoadReaders([fname1, fname2])
    assert "Could not find VCF file" in str(info.value)

# test unzipped, unindexed VCFs return 1
def test_UnzippedUnindexedFile():
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_unzipped1.vcf")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_unzipped2.vcf")
    with pytest.raises(ValueError) as info:
        mergeutils.LoadReaders([fname1, fname2])
    assert "is bgzipped and indexed" in str(info.value)

    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_unindexed1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr_unindexed2.vcf.gz")
    with pytest.raises(ValueError) as info:
        mergeutils.LoadReaders([fname1, fname2])
    assert "Could not find VCF index" in str(info.value)
