import os, sys
import numpy as np
import pytest
import vcf

import trtools.utils.mergeutils as mergeutils


@pytest.fixture
def mrgvcfdir(vcfdir):
	return os.path.join(vcfdir, "mergeSTR_vcfs")

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
    
    mergeutils.DebugPrintRecordLocations(dummy_records, [True, True])
    captured = capsys.readouterr()
    assert "Missing CHROM and POS in record" in captured.err

def test_CheckMin():
    assert mergeutils.CheckMin([True, False]) == False
    with pytest.raises(ValueError) as info:
        mergeutils.CheckMin([False, False])
    assert "Unexpected error. Stuck in infinite loop and exiting." in str(info.value)

def test_CheckVCFType(vcfdir):
    snps_path = os.path.join(vcfdir, "snps.vcf")
    gangstr_path = os.path.join(vcfdir, "test_gangstr.vcf")
    hipstr_path = os.path.join(vcfdir, "test_hipstr.vcf")
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    hipstr_vcf = vcf.Reader(filename=hipstr_path)
    snps_vcf = vcf.Reader(filename=snps_path)
    # TODO add tests to infer vcf type
    assert "gangstr" == mergeutils.GetAndCheckVCFType([gangstr_vcf], "gangstr")

    with pytest.raises(ValueError) as info:
        print(mergeutils.GetAndCheckVCFType([gangstr_vcf, hipstr_vcf], "auto"))
    assert "VCF files are of mixed types." in str(info.value)
    
    with pytest.raises(TypeError) as info:
        print(mergeutils.GetAndCheckVCFType([gangstr_vcf, snps_vcf], "auto"))
    assert "Could not identify the type of this vcf" in str(info.value)

# Test no such file or directory
def test_WrongFile(mrgvcfdir):
    # Try a dummy file name. Make sure it doesn't exist before we try
    fname1 = os.path.join(mrgvcfdir, "test_non_existent1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_non_existent2.vcf.gz")
    if os.path.exists(fname1):
        os.remove(fname1)
    if os.path.exists(fname2):
        os.remove(fname2)
    print (os.path.isfile(fname1))
    with pytest.raises(ValueError) as info:
        mergeutils.LoadReaders([fname1, fname2])
    assert "Could not find VCF file" in str(info.value)

# test unzipped, unindexed VCFs return 1
def test_UnzippedUnindexedFile(mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr_unzipped1.vcf")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr_unzipped2.vcf")
    with pytest.raises(ValueError) as info:
        mergeutils.LoadReaders([fname1, fname2])
    assert "is bgzipped and indexed" in str(info.value)

    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr_unindexed1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr_unindexed2.vcf.gz")
    with pytest.raises(ValueError) as info:
        mergeutils.LoadReaders([fname1, fname2])
    assert "Could not find VCF index" in str(info.value)
