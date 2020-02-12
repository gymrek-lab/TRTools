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

# Set up dummy class
class DummyRecord:
    def __init__(self, chrom, pos, ref, alts=[], info = {}):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALTS = alts
        self.INFO = info


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


def test_MissingFieldWarnings(capsys):
    args = base_argparse()
    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_missinginfo1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    main(args)
    captured = capsys.readouterr()
    assert "Expected info field STUTTERP not found" in captured.err

    fname1 = os.path.join(MRGVCFDIR, "test_file_gangstr_missingformat1.vcf.gz")
    fname2 = os.path.join(MRGVCFDIR, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    main(args)
    captured = capsys.readouterr()
    assert "Expected format field DP not found" in captured.err

def test_ConflictingRefs():
    # Set up dummy records
    dummy_records = [] 
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyRecord('chr1', 100, 'CAG'))

    with pytest.raises(ValueError) as info:
        GetRefAllele(dummy_records, [True, True, True])
    assert "Conflicting refs found at chr1:100" in str(info.value)

    retval = GetRefAllele(dummy_records, [True, True, False])
    assert retval == "CAGCAG"

def test_GetInfoItem(capsys):
    # Set up dummy records
    dummy_records = [] 
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 110}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={}))

    GetInfoItem(dummy_records, [True, True, True, False], 'END')
    captured = capsys.readouterr()
    assert "Incompatible info field value END" in captured.err

    with pytest.raises(ValueError) as info:
        GetInfoItem(dummy_records, [True, True, False, True], 'END')
    assert "Missing info field END" in str(info.value)

    retval = GetInfoItem(dummy_records, [True, True, False, False], 'END')
    assert retval == "END=120"

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
    
    PrintCurrentRecords(dummy_records, [True, True])
    captured = capsys.readouterr()
    assert "Missing CHROM and POS in record" in captured.err

def test_CheckMin():
    assert CheckMin([True, False]) == False
    with pytest.raises(ValueError) as info:
        CheckMin([False, False])
    assert "Unexpected error. Stuck in infinite loop and exiting." in str(info.value)

def test_CheckVCFType():
    snps_path = os.path.join(VCFDIR, "snps.vcf")
    gangstr_path = os.path.join(VCFDIR, "test_gangstr.vcf")
    hipstr_path = os.path.join(VCFDIR, "test_hipstr.vcf")
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    hipstr_vcf = vcf.Reader(filename=hipstr_path)
    snps_vcf = vcf.Reader(filename=snps_path)
    # TODO add tests to infer vcf type
    assert "gangstr" == GetVCFType([gangstr_vcf], "gangstr")

    with pytest.raises(ValueError) as info:
        print(GetVCFType([gangstr_vcf, hipstr_vcf], "auto"))
    assert "VCF files are of mixed types." in str(info.value)
    
    with pytest.raises(ValueError) as info:
        print(GetVCFType([gangstr_vcf, snps_vcf], "auto"))
    assert "Could not identify the type of this vcf" in str(info.value)

def test_GetSampleInfo():
    # TODO add more in depth tests (Create a better dummy class or import from vcf files)
    
    gangstr_path = os.path.join(VCFDIR, "test_gangstr.vcf")
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    record = next(gangstr_vcf)
    
    args = base_argparse()
    
    # TODO final percent!!!
    #for sample in record:
    #    with pytest.raises(ValueError) as info:
    #        GetSampleInfo(record, sample.gt_bases.split(sample.gt_phase_char()), ['UNKNOWNFORMAT'], args)
    #    print(info.traceback)
    #    assert "lolz" in str(info.value)
