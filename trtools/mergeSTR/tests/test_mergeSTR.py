import argparse
import os

import cyvcf2
import numpy as np
import pytest

from ..mergeSTR import *
from trtools.testsupport.utils import assert_same_vcf


# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = str(tmpdir / "test")
    args.update_sample_from_file = False 
    args.quiet = False
    args.verbose = False
    args.vcftype = "auto"
    return args


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


class DummyHarmonizedRecord:
    def __init__(self, chrom, pos, ref, alts=None, info=None):
        self.chrom = chrom
        self.pos = pos
        self.ref_allele = ref
        self.alt_alleles = alts if alts is not None else []
        self.info = info if info is not None else {}
        self.vcfrecord = DummyRecord(chrom, pos, ref, self.alt_alleles, self.info)

# Test right files or directory - GangSTR
def test_GangSTRRightFile(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2.vcf.gz")
    args.vcftype = "gangstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

# Test right files or directory - advntr
def test_AdVNTRRightFile(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_advntr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_advntr2.vcf.gz")
    args.vcftype = "advntr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

# Test right files or directory - hipstr
def test_hipSTRRightFile(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_hipstr2.vcf.gz")
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
def test_ExpansionHunterRightFile(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_eh1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh2.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

def test_GangSTRDuplicate(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    args.vcfs = fname1 + "," + fname1
    assert main(args)==1

# Test right files or directory - popstr
def test_PopSTRRightFile(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_popstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_popstr2.vcf.gz")
    args.vcftype = "popstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args)==0
    args.vcftype = "auto"
    assert main(args)==0
    args.update_sample_from_file = True
    assert main(args)==0
    args.verbose = True
    assert main(args)==0

def test_multiple_vcf_types(args, mrgvcfdir, capsys):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_popstr2.vcf.gz")
    args.vcftype = "auto"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 1
    assert 'mixed types' in capsys.readouterr().err

def test_duplicate_ids(args, mrgvcfdir, capsys):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr_dupID.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 1
    assert 'same sample' in capsys.readouterr().err

# test VCFs with different ref genome contigs return 1
def test_RecordChromsNotInContigs(args, mrgvcfdir, capsys):
    #both files have records with chroms not listed in contigs
    fname1 = os.path.join(mrgvcfdir, "test_file_contigmissing1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_contigmissing2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 1
    assert 'not found in the contig list' in capsys.readouterr().err

    #first file has records with chroms not listed in contigs
    #first file missing contig is the first record
    fname1 = os.path.join(mrgvcfdir, "test_file_contigmissing1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2_1contig.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 1
    assert 'not found in the contig list' in capsys.readouterr().err

    #second file has records with chroms not listed in contigs
    #second file missing contig is not the first record
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1_1contig.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_contigmissing2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 1
    assert 'not found in the contig list' in capsys.readouterr().err

def test_DifferentContigs(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_contigdifferent1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_contigdifferent2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "Different contigs found across VCF files." in str(info.value)

def test_DifferentContigLengths(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_contigdifflength.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "Different contigs found across VCF files." in str(info.value)

def test_SameContigsDifferentOrder(args, vcfdir, mrgvcfdir):
    fname1 = os.path.join(vcfdir, "one_sample_multiple_chroms.vcf.gz")
    fname2 = os.path.join(
        mrgvcfdir,
        "one_sample_multiple_chroms_diff_contig_order.vcf.gz"
    )
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0

def test_MissingFieldWarnings(capsys, args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr_missinginfo1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    main(args)
    captured = capsys.readouterr()
    assert "Expected info field STUTTERP not found" in captured.err

    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr_missingformat1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    main(args)
    captured = capsys.readouterr()
    assert "Expected format field DP not found" in captured.err

def test_alt_same_len_as_ref_different_flanking(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_hipstr2_alt_v_ref.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    main(args)

    vcf = cyvcf2.VCF(args.out + '.vcf')
    var = next(vcf)
    for alt in var.ALT:
        assert alt != var.REF

def test_ConflictingRefs():
    # Set up dummy records
    dummy_records = []
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAG'))

    retval = GetRefAllele(dummy_records, [True, True, True], None)
    assert retval is None

    retval = GetRefAllele(dummy_records, [True, True, False], None)
    assert retval == "CAGCAG"

def test_GetInfoItem(capsys):
    # Set up dummy records
    dummy_records = []
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG', info={'END': 110}))
    dummy_records.append(DummyHarmonizedRecord('chr1', 100, 'CAGCAG', info={}))

    GetInfoItem(dummy_records, [True, True, True, False], 'END')
    captured = capsys.readouterr()
    assert ("Incompatible values" in captured.err and
            "info field END" in captured.err)

    with pytest.raises(ValueError) as info:
        GetInfoItem(dummy_records, [True, True, False, True], 'END')
    assert "Missing info field END" in str(info.value)

    retval = GetInfoItem(dummy_records, [True, True, False, False], 'END')
    assert retval == "END=120"

# TODO write WriteSampleData tests

"""
These tests run mergeSTR and compare its output
to output that has been generated by a pervious version of
mergeSTR and saved in the repo. The results are expected
to be identical.

These tests are too strict and will often break because
mergeSTR output has been intentionally changed
However, the presence of these tests is important because
it should prevent any unexpected changes in output.
If you've reviewed the change in output and find it acceptable,
use trtools/testsupport/sample_vcfs/mergeSTR_vcfs/create_test_files.sh
to regenerate the test files with the new version of mergeSTR.
"""


def test_advntr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_advntr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_advntr2.vcf.gz")
    args.vcftype = "advntr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/advntr_merged.vcf")



def test_eh_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_eh1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh2.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/eh_merged.vcf")

def test_eh_no_alt(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_eh1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh_no_alt.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/eh_no_alt_merged.vcf")

def test_eh_mixed_ploidy_no_alt(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_eh_X1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh_X_no_alt.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/eh_X_no_alt_merged.vcf")

    # reverse order
    fname1 = os.path.join(mrgvcfdir, "test_file_eh_X_no_alt.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh_X1.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/eh_X_no_alt_merged_swap.vcf")



def test_gangstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2.vcf.gz")
    args.vcftype = "gangstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/gangstr_merged.vcf")


def test_hipstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_hipstr2.vcf.gz")
    args.vcftype = "hipstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/hipstr_merged.vcf")

def test_hipstr_output_flanking_pb_harmonization(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "hipstr-harmonized-merge-contains-flanking.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "hipstr-harmonized-merge-no-flanking.vcf.gz")
    args.vcftype = "hipstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/hipstr_flanking_harmonization_test_output.vcf")

def test_popstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_popstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_popstr2.vcf.gz")
    args.vcftype = "popstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/popstr_merged.vcf")


# TODO questions and issues to confirm, test or  address:
# confirm conflicting samples cause this to fail unless --update-sample-from-file is given
# required info fields cause failure when a record is missing them but not when the entire VCF is missing them
# we silently return if there are no info or header fields
# we don't fail if info fields are different
# how should we merge info fields with reqd = false
# what if there are multiple records at the same location in the same VCF
# if a genotype is no called but there is other format info,
# we aren't emitting it. Is that intended?
# TODO write test where sample and allele orderings change
