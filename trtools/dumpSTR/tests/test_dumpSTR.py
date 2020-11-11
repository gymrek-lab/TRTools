import argparse
import os

import pytest
import vcf

from ..dumpSTR import *


# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.vcftype = "auto"
    args.out = str(tmpdir / "test")
    args.min_locus_callrate = None
    args.min_locus_hwep = None
    args.min_locus_het = None
    args.max_locus_het = None
    args.use_length = False
    args.filter_regions = None
    args.filter_regions_names = None
    args.filter_hrun = False
    args.drop_filtered = False
    args.hipstr_min_call_DP = None
    args.hipstr_max_call_DP = None
    args.hipstr_min_call_Q = None
    args.hipstr_max_call_flank_indel = None
    args.hipstr_max_call_stutter = None
    args.hipstr_min_supp_reads = None
    args.gangstr_expansion_prob_het = None
    args.gangstr_expansion_prob_hom = None
    args.gangstr_expansion_prob_total = None
    args.gangstr_filter_span_only = False
    args.gangstr_filter_spanbound_only = False
    args.gangstr_filter_badCI = None
    args.gangstr_require_support = None
    args.gangstr_readlen = None
    args.gangstr_min_call_DP = None
    args.gangstr_max_call_DP = None
    args.gangstr_min_call_Q = None
    args.advntr_min_call_DP = None
    args.advntr_max_call_DP = None
    args.advntr_min_spanning = None
    args.advntr_min_flanking = None
    args.advntr_min_ML = None
    args.eh_min_ADFL = None
    args.eh_min_ADIR = None
    args.eh_min_ADSP = None
    args.eh_min_call_LC = None
    args.eh_max_call_LC = None
    args.popstr_min_call_DP = None
    args.popstr_max_call_DP = None
    args.popstr_require_support = None
    args.num_records = None
    args.die_on_warning = False
    args.verbose = False
    return args

@pytest.fixture
def testfiledir(vcfdir):
    return vcfdir + "/dumpSTR_vcfs"

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

# Test if basic inputs and threshold filters work for each file
def test_GangSTRFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.gangstr_min_call_DP = 10
    args.gangstr_max_call_DP = 20
    args.gangstr_min_call_Q = 0.99
    args.gangstr_filter_span_only = True
    args.gangstr_filter_spanbound_only = True
    args.gangstr_filter_badCI = True
    args.gangstr_require_support = 2
    args.gangstr_readlen = 100
    retcode = main(args)
    assert retcode==0

    # Test expansion options
    args.gangstr_expansion_prob_het = 0.8
    retcode = main(args)
    assert retcode==0

    args.gangstr_expansion_prob_het = None
    args.gangstr_expansion_prob_hom = 0.8
    retcode = main(args)
    assert retcode==0

    args.gangstr_expansion_prob_het = None
    args.gangstr_expansion_prob_hom = None
    args.gangstr_expansion_prob_total = 0.8
    retcode = main(args)
    assert retcode==0

def test_HipSTRFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_min_call_DP = 10
    args.hipstr_max_call_DP = 100
    args.hipstr_min_call_Q = 0.9
    args.hipstr_min_supp_reads = 2
    args.hipstr_max_call_flank_indel = 0.05
    args.hipstr_max_call_stutter = 0.01
    retcode = main(args)
    assert retcode==0

def test_AdVNTRFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_advntr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.advntr_min_call_DP = 10
    args.advntr_max_call_DP = 20
    args.advntr_min_spanning = 2
    args.advntr_min_flanking = 2
    args.advntr_min_ML = 0
    retcode = main(args)
    assert retcode==0

# TODO: uncomment. EH not implemented yet in TR Harmonizer
"""
def test_EHFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_ExpansionHunter.vcf")
    args.vcf = fname
    args.num_records = 10
    retcode = main(args)
    assert retcode==0
"""

def test_PopSTRFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_popstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.popstr_min_call_DP = 5
    args.popstr_max_call_DP = 100
    args.popstr_require_support = 2
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode==0

# Test invalid options
def test_InvalidOptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_popstr.vcf")
    args.vcf = fname
    # HWE
    args.min_locus_hwep = -1
    retcode = main(args)
    assert retcode==1
    args.min_locus_hwep = 2
    retcode = main(args)
    assert retcode==1
    # Het
    args.min_locus_hwep = None
    args.min_locus_het = -1
    retcode = main(args)
    assert retcode==1
    args.min_locus_het = 2
    retcode = main(args)
    assert retcode==1
    args.min_locus_het = None
    args.max_locus_het = -1
    retcode = main(args)
    assert retcode==1
    args.max_locus_het = 2
    retcode = main(args)
    assert retcode==1
    args.min_locus_het = 0.5
    args.max_locus_het = 0.2
    retcode = main(args)
    assert retcode==1

# Test locus-level filters
def test_LocusLevel(args, vcfdir):
    for tool in ["hipstr","gangstr","popstr","advntr"]:
        fname = os.path.join(vcfdir, "test_%s.vcf"%tool)
        args.vcf = fname
        args.num_records = 10
        args.min_locus_callrate = 0.8
        args.min_locus_hwep = 10e-4
        args.min_locus_het = 0.1
        args.max_locus_het = 0.3
        args.use_length = True
        args.drop_filtered = False
        args.filter_hrun = True
        assert main(args)==0
        args.drop_filtered = True
        assert main(args)==0

def test_RegionFilters(args, regiondir, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.num_records = 10
    # Correct filters
    args.filter_regions = os.path.join(regiondir, "test_regions1.bed.gz")
    retcode = main(args)
    assert retcode==0
    args.filter_regions_names = "test"
    retcode = main(args)
    assert retcode==0
    # Correct filters, multiple regions
    args.filter_regions = os.path.join(regiondir, "test_regions1.bed.gz") + "," + os.path.join(regiondir, "test_regions2.bed.gz")
    args.filter_regions_names = "test1,test2"
    retcode = main(args)
    assert retcode==0
    # Mismatch between region names and regions
    args.filter_regions_names = "test1"
    retcode = main(args)
    assert retcode==1
    # Nonexistent regions file
    args.filter_regions = os.path.join(regiondir, "test_nonexistent.bed")
    retcode = main(args)
    assert retcode==1
    # File missing tabix
    args.filter_regions = os.path.join(regiondir, "test_regions3.bed.gz")
    assert main(args)==1
    # File with no chr
    args.filter_regions = os.path.join(regiondir, "test_regions4.bed.gz")
    assert main(args)==0
    args.vcf = os.path.join(vcfdir, "test_gangstr_nochr.vcf")
    assert main(args)==0

def test_InvalidHipstrOptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_max_call_flank_indel = -1
    retcode = main(args)
    assert retcode==1
    args.hipstr_max_call_flank_indel = None
    args.hipstr_max_call_flank_indel = 2
    retcode = main(args)
    assert retcode==1
    args.hipstr_max_call_flank_indel = None
    args.hipstr_max_call_stutter = -1
    retcode = main(args)
    assert retcode==1
    args.hipstr_max_call_stutter = 2
    retcode = main(args)
    assert retcode==1
    args.hipstr_max_call_stutter = None
    args.hipstr_min_supp_reads = -1
    retcode = main(args)
    assert retcode==1
    args.hipstr_min_supp_reads = None
    args.hipstr_min_call_DP = -1
    assert main(args)==1
    args.hipstr_min_call_DP = None
    args.hipstr_max_call_DP = -1
    assert main(args)==1
    args.hipstr_min_call_DP = 5
    args.hipstr_max_call_DP = 2
    assert main(args)==1
    args.hipstr_min_call_DP = None
    args.hipstr_max_call_DP = None
    args.hipstr_min_call_Q = -1
    assert main(args)==1
    args.hipstr_min_call_Q = 2
    assert main(args)==1

def test_InvalidGangSTROptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.gangstr_min_call_DP = -1
    assert main(args)==1
    args.gangstr_min_call_DP = None
    args.gangstr_max_call_DP = -1
    assert main(args)==1
    args.gangstr_min_call_DP = 5
    args.gangstr_max_call_DP = 2
    assert main(args)==1
    args.gangstr_min_call_DP = None
    args.gangstr_max_call_DP = None
    args.gangstr_min_call_Q = -1
    assert main(args)==1
    args.gangstr_min_call_Q = 2
    assert main(args)==1
    args.gangstr_min_call_Q = None
    args.gangstr_expansion_prob_het = -1
    assert main(args)==1
    args.gangstr_expansion_prob_het = 2
    assert main(args)==1
    args.gangstr_expansion_prob_het = None
    args.gangstr_expansion_prob_hom = -1
    assert main(args)==1
    args.gangstr_expansion_prob_hom = 2
    assert main(args)==1
    args.gangstr_expansion_prob_hom = None
    args.gangstr_expansion_prob_total = -1
    assert main(args)==1
    args.gangstr_expansion_prob_total = 2
    assert main(args)==1
    args.gangstr_expansion_prob_total = None
    args.gangstr_require_support = -1
    assert main(args)==1
    args.gangstr_require_support = 2
    assert main(args)==1
    args.gangstr_readlen = 1
    assert main(args)==1

def test_InvalidAdVNTROptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_advntr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.advntr_min_call_DP = -1
    assert main(args)==1
    args.advntr_min_call_DP = None
    args.advntr_max_call_DP = -1
    assert main(args)==1
    args.advntr_min_call_DP = 5
    args.advntr_max_call_DP = 2
    assert main(args)==1
    args.advntr_min_call_DP = None
    args.advntr_max_call_DP = None
    args.advntr_min_ML = -1
    assert main(args)==1
    args.advntr_min_ML = None
    args.advntr_min_flanking = -1
    assert main(args)==1
    args.advntr_min_spanning = -1
    assert main(args)==1

"""
def test_InvalidEHOptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_ExpansionHunter.vcf")
    args.vcf = fname
    args.num_records = 10
    # TODO add once EH is implemented
"""

def test_InvalidPopSTROptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_popstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.popstr_min_call_DP = -1
    assert main(args)==1
    args.popstr_min_call_DP = None
    args.popstr_max_call_DP = -1
    assert main(args)==1
    args.popstr_min_call_DP = 5
    args.popstr_max_call_DP = 2
    assert main(args)==1
    args.popstr_min_call_DP = None
    args.popstr_max_call_DP = None
    args.popstr_require_support = -1
    assert main(args)==1

def test_InvalidGenotyperOptions(args, vcfdir):
    fname = os.path.join(vcfdir, "test_popstr.vcf")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_min_call_DP = 10
    assert main(args)==1
    args.hipstr_min_call_DP = None

    args.gangstr_min_call_DP = 10
    assert main(args)==1
    args.gangstr_min_call_DP = None

    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcf = fname
    args.popstr_min_call_DP = 10
    assert main(args)==1
    args.popstr_min_call_DP = None
    args.advntr_min_call_DP = 10
    assert main(args)==1
    args.advntr_min_call_DP = None
    args.eh_min_call_LC = 5
    assert main(args)==1
    args.eh_min_call_LC = None

def test_InvalidOutput(capsys, args, vcfdir, tmpdir):
    fname = os.path.join(vcfdir, "test_popstr.vcf")
    args.vcf = fname

    # Fail when trying to output inside a nonexistant directory
    args.out = str(tmpdir / "notadirectory" / "somefilename")
    assert main(args) == 1

    # To simulate a permissions issue: fail when trying to write a file in a location
    # that is already a directory
    capsys.readouterr()
    (tmpdir / "foo.vcf").mkdir()
    args.out = str(tmpdir / "foo")
    assert main(args) == 1
    # Make sure we produce a meaningful error message for this issue
    assert 'Is a directory:' in str(capsys.readouterr())

def test_TwoDumpSTRRounds(args, vcfdir, tmpdir):
    args.num_records = 10
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.min_locus_callrate = 0
    main(args) # produces DUMPDIR/test.vcf
    args.vcf = str(tmpdir / "test.vcf")
    args.out = str(tmpdir / "test2")
    assert main(args)==0

def test_BrokenVCF(args, vcfdir):
    args.num_records = 10
    fname = os.path.join(vcfdir, "test_broken.vcf")
    args.vcf = fname
    args.die_on_warning = True
    args.verbose = True
    assert main(args)==1




"""
These tests run dumpSTR and compare its output
to output that has been generated by a pervious version of 
dumpSTR and saved in the repo. The results are expected
to be identical.
"""

# fname1 should be the output file
# fname2 should be the control file
def _assert_same_file(fname1, fname2):
    with open(fname1) as file1, open(fname2) as file2:
        iter1 = iter(file1)
        iter2 = iter(file2)
        linenum = 0
        while True:
            file1ended = False
            file2ended = False
            try:
                line1 = next(iter1)
            except StopIteration:
                file1ended = True
            try:
                line2 = next(iter2)
            except StopIteration:
                file2ended = True
            if file1ended != file2ended:
                if file1ended:
                    raise ValueError(
                        'Output file ' + fname1 + ' differs from control file '
                        + fname2 + '. Output file ended after ' + str(linenum) +
                        ' lines, but control file  did not'
                    )
                else:
                    raise ValueError(
                        'Output file ' + fname1 + ' differs from control file '
                        + fname2 + '. Control file ended after ' + str(linenum) +
                        ' lines, but output file did not'
                    )
            if file1ended and file2ended:
                return

            if line1 != line2:
                print('Output file ' + fname1 + ' differs from control file'
                      + fname2 + ' at line ' + str(linenum) + '.\nLine in output'
                      ' file: ' + line1 + '\nLine in control file: ' + line2)


def test_output_locus_filters(args, testfiledir):
    args.vcf = testfiledir + '/trio_chr21_hipstr.sorted.vcf.gz'
    args.min_locus_callrate = 0.5
    args.min_locus_hwp = 0.5
    args.min_locus_het = 0.05
    args.max_locus_het = 0.45
    args.filter_regions_names = 'foo_region'
    args.filter_regions = testfiledir + '/sample_region.bed.gz'
    args.vcftype = 'hipstr'

    assert main(args) == 0
    for ext in '.samplog.tab', '.loclog.tab', '.vcf':
        _assert_same_file(args.out + ext, testfiledir + 'base_filters' + ext)

