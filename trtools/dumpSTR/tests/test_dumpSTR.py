import argparse
import gzip
import os

import pytest

from ..dumpSTR import *
from trtools.testsupport.utils import assert_same_vcf, assert_same_file


# Set up base argparser
@pytest.fixture(name='args')
def args_fixture(tmpdir):
    return args(tmpdir)

def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.vcftype = "auto"
    args.out = str(tmpdir / "test")
    args.zip = False
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
    #args.gangstr_require_support = None
    #args.gangstr_readlen = None
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
def testDumpSTRdir(vcfdir):
    return vcfdir + "/dumpSTR_vcfs"

# Test no such file or directory
def test_WrongFile(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1


# Test a file that already has Filter IDs defined
# that we want to use that are of either the wrong number of type.
# Since cyvcf2 currently won't allow us to overwrite them,
# error out
def test_BadPreexistingFields(args, testDumpSTRdir, capsys):
    fname = os.path.join(testDumpSTRdir, "bad_preexisting_hrun.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode == 1
    captured = capsys.readouterr()
    assert "HRUN" in captured.err

    fname = os.path.join(testDumpSTRdir, "bad_preexisting_het_hwep.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode == 1
    captured = capsys.readouterr()
    assert "HWEP" in captured.err and "HET" in captured.err

    fname = os.path.join(testDumpSTRdir, "bad_preexisting_filter_ac_refac.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode == 1
    captured = capsys.readouterr()
    assert ("FILTER" in captured.err and "AC" in captured.err
            and "REFAC" in captured.err)


# Test a file that already has a HWE Filter ID defined
# if the field is of the correct type and number, as in this case
# we overwrite it and emit a warning instead of failing
# this allows dumpSTR to be run multiple times in succession
# on the same file
def test_WorrisomePreexistingFilter(args, testDumpSTRdir, capsys):
    fname = os.path.join(testDumpSTRdir, "worrisome_preexisting_filter.vcf")
    args.vcf = fname
    args.min_locus_hwep = 0.5
    retcode = main(args)
    assert retcode == 0
    captured = capsys.readouterr()
    assert 'HWE0.5' in captured.err

# Test if basic inputs and threshold filters work for each file
def test_GangSTRFile(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.gangstr_min_call_DP = 10
    args.gangstr_max_call_DP = 20
    args.gangstr_min_call_Q = 0.99
    args.gangstr_filter_span_only = True
    args.gangstr_filter_spanbound_only = True
    args.gangstr_filter_badCI = True
    #args.gangstr_require_support = 2
    #args.gangstr_readlen = 100
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

def test_HipSTRFile(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "trio_chr21_hipstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_min_call_DP = 10
    args.hipstr_max_call_DP = 100
    args.hipstr_min_call_Q = 0.9
    args.hipstr_min_supp_reads = 2
    args.hipstr_max_call_flank_indel = 0.05
    args.hipstr_max_call_stutter = 0.01
    args.vcftype = 'hipstr'
    retcode = main(args)
    assert retcode==0

def test_AdVNTRFile(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_advntr.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.advntr_min_call_DP = 10
    args.advntr_max_call_DP = 20
    args.advntr_min_spanning = 2
    args.advntr_min_flanking = 2
    args.advntr_min_ML = 0
    retcode = main(args)
    assert retcode==0

def test_EHFile(args, testDumpSTRdir):
    # TODO add EH options
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_eh.sorted.vcf.gz")
    args.vcf = fname
    args.use_length = True
    args.num_records = 10
    retcode = main(args)
    assert retcode==0

def test_PopSTRFile(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_popstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.use_length = True
    args.popstr_min_call_DP = 5
    args.popstr_max_call_DP = 100
    args.popstr_require_support = 2
    retcode = main(args)
    assert retcode==0

# confirm that producing zipped output doesn't crash
def test_zippedOutput(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.gangstr_min_call_DP = 10
    args.gangstr_max_call_DP = 20
    args.gangstr_min_call_Q = 0.99
    args.gangstr_filter_span_only = True
    args.gangstr_filter_spanbound_only = True
    args.gangstr_filter_badCI = True
    #args.gangstr_require_support = 2
    #args.gangstr_readlen = 100
    args.zip = True
    retcode = main(args)
    assert retcode==0

# Test invalid options
def test_InvalidOptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_popstr.sorted.vcf.gz")
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
def test_LocusLevel(args, testDumpSTRdir):
    tool_files = [
        "trio_chr21_hipstr.sorted.vcf.gz",
        "trio_chr21_gangstr.sorted.vcf.gz",
        "NA12878_chr21_eh.sorted.vcf.gz",
        "NA12878_chr21_popstr.sorted.vcf.gz",
        "NA12878_chr21_popstr.sorted.vcf.gz",
        "NA12878_chr21_advntr.sorted.vcf.gz"
    ]
    for fname in tool_files:
        args.vcf = os.path.join(testDumpSTRdir, fname)
        args.num_records = 10
        args.min_locus_callrate = 0.8
        args.min_locus_hwep = 10e-4
        args.min_locus_het = 0.1
        args.max_locus_het = 0.3
        args.use_length = True
        args.drop_filtered = False
        args.filter_hrun = True
        if 'hipstr' in fname:
            args.vcftype = 'hipstr'
        else:
            args.vcftype = 'auto'
        assert main(args)==0
        args.drop_filtered = True
        assert main(args)==0

def test_RegionFilters(args, regiondir, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_gangstr.vcf.gz")
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
    args.vcf = os.path.join(testDumpSTRdir, "test_gangstr_nochr.vcf.gz")
    assert main(args)==0

def test_InvalidHipstrOptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "trio_chr21_hipstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_max_call_flank_indel = -1
    args.vcftype = 'hipstr'
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

def test_InvalidGangSTROptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_gangstr.vcf.gz")
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
    '''
    args.gangstr_require_support = -1
    assert main(args)==1
    args.gangstr_require_support = 2
    assert main(args)==1
    args.gangstr_readlen = 1
    assert main(args)==1
    '''

def test_InvalidAdVNTROptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_advntr.vcf.gz")
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
def test_InvalidEHOptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "test_ExpansionHunter.vcf")
    args.vcf = fname
    args.num_records = 10
    # TODO add once EH is implemented
"""

def test_InvalidPopSTROptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_popstr.sorted.vcf.gz")
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

def test_InvalidGenotyperOptions(args, testDumpSTRdir):
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_popstr.sorted.vcf.gz")
    args.vcf = fname
    args.num_records = 10
    args.hipstr_min_call_DP = 10
    assert main(args)==1
    args.hipstr_min_call_DP = None

    args.gangstr_min_call_DP = 10
    assert main(args)==1
    args.gangstr_min_call_DP = None

    fname = os.path.join(testDumpSTRdir, "trio_chr21_hipstr.sorted..vcf.gz")
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

def test_InvalidOutput(capsys, args, testDumpSTRdir, tmpdir):
    fname = os.path.join(testDumpSTRdir, "NA12878_chr21_popstr.sorted.vcf.gz")
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
    assert 'is a directory' in str(capsys.readouterr())

def test_TwoDumpSTRRounds(args, testDumpSTRdir, tmpdir):
    args.num_records = 10
    fname = os.path.join(testDumpSTRdir, "test_gangstr.vcf.gz")
    args.vcf = fname
    args.min_locus_callrate = 0
    args.zip = True
    main(args) # produces DUMPDIR/test.vcf
    args.vcf = str(tmpdir / "test.vcf.gz")
    args.out = str(tmpdir / "test2")
    assert main(args)==0

def test_BrokenVCF(args, testDumpSTRdir):
    args.num_records = 10
    fname = os.path.join(testDumpSTRdir, "test_broken.vcf.gz")
    args.vcf = fname
    args.die_on_warning = True
    args.verbose = True
    assert main(args)==1

def test_beagle_allowed_locus_filters(args, vcfdir, regiondir):#(tmpdir):
    args.min_locus_hwep = 0.1
    args.min_locus_het = 0.1
    args.max_locus_het = 0.9
    args.filter_regions = os.path.join(regiondir, "test_regions1.bed.gz")
    for caller in 'advntr', 'eh', 'gangstr', 'hipstr':
        args.vcf = os.path.join(vcfdir, f'beagle/{caller}_imputed.vcf.gz')
        assert main(args) == 0

def test_beagle_disallowed_locus_filters(tmpdir, vcfdir):
    args_obj = args(tmpdir)
    args_obj.min_locus_callrate = 0.1
    for caller in 'advntr', 'eh', 'gangstr', 'hipstr':
        args_obj.vcf = os.path.join(vcfdir, f'beagle/{caller}_imputed.vcf.gz')
        assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.filter_hrun = True
    assert main(args_obj) == 1

def test_beagle_disallowed_call_filters(tmpdir, vcfdir):
    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_min_call_DP = 5
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_max_call_DP = 1000
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_min_call_Q = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_max_call_flank_indel = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_max_call_stutter = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/hipstr_imputed.vcf.gz')
    args_obj.hipstr_min_supp_reads = 2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_expansion_prob_het = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_expansion_prob_hom = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_expansion_prob_total = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_filter_span_only = True
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_filter_spanbound_only = True
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_filter_badCI = True
    assert main(args_obj) == 1

    #args_obj = args(tmpdir)
    #args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    #args_obj.gangstr_require_support = None
    #args_obj.gangstr_readlen = None
    #assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_min_call_DP = 2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_max_call_DP = 1000
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/gangstr_imputed.vcf.gz')
    args_obj.gangstr_min_call_Q = 0.2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/advntr_imputed.vcf.gz')
    args_obj.advntr_min_call_DP = 2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/advntr_imputed.vcf.gz')
    args_obj.advntr_max_call_DP = 1000
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/advntr_imputed.vcf.gz')
    args_obj.advntr_min_spanning = 2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/advntr_imputed.vcf.gz')
    args_obj.advntr_min_flanking = 2
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/advntr_imputed.vcf.gz')
    args_obj.advntr_min_ML = 0.1
    assert main(args_obj) == 1

    '''
    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/eh_imputed.vcf.gz')
    args_obj.eh_min_ADFL = None
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/eh_imputed.vcf.gz')
    args_obj.eh_min_ADIR = None
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/eh_imputed.vcf.gz')
    args_obj.eh_min_ADSP = None
    assert main(args_obj) == 1
    '''

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/eh_imputed.vcf.gz')
    args_obj.eh_min_call_LC = 1
    assert main(args_obj) == 1

    args_obj = args(tmpdir)
    args_obj.vcf = os.path.join(vcfdir, 'beagle/eh_imputed.vcf.gz')
    args_obj.eh_max_call_LC = 50
    assert main(args_obj) == 1


"""
These tests run dumpSTR and compare its output
to output that has been generated by a pervious version of 
dumpSTR and saved in the repo. The results are expected
to be identical.

These tests are too strict and will often break because
dumpSTR output has been intentionally changed
However, the presence of these tests is important because
it should prevent any unexpected changes in output.
If you've reviewed the change in output and find it acceptable, 
use trtools/testsupport/sample_vcfs/dumpSTR_vcfs/create_test_files.sh
to regenerate the tests files with the new output.
"""


def test_output_locus_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/trio_chr21_hipstr.sorted.vcf.gz'
    args.min_locus_callrate = 0.5
    args.min_locus_hwep = 0.5
    args.min_locus_het = 0.05
    args.max_locus_het = 0.45
    args.filter_regions_names = 'foo_region'
    args.filter_regions = testDumpSTRdir + '/sample_region.bed.gz'
    args.vcftype = 'hipstr'

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    # there are also rounding errors with HipSTR field GLDIFF
    # that aren't worth worrying about
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/locus_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/locus_filters' + ext,
                          ext)


# make sure locus level filters produce the same output when
# --drop-filtered is set
def test_output_drop_filtered(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/trio_chr21_hipstr.sorted.vcf.gz'
    args.min_locus_callrate = 0.5
    args.min_locus_hwep = 0.5
    args.min_locus_het = 0.05
    args.max_locus_het = 0.45
    args.filter_regions_names = 'foo_region'
    args.filter_regions = testDumpSTRdir + '/sample_region.bed.gz'
    args.vcftype = 'hipstr'
    args.drop_filtered = True

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    # there are also rounding errors with HipSTR field GLDIFF
    # that aren't worth worrying about
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/drop_filtered.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/locus_filters' + ext,
                          ext)


# test advntr call level filters
def test_output_advntr_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/NA12878_chr21_advntr.sorted.vcf.gz'
    args.advntr_min_call_DP = 50
    args.advntr_max_call_DP = 2000
    args.advntr_min_spanning = 1
    args.advntr_min_flanking = 20
    args.advntr_min_ML = 0.95

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/advntr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/advntr_filters' + ext,
                          ext)


# test hipstr call and locus level filters
def test_output_hipstr_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/trio_chr21_hipstr.sorted.vcf.gz'
    args.filter_hrun = True
    args.use_length = True
    args.max_locus_het = 0.45
    args.min_locus_het = 0.05
    args.min_locus_hwep = 0.5
    args.hipstr_max_call_flank_indel = 0.05
    args.hipstr_max_call_stutter = 0.3
    args.hipstr_min_supp_reads = 10
    args.hipstr_min_call_DP = 30
    args.hipstr_max_call_DP = 200
    args.hipstr_min_call_Q = 0.9
    args.vcftype = 'hipstr'

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    # there are also rounding errors with HipSTR field GLDIFF
    # that aren't worth worrying about
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/hipstr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/hipstr_filters' + ext,
                          ext)


# test gangstr call level filters that don't begin
# with 'expansion' - those are tested on another file
def test_output_gangstr_most_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/trio_chr21_gangstr.sorted.vcf.gz'
    args.gangstr_min_call_DP = 10
    args.gangstr_max_call_DP = 100
    args.gangstr_min_call_Q = 0.9
    args.gangstr_filter_span_only = True
    args.gangstr_filter_spanbound_only = True
    args.gangstr_filter_badCI = True
    # args.gangstr_require_support = 10
    # args.gangstr_readlen = 150

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/gangstr_filters_most.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/gangstr_filters_most' + ext,
                          ext)


# test gangstr call level filters that begin with
# 'expansion' - the other gangstr call level filters
# are tested on another file
def test_output_gangstr_expansion_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/test_gangstr.vcf.gz'
    args.gangstr_expansion_prob_het = 0.001
    args.gangstr_expansion_prob_hom = 0.0005
    args.gangstr_expansion_prob_total =  0.001

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/gangstr_filters_expansion.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/gangstr_filters_expansion' + ext,
                          ext)


# test popstr call level filters
def test_output_popstr_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/NA12878_chr21_popstr.sorted.vcf.gz'
    args.popstr_min_call_DP = 30
    args.popstr_max_call_DP = 200
    args.popstr_require_support = 15
    args.use_length = True

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/popstr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        assert_same_file(args.out + ext,
                          testDumpSTRdir + '/popstr_filters' + ext,
                          ext)

