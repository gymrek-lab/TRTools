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
def test_Filters(args, vcfdir):
    fname = os.path.join(vcfdir, "artificial_gangstr.vcf")
    args.vcf = fname
    args.vcftype = "gangstr"
    artificial_vcf = vcf.Reader(filename=args.vcf)

    ## Test1: call passes with no filter
    vcfcall = vcf.model._Call # Blank call
    call_filters = []
    reasons1 = FilterCall(vcfcall, call_filters)
    assert reasons1==[]

    # Line 1 of artificial vcf
    record1 = next(artificial_vcf)
    call1 = record1.samples[0]

    ## Check call1 attributes:
    assert record1.CHROM=='chr1'
    assert record1.POS==3004986
    assert record1.REF=='tctgtctgtctg'
    assert record1.INFO['RU']=='tctg'
    assert call1['DP']==31
    assert call1['Q']==0.999912
    assert call1['QEXP']==[0.0, 5.1188e-05, 0.999949]
    assert call1['RC']=='17,12,0,2'

    ## Test2: call filter: LowCallDepth
    vcfcall = call1
    args = base_argparse()
    args.min_call_DP = 50
    call_filters = BuildCallFilters(args)
    reasons2 = FilterCall(vcfcall, call_filters)
    assert reasons2==['LowCallDepth']

    ## Test3: call filter: HighCallDepth
    vcfcall = call1
    args = base_argparse()
    args.max_call_DP = 10
    call_filters = BuildCallFilters(args)
    reasons3 = FilterCall(vcfcall, call_filters)
    assert reasons3==['HighCallDepth']

    ## Test4: call filter: LowCallQ
    vcfcall = call1
    args = base_argparse()
    args.min_call_Q = 1
    call_filters = BuildCallFilters(args)
    reasons4 = FilterCall(vcfcall, call_filters)
    assert reasons4==['LowCallQ']

    ## Test4: call filter: LowCallQ
    vcfcall = call1
    args = base_argparse()
    args.min_call_Q = 1
    call_filters = BuildCallFilters(args)
    reasons4 = FilterCall(vcfcall, call_filters)
    assert reasons4==['LowCallQ']

    ## Test5: call filter: ProbHom
    vcfcall = call1
    args = base_argparse()
    args.expansion_prob_hom = 1
    call_filters = BuildCallFilters(args)
    reasons5 = FilterCall(vcfcall, call_filters)
    assert reasons5==['ProbHom']

    ## Test6: call filter: ProbHet
    vcfcall = call1
    args = base_argparse()
    args.expansion_prob_het = 0.8
    call_filters = BuildCallFilters(args)
    reasons6 = FilterCall(vcfcall, call_filters)
    assert reasons6==['ProbHet']

    # Line 2 of artificial vcf
    record2 = next(artificial_vcf)
    call2 = record2.samples[0]

    ## Check call2 attributes:
    assert record2.CHROM=='chr1'
    assert record2.POS==3005549
    assert record2.REF=='aaaacaaaacaaaacaaaac'
    assert record2.INFO['RU']=='aaaac'
    assert call2['DP']==11
    assert call2['Q']==1
    assert call2['QEXP']==[0.8, 0.2, 0]
    assert call2['RC']=='0,11,0,0'

    ## Test7: call filter: ProbTotal
    vcfcall = call2
    args = base_argparse()
    args.expansion_prob_total = 1
    call_filters = BuildCallFilters(args)
    reasons7 = FilterCall(vcfcall, call_filters)
    assert reasons7==['ProbTotal']

    ## Test8: call filter: filter span only
    vcfcall = call1
    args = base_argparse()
    args.filter_span_only = True
    call_filters = BuildCallFilters(args)
    reasons71 = FilterCall(vcfcall, call_filters)
    assert reasons71==[]

    vcfcall = call2
    args = base_argparse()
    args.filter_span_only = True
    call_filters = BuildCallFilters(args)
    reasons72 = FilterCall(vcfcall, call_filters)
    assert reasons72==['SpanOnly']

    # Line 3 of artificial vcf
    record3 = next(artificial_vcf)
    call3 = record3.samples[0]

    ## Check call2 attributes:
    assert record3.CHROM=='chr1'
    assert record3.POS==3009351
    assert record3.REF=='tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg'
    assert record3.INFO['RU']=='tg'
    assert call3['DP']==20
    assert call3['RC']=='0,10,0,10'

    ## Test8: call filter: filter span-bound only
    vcfcall = call1
    args = base_argparse()
    args.filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons81 = FilterCall(vcfcall, call_filters)
    assert reasons81==[]

    vcfcall = call2
    args = base_argparse()
    args.filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons82 = FilterCall(vcfcall, call_filters)
    assert reasons82==['SpanBoundOnly']

    vcfcall = call3
    args = base_argparse()
    args.filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons83 = FilterCall(vcfcall, call_filters)
    assert reasons83==['SpanBoundOnly']
"""
