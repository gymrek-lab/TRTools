import argparse
import os,sys
import pytest

from ..dumpSTR import *
from ..filters import *

def base_argparse(tmpdir):
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
    args.zip = False
    return args


class EmptyLocInfo:
    def __getitem__(self, key):
        # return something that can be incremented
        return 0

    def __setitem__(self, key, value):
        pass

class VCFRec:
    def __init__(self):
        self.FILTER = ''

class DummyRecBase:
    def __init__(self):
        self.vcfrecord = VCFRec()
        self.info = {}
        self.format = {}
    def GetCalledSamples(self):
        return np.array([True, True, False])
    def GetNumSamples(self):
        return 3

def test_CallRateFilter(tmpdir):
    class TestCallRateRec(DummyRecBase):
        def GetCallRate(self):
            return 0.5

    args = base_argparse(tmpdir)
    args.min_locus_callrate = 0.7
    locus_filters = BuildLocusFilters(args)
    assert ApplyLocusFilters(
        TestCallRateRec(), locus_filters, EmptyLocInfo(), False
    )

    args = base_argparse(tmpdir)
    args.min_locus_callrate = 0.3
    locus_filters = BuildLocusFilters(args)
    assert not ApplyLocusFilters(
        TestCallRateRec(), locus_filters, EmptyLocInfo(), False
    )


def test_HWEFilter(tmpdir):
    class TestHWERec(DummyRecBase):
        def __init__(self):
            super().__init__()
        def GetGenotypeCounts(self, uselength=False):
            if not uselength:
                return {
                    ('ATATAT', 'ATATAT'): 2,
                    ('ATATAT', 'ATAAAT'):  2,
                    ('ATATAT', 'ATATATAT'):  1,
                    ('ATAAAT', 'ATAAAT'):  2,
                    ('ATAAAT', 'ATATATAT'):  1,
                    ('ATATATAT', 'ATATATAT'):  2
                }
                # acutal het = 0.4
                # exp het : 1-2*(7/20)**2-(6/20)**2 = 0.665
                # p = 0.95
            else:
                return {
                    (3, 3):  6,
                    (3, 4):  2,
                    (4, 4):  2
                }
                # actual het = 0.2
                # exp het : 1 - .7**2 - .3**2 = 0.42
                # p = 0.21
        def GetAlleleFreqs(self, uselength=False):
            if not uselength:
                return {'ATATAT': .35, 'ATAAAT': .35, 'ATATATAT': .3}
            else:
                return {3: .7, 4: .3}

    def run_test(thresh, passes, uselength=False):
        args = base_argparse(tmpdir)
        args.min_locus_hwep = thresh
        args.use_length =  uselength
        locus_filters = BuildLocusFilters(args)
        assert passes != ApplyLocusFilters(
            TestHWERec(), locus_filters, EmptyLocInfo(), False
        )

    run_test(0.05, True, uselength = True)
    run_test(0.1, True, uselength = True)
    run_test(0.3, False, uselength = True)
    run_test(0.05, True)
    run_test(0.1, False)
    run_test(0.3, False)

def test_HETFilter(tmpdir):
    class TestHETRec(DummyRecBase):
        def __init__(self, count_3_1, count_3_2, count_4, count_5):
            super().__init__()
            self.count_3_1  = count_3_1
            self.count_3_2  = count_3_2
            self.count_4  = count_4
            self.count_5  = count_5
        def GetAlleleFreqs(self, uselength=False):
            total = (self.count_3_1 + self.count_3_2 + self.count_4 +
                self.count_5)
            if not uselength:
                return {
                    'ATATAT': self.count_3_1/total,
                    'ATAAAT': self.count_3_2/total,
                    'ATATATAT': self.count_4/total,
                    'ATATATATAT': self.count_5/total
                }
            else:
                return {
                    3: (self.count_3_1 + self.count_3_2)/total,
                    4: self.count_4/total,
                    5: self.count_5/total
                }


    def run_test(freq_list, thresh, higher, uselength=False):
        args = base_argparse(tmpdir)
        args.min_locus_het = thresh
        args.use_length = uselength
        locus_filters = BuildLocusFilters(args)
        assert higher != ApplyLocusFilters(
            TestHETRec(*freq_list), locus_filters, EmptyLocInfo(), False
        )

        args = base_argparse(tmpdir)
        args.max_locus_het = thresh
        args.use_length = uselength
        locus_filters = BuildLocusFilters(args)
        assert higher == ApplyLocusFilters(
            TestHETRec(*freq_list), locus_filters, EmptyLocInfo(), False
        )

    run_test([0.25, 0.25, 0.25, 0.25], 0.7, True)
    run_test([0.25, 0.25, 0.25, 0.25], 0.7, False, uselength=True)
    run_test([0.25, 0.25, 0.25, 0.25], 0.8, False)


def test_RegionFilter(tmpdir, vcfdir):
    class TestRegionRec(DummyRecBase):
        def __init__(self, chrom, pos):
            super().__init__()
            self.pos = pos
            self.chrom = chrom
            self.ref_allele_length = 10
    args = base_argparse(tmpdir)
    args.filter_regions = (vcfdir + '/dumpSTR_vcfs/sample_region.bed.gz,'
        + vcfdir + '/dumpSTR_vcfs/sample_region2.bed.gz')
    args.filter_regions_names = 'foo,bar'
    locus_filters = BuildLocusFilters(args)

    test_foo = TestRegionRec('chr21', 9487191)
    assert ApplyLocusFilters(test_foo, locus_filters, EmptyLocInfo(), False)
    assert test_foo.vcfrecord.FILTER == 'foo'

    test_not = TestRegionRec('chr21', 9487171)
    assert not ApplyLocusFilters(test_not, locus_filters, EmptyLocInfo(), False)
    assert test_not.vcfrecord.FILTER == 'PASS'

    test_both = TestRegionRec('chr21', 9487291)
    assert ApplyLocusFilters(test_both, locus_filters, EmptyLocInfo(), False)
    assert test_both.vcfrecord.FILTER == 'foo;bar'

    test_bar1 = TestRegionRec('chr20', 30)
    assert ApplyLocusFilters(test_bar1, locus_filters, EmptyLocInfo(), False)
    assert test_bar1.vcfrecord.FILTER == 'bar'

    test_bar2 = TestRegionRec('chr20', 230)
    assert ApplyLocusFilters(test_bar2, locus_filters, EmptyLocInfo(), False)
    assert test_bar2.vcfrecord.FILTER == 'bar'

    test_not_bar = TestRegionRec('chr20', 130)
    assert not ApplyLocusFilters(test_not_bar, locus_filters, EmptyLocInfo(), False)
    assert test_not_bar.vcfrecord.FILTER == 'PASS'


def test_HRUNFilter(tmpdir):
    class TestHRUNRec(DummyRecBase):
        def __init__(self, ref, period, full=None):
            super().__init__()
            self.ref_allele = ref
            if full is not None:
                self.full_alleles = (full, None)
            self.full = full is not None
            self.full = full
            self.info['PERIOD'] = period
        def HasFullStringGenotypes(self):
            return self.full

    args = base_argparse(tmpdir)
    args.filter_hrun = True
    locus_filters = BuildLocusFilters(args)
    for bp in {'A', 'T', 'G', 'C'}:
        assert ApplyLocusFilters(
            TestHRUNRec(bp*5, 5), locus_filters, EmptyLocInfo(), False
        )
        assert not ApplyLocusFilters(
            TestHRUNRec(bp*5, 6), locus_filters, EmptyLocInfo(), False
        )
        assert ApplyLocusFilters(
            TestHRUNRec(bp*6, 6), locus_filters, EmptyLocInfo(), False
        )
    assert not ApplyLocusFilters(
        TestHRUNRec('TTTTATTTT', 5), locus_filters, EmptyLocInfo(), False
    )
    assert ApplyLocusFilters(
        TestHRUNRec('ATTTTATTTTATTTTATTTTTATTTTATTTTATTTT', 5), locus_filters,
        EmptyLocInfo(), False
    )
    assert ApplyLocusFilters(
        TestHRUNRec('TTTTATTTTATTTTA', 5, full='TTTTTATTTTATTTTA'),
        locus_filters, EmptyLocInfo(), False
    )


def test_HipstrMaxCallFlankIndel(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['DFLANKINDEL'] = np.array([10, 5, np.nan]).reshape(-1, 1)
            self.format['DP'] = np.array([20, 20, np.nan]).reshape(-1,1)

    args = base_argparse(tmpdir)
    args.hipstr_max_call_flank_indel = 0.4
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    print(out)
    assert out[0] == pytest.approx(0.5)
    assert np.isnan(out[1])
    assert np.isnan(out[2]) # don't apply filter to nocalls


def test_HipstrMaxCallStutter(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['DSTUTTER'] = np.array([10, 5, np.nan]).reshape(-1, 1)
            self.format['DP'] = np.array([20, 20, np.nan]).reshape(-1,1)

    args = base_argparse(tmpdir)
    args.hipstr_max_call_stutter = 0.4
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert out[0] == pytest.approx(0.5)
    assert np.isnan(out[1])
    assert np.isnan(out[2]) # don't apply filter to nocalls


def test_HipstrMinSuppReads(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['ALLREADS'] = np.array([
                '0|23;1|123;2|5', '0|15;1|23;2|7',
                '0|23;1|444;2|12', '0|23;1|32;2|66',
                '0|867;1|23;2|13', '0|848;1|92;2|483',
                '', '', '.'])
            self.format['GB'] = np.array(['1|1', '1|1', '1|2', '2|1', '2|0',
                                          '0|2', '1|1', '0|0', '1|0'])
        def GetNumSamples(self):
            return 9
        def GetCalledSamples(self):
            return np.array([True, True, True, True, True, True, True, False, False])

    args = base_argparse(tmpdir)
    args.hipstr_min_supp_reads = 50
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert np.isnan(out[0])
    assert out[1] == 23
    assert out[2] == 12
    assert out[3] == 32
    assert out[4] == 13
    assert np.isnan(out[5])
    assert out[6] == 0 # If ALLREADS is missing, filter
    assert np.isnan(out[7]) # don't apply filter to nocalls
    assert np.isnan(out[8]) # don't apply filter to nocalls

def test_HipstrMinSuppReads_no_called_samples_with_reads(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['ALLREADS'] = np.array([
                '0|23;1|123;2|5', '0|15;1|23;2|7',
                '0|23;1|444;2|12', '0|23;1|32;2|66',
                '0|867;1|23;2|13', '0|848;1|92;2|483',
                '', '', '.'])
            self.format['GB'] = np.array(['1|1', '1|1', '1|2', '2|1', '2|0',
                                          '0|2', '1|1', '0|0', '1|0'])
        def GetNumSamples(self):
            return 9
        def GetCalledSamples(self):
            return np.array([False, False, False, False, False, False, True, False, True])

    args = base_argparse(tmpdir)
    args.hipstr_min_supp_reads = 50
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert out.shape == (9,)
    assert np.all(out[[6,8]] == 0)
    assert np.all(np.isnan(out[[0,1,2,3,4,5,7]]))

def test_HipstrDP(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['DP'] = np.array([10,20, np.nan]).reshape(-1,1)

    args = base_argparse(tmpdir)
    args.hipstr_min_call_DP = 15
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert out[0] == 10
    assert np.isnan(out[1])
    assert np.isnan(out[2]) # don't apply filter to nocalls

    args = base_argparse(tmpdir)
    args.hipstr_max_call_DP = 15
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert out[1] == 20
    assert np.isnan(out[0])
    assert np.isnan(out[2]) # don't apply filter to nocalls


def test_HipstrMinCallQ(tmpdir):
    class TestRec(DummyRecBase):
        def __init__(self):
            super().__init__()
            self.format['Q'] = np.array([.5, .9, np.nan]).reshape(-1, 1)

    args = base_argparse(tmpdir)
    args.hipstr_min_call_Q = 0.6
    call_filters = BuildCallFilters(args)
    assert len(call_filters) == 1
    out = call_filters[0](TestRec())
    assert out[0] == pytest.approx(0.5)
    assert np.isnan(out[1])
    assert np.isnan(out[2]) # don't apply filter to nocalls

'''
Test needs to be updated from old pyvcf
syntax to new cyvcf2 syntax
def test_GangSTRFilters(tmpdir, vcfdir): 
    args = base_argparse(tmpdir) 
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

    ## Test2: call filter: GangSTRCallMinDepth
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.gangstr_min_call_DP = 50
    call_filters = BuildCallFilters(args)
    reasons2 = FilterCall(vcfcall, call_filters)
    assert reasons2==['GangSTRCallMinDepth']
    
    ## Test3: call filter: GangSTRCallMaxDepth
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.gangstr_max_call_DP = 10
    call_filters = BuildCallFilters(args)
    reasons3 = FilterCall(vcfcall, call_filters)
    assert reasons3==['GangSTRCallMaxDepth']

    ## Test4: call filter: GangSTRCallMinQ
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.gangstr_min_call_Q = 1
    call_filters = BuildCallFilters(args)
    reasons4 = FilterCall(vcfcall, call_filters)
    assert reasons4==['GangSTRCallMinQ']

    ## Test5: call filter: GangSTRCallExpansionProbHom
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.gangstr_expansion_prob_hom = 1
    call_filters = BuildCallFilters(args)
    reasons5 = FilterCall(vcfcall, call_filters)
    assert reasons5==['GangSTRCallExpansionProbHom']

    ## Test6: call filter: ProbHet
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.gangstr_expansion_prob_het = 0.8
    call_filters = BuildCallFilters(args)
    reasons6 = FilterCall(vcfcall, call_filters)
    assert reasons6==['GangSTRCallExpansionProbHet']

    args = base_argparse(tmpdir)
    args.gangstr_expansion_prob_het = 0.0
    call_filters = BuildCallFilters(args)
    reasons61 = FilterCall(vcfcall, call_filters)
    assert reasons61==[]

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

    ## Test7: call filter: GangSTRCallExpansionProTotal
    vcfcall = call2
    args = base_argparse(tmpdir)
    args.gangstr_expansion_prob_total = 1
    call_filters = BuildCallFilters(args)
    reasons7 = FilterCall(vcfcall, call_filters)
    assert reasons7==['GangSTRCallExpansionProbTotal']

    ## Test8: call filter: filter span only
    vcfcall = call1
    args = base_argparse(tmpdir)
    args.filter_span_only = True
    call_filters = BuildCallFilters(args)
    reasons71 = FilterCall(vcfcall, call_filters)
    assert reasons71==[]

    vcfcall = call2
    args = base_argparse(tmpdir)
    args.gangstr_filter_span_only = True
    call_filters = BuildCallFilters(args)
    reasons72 = FilterCall(vcfcall, call_filters)
    assert reasons72==['GangSTRCallSpanOnly']

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
    args = base_argparse(tmpdir)
    args.gangstr_filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons81 = FilterCall(vcfcall, call_filters)
    assert reasons81==[]

    vcfcall = call2
    args = base_argparse(tmpdir)
    args.gangstr_filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons82 = FilterCall(vcfcall, call_filters)
    assert reasons82==['GangSTRCallSpanBoundOnly']

    vcfcall = call3    
    args = base_argparse(tmpdir)
    args.gangstr_filter_spanbound_only = True
    call_filters = BuildCallFilters(args)
    reasons83 = FilterCall(vcfcall, call_filters)
    assert reasons83==['GangSTRCallSpanBoundOnly']

    # Line 4 of artificial vcf
    record4 = next(artificial_vcf)
    call4 = record4.samples[0]
    assert record4.POS==3011925
    
    # Check min supporting reads
    # long compared to read length
    vcfcall = call4
    args = base_argparse(tmpdir)
    args.gangstr_require_support = 2
    args.gangstr_readlen = 20
    call_filters = BuildCallFilters(args)
    reason84 = FilterCall(vcfcall, call_filters)
    assert reason84==['GangSTRCallRequireSupport']

    # intermediate range compared to read length
    vcfcall = call4
    args = base_argparse(tmpdir)
    args.gangstr_require_support = 2
    args.gangstr_readlen = 400
    call_filters = BuildCallFilters(args)
    reason85 = FilterCall(vcfcall, call_filters)
    assert reason85==['GangSTRCallRequireSupport']


    # should see enclosing for everything. ultra-long read
    vcfcall = call4
    args = base_argparse(tmpdir)
    args.gangstr_require_support = 2
    args.gangstr_readlen = 10000
    call_filters = BuildCallFilters(args)
    reason85 = FilterCall(vcfcall, call_filters)
    assert reason85==['GangSTRCallRequireSupport']
'''
