import argparse
import gzip
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

def test_BadPreexistingFields(args, testDumpSTRdir, capsys):
    fname = os.path.join(testDumpSTRdir, "bad_preexisting_hrun.vcf")
    args.vcf = fname
    retcode = main(args)
    captured = capsys.readouterr()
    assert "HRUN" in captured.err

    fname = os.path.join(testDumpSTRdir, "bad_preexisting_het_hwep.vcf")
    args.vcf = fname
    retcode = main(args)
    captured = capsys.readouterr()
    assert "HWEP" in captured.err and "HET" in captured.err

    fname = os.path.join(testDumpSTRdir, "bad_preexisting_filter_ac_refac.vcf")
    args.vcf = fname
    retcode = main(args)
    captured = capsys.readouterr()
    assert ("FILTER" in captured.err and "AC" in captured.err
            and "REFAC" in captured.err)

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
    args.gangstr_readlen = 100
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

def _make_info_dict(info):
    d = {}
    for pair in info.split(';'):
        k, v = pair.split('=')
        vals = np.array(v.split(','))
        try:
            vals = vals.astype(float)
        except ValueError: # not a numeric field
            pass
        d[k] = vals
    return d

def _make_format_list(fmt):
    l = []
    for val in fmt.split(':'):
        vals = np.array(val.split(','))
        try:
            vals = vals.astype(float)
        except ValueError: # not a numeric field
            pass
        l.append(vals)
    return l


# fname1 should be the output file
# fname2 should be the control file
# allow reordering of header lines
def _assert_same_vcf(fname1, fname2, info_ignore = set(),
                     format_ignore = set()):
    open_fn = open
    if fname1[-3:] == '.gz':
        open_fn = gzip.open
    print(fname1, fname2)
    headers1 = set()
    headers2 = set()
    with open_fn(fname1, mode='rt') as file1, open_fn(fname2, mode='rt') as file2:
        iter1 = iter(file1)
        iter2 = iter(file2)
        failed = False
        while True:
            line1, line2 = _grab_line_for_assertion(iter1, iter2)
            if line1[0] != '#':
                raise ValueError('Output VCF header truncated abruptly')
            if line1[1] != '#' and line2[1] == '#':
                print('Output VCF header has fewer lines than the control header')
                failed = True
                headers2.add(line2)
                break
            if line1[1] == '#' and line2[1] != '#':
                print('Output VCF header has more lines than the control header')
                failed = True
                headers1.add(line1)
                break
            if line1[1] != '#' and line2[1] != '#':
                break # these are the sample lines
            headers1.add(line1)
            headers2.add(line2)

        # ignore command lines, they will be different
        found_command = False
        for line in headers1:
            if '##command-DumpSTR=' in line:
                found_command = True
                headers1.remove(line)
                break
        for line in headers2:
            if '##command-DumpSTR=' in line:
                headers2.remove(line)
                break

        # compare headers
        for line in headers1:
            if line not in headers2:
                print("Found header line '" + line + "' in the output vcf "
                      "but not in the control vcf.")
                failed = True
        for line in headers2:
            if line not in headers1:
                print("Found header line '" + line + "' in the control vcf "
                      "but not in the output vcf.")
                failed = True
        if failed:
            raise ValueError("VCF headers not identical. See output")

        if not found_command:
            raise ValueError('Output VCF does not have a command-DumpSTR'
                             ' line.')
        # compare sample lines
        if line1 != line2:
            raise ValueError('Output vcf sample line differs from control vcf'
                  ' sample line.\nSample line in output vcf'
                  ': ' + line1 + '\nSample line in control vcf: ' + line2)

        # compare loci
        iter1 = iter(file1)
        iter2 = iter(file2)
        linenum = 2 + len(headers1)
        format_ignore_idxs = set()
        while True:
            linenum += 1
            lines = _grab_line_for_assertion(iter1, iter2)
            if lines is None:
                return
            line1, line2 = lines
            for idx in range(len(line1)):
                if idx <= 6 or idx == 8:
                    if idx == 8:
                        fmt = line1[idx].split(':')
                        for val in format_ignore:
                            format_ignore_idxs.add(fmt.index(val))
                    if line1[idx] == line2[idx]:
                        continue
                    if idx == 0:
                        field_name = 'CHROM'
                    if idx == 1:
                        field_name = 'POS'
                    if idx == 2:
                        field_name = 'ID'
                    if idx == 3:
                        field_name = 'REF'
                    if idx == 4:
                        field_name = 'ALT'
                    if idx == 5:
                        field_name = 'QUAL'
                    if idx == 6:
                        field_name = 'FILTER'
                    if idx == 8:
                        field_name = 'FORMAT'
                    raise ValueError('Output file differs from control file'
                                     ' at line ' + str(linenum) + ' at field '
                                     + field_name + '\nOutput line: ' + line1[idx] +
                                     '\nControl line: ' + line2[idx])
                elif idx == 7:
                    # INFO field, allow permuations and changes from int to
                    # double
                    info1 = _make_info_dict(line1[7])
                    info2 = _make_info_dict(line2[7])
                    if info1.keys() != info2.keys():
                        raise ValueError(
                            'Output file differs from control file'
                             ' at line ' + str(linenum) + ' where they '
                            'have different INFO fields' + '\nOutput: ' +
                            line1[7] + '\nControl: ' + line2[7])
                    for k in info1:
                        if k in info_ignore:
                            continue
                        if not np.all(info1[k] == info2[k]):
                            raise ValueError(
                                'Output file differs from control file'
                                ' at line ' + str(linenum) + ' at INFO '
                                'field ' + k + '\nOutput: ' + str(info1[k]) +
                                '\nControl: ' + str(info2[k]))
                elif idx > 8:
                    # FORMAT field, allow changes from '.' to '.,.' and
                    # allow changes from int to double
                    format1 = _make_format_list(line1[idx])
                    format2 = _make_format_list(line2[idx])
                    sample_num = str(idx - 8)
                    if len(format1) != len(format2):
                        raise ValueError(
                            'Output file differs from control file'
                             ' at line ' + str(linenum) + ' where they '
                             'have different numbers of fields for sample #' +
                             sample_num + '\n' + line1[idx] + '\n' + line2[idx])
                    for count, (f1, f2) in enumerate(zip(format1, format2)):
                        if count in format_ignore_idxs:
                            continue
                        if (f1.dtype.kind == 'U' and
                                np.all(f1 == '.') and
                                np.all(f2 == '.')):
                            continue
                        if not np.all(f1 == f2):
                            raise ValueError(
                                'Output file differs from control file'
                                ' at line ' + str(linenum) + ' at sample #'
                                + sample_num + ' at field ' + str(count + 1) +
                                '\nOutput: ' + str(f1) + '\nControl: ' + str(f2))


# fname1 should be the output file
# fname2 should be the control file
# files should not be gzipped
def _assert_same_file(fname1, fname2, simple_name):
    print(fname1, fname2)
    with open(fname1) as file1, open(fname2) as file2:
        iter1 = iter(file1)
        iter2 = iter(file2)
        linenum = 0
        while True:
            linenum += 1
            lines = _grab_line_for_assertion(iter1, iter2)
            if lines is None:
                return
            
            line1, line2 = lines
            if line1 != line2:
                raise ValueError(
                    'Output ' + simple_name + ' file differs from control file'
                    ' at line ' + str(linenum) + '.\nLine in output'
                    ' file: ' + line1 + '\nLine in control file: ' + line2)

# return a pair of lines
# or raise an error if one iterator returned before the other
# iter1 should be the output file
# iter2 should be the control file
def _grab_line_for_assertion(iter1, iter2):
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
                'Output file has fewer lines than control file. '
            )
        else:
            raise ValueError(
                'Output file has more lines than control file. '
            )
    if file1ended and file2ended:
        return None

    return line1.strip(), line2.strip()


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
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/locus_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/locus_filters' + ext,
                          ext)

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
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/drop_filtered.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/locus_filters' + ext,
                          ext)


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
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/advntr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/advntr_filters' + ext,
                          ext)


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
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/hipstr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'},
                     format_ignore= {'GLDIFF'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/hipstr_filters' + ext,
                          ext)


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
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/gangstr_filters_most.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/gangstr_filters_most' + ext,
                          ext)

def test_output_gangstr_expansion_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/test_gangstr.vcf.gz'
    args.gangstr_expansion_prob_het = 0.001
    args.gangstr_expansion_prob_hom = 0.0005
    args.gangstr_expansion_prob_total =  0.001

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/gangstr_filters_expansion.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/gangstr_filters_expansion' + ext,
                          ext)


def test_output_popstr_filters(args, testDumpSTRdir):
    args.vcf = testDumpSTRdir + '/NA12878_chr21_popstr.sorted.vcf.gz'
    args.popstr_min_call_DP = 30
    args.popstr_max_call_DP = 200
    args.popstr_require_support = 15
    args.use_length = True

    assert main(args) == 0
    # expect changes in precision for HET and HWEP
    # that will make them too much of a pain to compare
    _assert_same_vcf(args.out + '.vcf',
                     testDumpSTRdir + '/popstr_filters.vcf',
                     info_ignore = {'AC', 'REFAC', 'HET', 'HWEP'})
    for ext in '.samplog.tab', '.loclog.tab':
        _assert_same_file(args.out + ext,
                          testDumpSTRdir + '/popstr_filters' + ext,
                          ext)



