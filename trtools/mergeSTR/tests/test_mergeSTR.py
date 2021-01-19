import argparse
import os
import re

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
    args.trim = False
    return args


@pytest.fixture
def mrgvcfdir(vcfdir):
	return os.path.join(vcfdir, "mergeSTR_vcfs")


# Set up dummy class
class DummyRecord:
    def __init__(self, chrom, pos, ref, alts=[], info = {}):
        self.chrom = chrom
        self.pos = pos
        self.info = info
        self.end_pos = pos + len(ref) - 1
        self.vcfrecord = argparse.Namespace()
        self.vcfrecord.REF = ref
        self.vcfrecord.ALT = alts

    def __str__(self):
        return self.pos

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
    assert "Different contigs (or contig orderings) found across VCF files." in str(info.value)

def test_DifferentContigLengths(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_contigdifflength.vcf.gz")
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "Different contigs (or contig orderings) found across VCF files." in str(info.value)

def test_SameContigsDifferentOrder(args, vcfdir, mrgvcfdir):
    fname1 = os.path.join(vcfdir, "one_sample_multiple_chroms.vcf.gz")
    fname2 = os.path.join(
        mrgvcfdir,
        "one_sample_multiple_chroms_diff_contig_order.vcf.gz"
    )
    args.vcfs = fname1 + "," + fname2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "Different contigs (or contig orderings) found across VCF files." in str(info.value)

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

def test_ConflictingRefs():
    # Set up dummy records
    dummy_records = [] 
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG'))
    dummy_records.append(DummyRecord('chr1', 100, 'CAG'))

    retval = GetRefAllele(dummy_records, [True, True, True])
    assert retval is None

    retval = GetRefAllele(dummy_records, [True, True, False])
    assert retval == "CAGCAG"


def test_conflicting_refs_merge(capsys, args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "diff_ref2.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "diff_ref.vcf.gz")
    # confirm that output is independent of file order
    for fnames in fname1 + ',' + fname2, fname2 + ',' + fname1:
        args.vcfs = fnames
        main(args)

        captured = capsys.readouterr().err

        # looking for four skipped records
        # the first has different start coords
        # the second has different end coords
        # the third has smaller different start and end coords
        # the fourth has a greater start and smaller end coord
        assert re.search('3029300.*overlaps', captured) is not None
        assert re.search('3068357.*overlaps', captured) is not None
        assert re.search('3074665.*overlaps', captured) is not None
        assert re.search('3086956.*overlaps', captured) is not None
        # looking for one merged record - same refs
        vcf = cyvcf2.VCF(args.out + ".vcf")
        var = next(vcf)
        assert var.POS == 3045469
        with pytest.raises(StopIteration):
            next(vcf)

def test_conflicting_refs_same_file():
    recs1 = [DummyRecord('chr1', 100, 'CAG', 'CAGCAG'),
             DummyRecord('chr1', 100, 'CAG', 'CAGCAGCAG')]
    recs2 = [DummyRecord('chr1', 100, 'CAG', 'CAGCAG'),
             DummyRecord('chr1', 1000, 'CAG', 'CAGCAGCAG')]
    for test_count, readers in enumerate([(iter(recs1), iter(recs2)),
                                          (iter(recs2), iter(recs1))]):
        itr = RecordIterator(readers, ['chr1'])

        recs, mins = next(itr)
        if test_count == 1:
            # unreverse
            recs = [recs[1], recs[0]]
            mins = [mins[1], mins[0]]
        assert mins == [False, True]
        assert recs[1].pos == 1000
        with pytest.raises(StopIteration):
            next(itr)


def test_conflicting_refs_threeway():
    recs1 = [DummyRecord('chr1', 100, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 103, 'CAGCAG'),
             DummyRecord('chr1', 124, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 140, 'CAGCAGCAGCAG')]
    recs2 = [DummyRecord('chr1', 109, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 160, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 180, 'CAGCAGCAGCAG')]
    recs3 = [DummyRecord('chr1', 115, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 140, 'CAGCAGCAGCAG'),
             DummyRecord('chr1', 189, 'CAGCAGCAGCAG')]

    def order(l, idxs):
        new_l = []
        new_l.append(l[int(np.where(np.array(idxs) == 0)[0])])
        new_l.append(l[int(np.where(np.array(idxs) == 1)[0])])
        new_l.append(l[int(np.where(np.array(idxs) == 2)[0])])
        return new_l

    list_of_recs = [recs1, recs2, recs3]
    # repeat this test for any possible ordering of the three files
    # then reorder the results to be in the same order as displayed
    # above
    # only the 140 and 160 records should be returned, everything else
    # overlaps something else
    for idxs in [
        (0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)
    ]:
        readers = [iter(list_of_recs[idxs[0]]),
                   iter(list_of_recs[idxs[1]]),
                   iter(list_of_recs[idxs[2]])]
        itr = RecordIterator(readers, ['chr1'])

        recs, mins = next(itr)
        recs = order(recs, idxs)
        mins = order(mins, idxs)
        assert mins == [True, False, True]
        assert recs[0].pos == 140
        assert recs[2].pos == 140

        recs, mins = next(itr)
        recs = order(recs, idxs)
        mins = order(mins, idxs)
        assert mins == [False, True, False]
        assert recs[1].pos == 160

        with pytest.raises(StopIteration):
            next(itr)


def test_GetInfoItem(capsys):
    # Set up dummy records
    dummy_records = []
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 120}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={'END': 110}))
    dummy_records.append(DummyRecord('chr1', 100, 'CAGCAG', info={}))

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

def _check_merged_output(outputvcffile, inputvcffiles, vcftype):
    # currently doesn't handle changing chrom
    # doesn't support changed sample names
    # don't confirm that the choice of which records
    # have been merged is correct
    # Does confirm: for each merged record, confirm that all the data
    # transferred over is the same as before for each sample
    outputvcf = None
    inputvcfs = []
    try:
        outputvcf = trh.TRRecordHarmonizer(cyvcf2.VCF(outputvcffile), vcftype)
        inputvcfs = [
            trh.TRRecordHarmonizer(cyvcf2.VCF(inputvcffile), vcftype) for
            inputvcffile in inputvcffiles
        ]
        out_samps = outputvcf.vcffile.samples
        in_samps = [inputvcf.vcffile.samples for inputvcf in inputvcfs]
        inrecs = [next(inputvcf) for inputvcf in inputvcfs]
        for rec in outputvcf:
            for idx in range(len(inputvcfs)):
                while inrecs[idx] is not None and inrecs[idx].pos < rec.pos:
                    try:
                        inrecs[idx] = next(inputvcfs[idx])
                    except StopIteration:
                        inrecs[idx] = None
            for idx, inrec in enumerate(inrecs):
                if inrec.pos != rec.pos:
                    continue
                for info in INFOFIELDS[vcftype]:
                    if not info[1]:
                        continue
                    assert inrec.info[info[0]] == rec.info[info[0]]
                for old_samp_idx, samp in enumerate(in_samps[idx]):
                    new_samp_idx = out_samps.index(samp)
                    assert (rec.GetCalledSamples()[new_samp_idx] ==
                            inrec.GetCalledSamples()[old_samp_idx])
                    if not rec.GetCalledSamples()[new_samp_idx]:
                        continue
                    assert np.all(inrec.GetStringGenotypes()[old_samp_idx, :] ==
                                  rec.GetStringGenotypes()[new_samp_idx, :])
                    for fmt in FORMATFIELDS[vcftype]:
                        infmt = inrec.format[fmt]
                        if len(infmt.shape) == 1:
                            assert np.all(infmt[old_samp_idx] ==
                                          rec.format[fmt][new_samp_idx])
                        else:
                            nans = np.isnan(infmt[old_samp_idx, :])
                            assert np.all(nans ==
                                          np.isnan(rec.format[fmt][new_samp_idx, :]))
                            assert np.all(infmt[old_samp_idx, ~nans] ==
                                          rec.format[fmt][new_samp_idx, ~nans])
    finally:
        if outputvcf is not None:
            outputvcf.close()
        for inputvcf in inputvcfs:
            if inputvcf is not None:
                inputvcf.close()


def test_advntr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_advntr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_advntr2.vcf.gz")
    args.vcftype = "advntr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/advntr_merged.vcf")
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/test_file_advntr1.vcf.gz",
         mrgvcfdir + "/test_file_advntr2.vcf.gz"],
        trh.VcfTypes.advntr
    )


def test_eh_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_eh1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_eh2.vcf.gz")
    args.vcftype = "eh"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/eh_merged.vcf")
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/test_file_eh1.vcf.gz",
         mrgvcfdir + "/test_file_eh2.vcf.gz"],
        trh.VcfTypes.eh
    )


def test_gangstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_gangstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_gangstr2.vcf.gz")
    args.vcftype = "gangstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/test_file_gangstr1.vcf.gz",
         mrgvcfdir + "/test_file_gangstr2.vcf.gz"],
        trh.VcfTypes.gangstr
    )
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/gangstr_merged.vcf")


def test_hipstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_hipstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_hipstr2.vcf.gz")
    args.vcftype = "hipstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/test_file_hipstr1.vcf.gz",
         mrgvcfdir + "/test_file_hipstr2.vcf.gz"],
        trh.VcfTypes.hipstr
    )
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/hipstr_merged.vcf")


def test_popstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "test_file_popstr1.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "test_file_popstr2.vcf.gz")
    args.vcftype = "popstr"
    args.vcfs = fname1 + "," + fname2
    assert main(args) == 0
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/test_file_popstr1.vcf.gz",
         mrgvcfdir + "/test_file_popstr2.vcf.gz"],
        trh.VcfTypes.popstr
    )
    assert_same_vcf(args.out + '.vcf', mrgvcfdir + "/popstr_merged.vcf")


def test_trimmed_hipstr_output(args, mrgvcfdir):
    fname1 = os.path.join(mrgvcfdir, "NA12878_chr21_hipstr.sorted.vcf.gz")
    fname2 = os.path.join(mrgvcfdir, "NA12891_chr21_hipstr.sorted.vcf.gz")
    fname3 = os.path.join(mrgvcfdir, "NA12892_chr21_hipstr.sorted.vcf.gz")
    args.vcftype = "hipstr"
    args.trim = True
    args.vcfs = fname1 + "," + fname2 + "," + fname3
    assert main(args) == 0
    _check_merged_output(
        args.out + '.vcf',
        [mrgvcfdir + "/NA12878_chr21_hipstr.sorted.vcf.gz",
         mrgvcfdir + "/NA12891_chr21_hipstr.sorted.vcf.gz",
         mrgvcfdir + "/NA12892_chr21_hipstr.sorted.vcf.gz"],
        trh.VcfTypes.hipstr
    )


# TODO questions and issues to confirm, test or  address:
# confirm conflicting samples cause this to fail unless --update-sample-from-file is given
# required info fields cause failure when a record is missing them but not when the entire VCF is missing them
# we silently return if there are no info or header fields
# we don't fail if info fields are different
# how should we merge info fields with reqd = false
# what if there are multiple records at the same location in the same VCF
# if a genotype is no called but there is other format info,
# we aren't emitting it. Is that intended?
