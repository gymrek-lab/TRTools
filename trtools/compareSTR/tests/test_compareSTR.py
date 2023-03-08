import argparse
import os
from typing import List

import numpy as np
import pytest
from trtools.utils.tests.test_mergeutils import DummyHarmonizedRecord
from ..compareSTR import *


# Set up base argparser
def base_argparse(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf1 = None
    args.vcf2 = None
    args.out = str(tmpdir / "test_compare")
    args.vcftype = "auto"
    args.samples = None
    args.numrecords = None
    args.period = None
    args.region = "chr1"
    args.stratify_file = 0
    args.stratify_fields = None
    args.stratify_binsizes = None
    args.vcftype1 = "auto"
    args.vcftype2 = "auto"
    args.verbose = False
    args.noplot = False
    args.ignore_phasing = False
    args.bubble_min = -5
    args.bubble_max = 5
    return args


def test_main(tmpdir, vcfdir):
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    GangSTR_VCF1 = os.path.join(vcfcomp, "test_gangstr1.vcf.gz")
    GangSTR_VCF2 = os.path.join(vcfcomp, "test_gangstr2.vcf.gz")
    HipSTR_VCF = os.path.join(vcfcomp, "test_hipstr.vcf.gz")
    PopSTR_VCF = os.path.join(vcfcomp, "test_popstr.vcf.gz")
    EH_VCF = os.path.join(vcfcomp, "test_eh.vcf.gz")

    # Two Gangstr VCFs
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    retcode = main(args)
    assert retcode == 0

    # Two GangSTR VCFs, types set to auto
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcf2 = GangSTR_VCF2
    retcode = main(args)
    assert retcode == 0
    
    # Use samples
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.samples = os.path.join(vcfcomp, 'sample_list.txt')
    retcode = main(args)
    assert retcode == 0

    # Empty sample list
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.samples = os.path.join(vcfcomp, 'empty_list.txt')
    retcode = main(args)
    assert retcode == 1

    # TODO Two HipSTR VCFs
    # TODO Two popSTR VCFs
    # TODO Two EH VCFs
    # TODO Two advntr VCFs
    # TODO EH and GangSTR

    # TODO all comparisons between two different tools

    #args.vcf1 = GangSTR_VCF1
    #args.vcftype1 = 'gangstr'
    #args.vcf2 = EH_VCF
    #args.vcftype2 = 'eh'
    #retcode = main(args)
    #assert retcode == 0


    # Testing custom stratification:
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.stratify_fields = 'DP'
    args.stratify_binsizes = '0:100:10'
    retcode = main(args)
    assert retcode == 0

    # Incorrect stratification formats
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.stratify_fields = 'DP,ML'
    args.stratify_binsizes = '0:100:10'
    with pytest.raises(ValueError) as info:
        main(args)
    assert "--stratify-formats must be same length as --stratify-binsizes" in str(info.value)
    
    # Incorrect stratify file value
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.stratify_file = 4
    retcode = main(args)
    assert retcode == 1


    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.stratify_fields = 'NONEXISTENTFIELD'
    args.stratify_binsizes = '0:100:10'
    args.stratify_file = 0
    with pytest.raises(ValueError) as info:
        main(args)
    assert "FORMAT field NONEXISTENTFIELD must be present in both VCFs if --stratify-file=0" in str(info.value)
    args.stratify_file = 1
    with pytest.raises(ValueError) as info:
        main(args)
    assert "FORMAT field NONEXISTENTFIELD must be present in --vcf1 if --stratify-file=1" in str(info.value)
    args.stratify_file = 2
    with pytest.raises(ValueError) as info:
        main(args)
    assert "FORMAT field NONEXISTENTFIELD must be present in --vcf2 if --stratify-file=2" in str(info.value)
     # correct stratification
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'
    args.stratify_fields = 'DP'
    args.stratify_binsizes = '0:100:10'
    args.stratify_file = 0
    retcode = main(args)
    assert retcode == 0

    args.stratify_file = 1
    retcode = main(args)
    assert retcode == 0
    args.stratify_file = 2
    retcode = main(args)
    assert retcode == 0


    # No shared samples
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = os.path.join(vcfcomp, "test_gangstr2_wrongsamp.vcf.gz")
    args.vcftype2 = 'gangstr'
    retcode = main(args)
    assert retcode == 1

def test_no_comparable_records(tmpdir, vcfdir, capfd):
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    hipstr_vcf_1 = os.path.join(vcfcomp, "test_no_comparable_records_1.vcf.gz")
    hipstr_vcf_2 = os.path.join(vcfcomp, "test_no_comparable_records_2.vcf.gz")

    args = base_argparse(tmpdir)
    args.vcf1 = hipstr_vcf_1
    args.region = ""
    args.vcf2 = hipstr_vcf_2

    ret = main(args)
    _, err = capfd.readouterr()
    assert ret == 1
    assert err == "No comparable records were found, exiting!\n"


def test_better_comparability_calculation(tmpdir, vcfdir, capfd):
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")

    test_vcf_1 = os.path.join(vcfcomp, "test_better_comparability_calculation_1.vcf.gz")
    test_vcf_2 = os.path.join(vcfcomp, "test_better_comparability_calculation_2.vcf.gz")

    args = base_argparse(tmpdir)
    args.vcf1 = test_vcf_1
    args.region = ""
    args.vcftype1 = 'hipstr'
    args.vcf2 = test_vcf_2
    args.vcftype2 = 'hipstr'
    retcode = main(args)
    assert retcode == 0
    with open(tmpdir + "/test_compare-locuscompare.tab", "r") as out_overall:
        lines = out_overall.readlines()
        # Two of the records wont be compared
        assert len(lines) == 2
        _, err = capfd.readouterr()
                        ## first output is about records that have the same starting position but different end pos
        assert err == ("Records STR_40 and STR_40 overlap:\n"
                        "STR_40: (112695, 112700)\n"
                        "STR_40: (112695, 112702),\n"
                        "but are NOT comparable!\n"
                       ## second is more general, they just sort of overlap each other
                        "Records STR_41 and STR_41 overlap:\n"
                        "STR_41: (113695, 113700)\n"
                        "STR_41: (113693, 113702),\n"
                        "but are NOT comparable!\n"
                       ## ends are the same but start positions are different 
                       "Records STR_42 and STR_42 overlap:\n"
                       "STR_42: (114695, 114700)\n"
                       "STR_42: (114693, 114700),\n"
                       "but are NOT comparable!\n"
                       )
def test_comparability_handler(tmpdir, vcfdir, capfd):

    ### Tests without arguments
    handler = handle_overlaps

    records = [None, None]
    chrom_idxs = [np.inf, np.inf]
    min_idx = np.inf

    assert not handler(records, chrom_idxs, min_idx)

    records = [DummyHarmonizedRecord("chr1", 10), None]
    min_idx = 0

    assert not handler(records, chrom_idxs, min_idx)

    records = [None, DummyHarmonizedRecord("chr1", 10, 4, "AC")]
    chrom_idxs = [np.inf, 0]
    assert not handler(records, chrom_idxs, min_idx)

    records = [DummyHarmonizedRecord("chr2", 10, 4, "AC", end_pos=17), DummyHarmonizedRecord("chr1", 10, 4, "AC", end_pos=17)]
    chrom_idxs = [1, 0]
    # records from different chromosomes aren't comparable
    assert not handler(records, chrom_idxs, min_idx)

    chrom_idxs = [0, 0]
    assert handler(records, chrom_idxs, min_idx)

    records = [DummyHarmonizedRecord("chr1", 10, 5, "AC", "rec1", end_pos=19), DummyHarmonizedRecord("chr1", 10, 4, "AC", "rec2", end_pos=17)]
    assert not handler(records, chrom_idxs, min_idx)
    _, err = capfd.readouterr()
    assert err != ""

    records = [DummyHarmonizedRecord("chr1", 10, 4, "AC", end_pos=17), DummyHarmonizedRecord("chr1", 10, 4, "TG", end_pos=17)]
    assert handler(records, chrom_idxs, min_idx)

    records = [DummyHarmonizedRecord("chr1", 8, 5, "AC", end_pos=17), DummyHarmonizedRecord("chr1", 10, 4, "AC", end_pos=17)]
    assert not handler(records, chrom_idxs, min_idx)


def test_hipstr_position_harmonisation(tmpdir, vcfdir):
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")

    hipstr_vcf_1 = os.path.join(vcfcomp, "test_hipstr_flanking_bp_flanking.vcf.gz")
    hipstr_vcf_2 = os.path.join(vcfcomp, "test_hipstr_flanking_bp_non_flanking.vcf.gz")

    args = base_argparse(tmpdir)
    args.vcf1 = hipstr_vcf_1
    args.region = ""
    args.vcftype1 = 'hipstr'
    args.vcf2 = hipstr_vcf_2
    args.vcftype2 = 'hipstr'
    retcode = main(args)
    assert retcode == 0


    with open(tmpdir + "/test_compare-locuscompare.tab", "r") as out_overall:
        lines = out_overall.readlines()
        ## vcf1 : flank bp at start of record, vcf2: no flanking bp
        assert lines[1] == "1	101675	1.0	1.0	1\n"
        ## vcf1 : no flanking bp, vcf2: no flanking bp
        assert lines[2] == "1	111675	1.0	1.0	1\n"
        ## vcf1 : flanking bp at the end, vcf2: no flanking bp
        assert lines[3] == "1	112655	1.0	1.0	1\n"
        ## vcf1 : flanking bp at both sides, vcf2: no flanking bp
        assert lines[4] == "1	125557	1.0	1.0	1\n"

def test_wrong_vcftype(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    args.vcf1 = os.path.join(vcfcomp, "test_gangstr1.vcf.gz")
    args.vcf2 = os.path.join(vcfcomp, "test_gangstr2.vcf.gz")

    args.vcftype1 = 'eh'
    args.vcftype2 = 'gangstr'
    retcode = main(args)
    assert retcode == 1
    assert 'not one of those types' in capsys.readouterr().err

    args.vcftype1 = 'gangstr'
    args.vcftype2 = 'eh'
    retcode = main(args)
    assert retcode == 1
    assert 'not one of those types' in capsys.readouterr().err

def test_region(tmpdir, vcfdir, capsys):
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    GangSTR_VCF1 = os.path.join(vcfcomp, "test_gangstr1.vcf.gz")
    GangSTR_VCF2 = os.path.join(vcfcomp, "test_gangstr2.vcf.gz")
    args = base_argparse(tmpdir)
    args.vcf1 = GangSTR_VCF1
    args.vcftype1 = 'gangstr'
    args.vcf2 = GangSTR_VCF2
    args.vcftype2 = 'gangstr'

    # test correct region strings
    args.region = 'chr1'
    retcode = main(args)
    assert retcode == 0

    args.region = 'chr1:5000000000-'
    retcode = main(args)
    assert retcode == 0

    args.region = 'chr1:29-42'
    retcode = main(args)
    assert retcode == 0

    # test incorrect region strings
    args.region = '1'
    retcode = main(args)
    assert retcode == 1

    args.region = '1:-42'
    retcode = main(args)
    assert retcode == 1

def test_GetBubbleLegend():
    # only 3 values
    sample_counts = [1,2,3]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 2, 3]
    assert all([a == b for a, b in zip(actual, expected)])

    # More than three, linear plot (max(val)/min(val) < 10)
    sample_counts = [1,2,3,4,5]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 3, 5]
    assert all([a == b for a, b in zip(actual, expected)])

    # More than three, log plot (max(val)/min(val) > 10)
    sample_counts = [1,5,10,14,100]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 10, 100]
    assert all([a == b for a, b in zip(actual, expected)])

