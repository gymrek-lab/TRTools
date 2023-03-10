import argparse
import os

import numpy as np
import pandas as pd
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
    args.balanced_accuracy = False
    args.fraction_concordant_len_sum = False
    args.vcf2_beagle_dosages = False
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

def test_ten_sample_output(tmpdir, vcfdir):
    # also tests sample reordering and non equivalent sample lists
    # also tests missingness for one or the other or both samples
    # also tests samples with no shared calls (but does not test loci with no shared calls)
    # does not test phasedness
    # currently only tests the locus file output
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    VCF1 = os.path.join(vcfcomp, "test_ten_sample1.vcf.gz")
    VCF2 = os.path.join(vcfcomp, "test_ten_sample2.vcf.gz")
    args = base_argparse(tmpdir)
    args.region = None
    args.vcf1 = VCF1
    args.vcf2 = VCF2
    args.ignore_phasing = True
    args.balanced_accuracy = True
    args.fraction_concordant_len_sum = True
    retcode = main(args)
    assert retcode == 0

    assert open(tmpdir + "/test_compare-vcf1-omitted-samples.tab", "r").read().strip() == 'HG00117'
    assert open(tmpdir + "/test_compare-vcf2-omitted-samples.tab", "r").read().strip() == 'HG00119'

    locus = pd.read_csv(tmpdir + "/test_compare-locuscompare.tab", sep='\t')
    assert np.all(locus['start'] == [31720, 33450])
    assert np.allclose(locus['fraction_concordant_len'], [1/3, 4/9])
    assert np.allclose(locus['fraction_concordant_len_sum'], [4/9, 4/9])
    assert np.allclose(locus['balanced_accuracy'], [(1 + 3/8)/2, (3/4 + 1/4) / 3])
    assert eval(locus['len_sum_frequencies'][0]) == {28: 8/9, 27: 1/9}
    assert eval(locus['len_sum_frequencies'][1]) == {30: 4/9, 31: 4/9, 28: 1/9}
    assert eval(locus['len_sum_accuracies'][0]) == {28: 3/8, 27: 1}
    assert eval(locus['len_sum_accuracies'][1]) == {30: 3/4, 31: 1/4, 28: 0}
    assert np.allclose(locus['mean_absolute_difference'], [5/9, 2])
    assert np.allclose(locus['r'], [0.35, 0.3228208661506])
    assert np.all(locus['numcalls']  == [9, 9])
    assert np.all(locus['n_missing_only_vcf1']  == [0, 1])
    assert np.all(locus['n_missing_only_vcf2']  == [1, 0])
    assert np.all(locus['n_missing_both']  == [1, 1])

def dicts_close(d1, d2):
    keylist = d1.keys()
    if not set(keylist) == set(d2.keys()):
        return False
    array_1 = np.array([d1[key] for key in keylist])
    array_2 = np.array([d2[key] for key in keylist])

    return np.allclose(array_1, array_2)

def test_ten_sample_output_beagle(tmpdir, vcfdir):
    # also tests sample reordering and non equivalent sample lists
    # also tests missingness for one or the other or both samples
    # also tests samples with no shared calls (but does not test loci with no shared calls)
    # does not test phasedness
    # currently only tests the locus file output
    vcfcomp = os.path.join(vcfdir, "compareSTR_vcfs")
    VCF1 = os.path.join(vcfcomp, "test_ten_sample1.vcf.gz")
    VCF2 = os.path.join(vcfcomp, "test_ten_sample2.vcf.gz")
    args = base_argparse(tmpdir)
    args.region = None
    args.vcf1 = VCF1
    args.vcf2 = VCF2
    args.ignore_phasing = True
    args.balanced_accuracy = True
    args.fraction_concordant_len_sum = True
    args.vcf2_beagle_dosages = True
    args.samples = os.path.join(vcfcomp, "test_ten_sample_beagle_samples.txt")
    retcode = main(args)
    assert retcode == 0

    assert "test_compare-vcf1-omitted-samples.tab" in os.listdir(tmpdir)
    assert "test_compare-vcf2-omitted-samples.tab" not in os.listdir(tmpdir)
    assert open(tmpdir + "/test_compare-vcf1-omitted-samples.tab", "r").read().strip() == 'HG00117'

    locus = pd.read_csv(tmpdir + "/test_compare-locuscompare.tab", sep='\t')
    assert np.all(locus['start'] == [31720, 33450])
    assert np.allclose(locus['fraction_concordant_len'], [(1 + .1 + .08)/5, (.42 + 1 + .06 + .6)/5])
    assert np.allclose(locus['fraction_concordant_len_sum'], [(1+.32+.32)/5, (.42 + 1 + .06 + .6)/5])
    assert np.allclose(locus['balanced_accuracy'], [(1+.32+.32)/5,  ((.42 + .6)/2 + (1 + .06)/3) / 2])
    assert dicts_close(eval(locus['len_sum_frequencies'][0]), {28: 1})
    assert dicts_close(eval(locus['len_sum_frequencies'][1]), {30: 2/5, 31: 3/5})
    assert dicts_close(eval(locus['len_sum_accuracies'][0]), {28: (1+.32+.32)/5})
    assert dicts_close(eval(locus['len_sum_accuracies'][1]), {30: (.42 + .6)/2, 31: (1 + .06)/3})
    #assert np.allclose(locus['mean_absolute_difference'], [5/9, 2])
    #assert np.allclose(locus['r'], [0.35, 0.3228208661506])
    assert np.all(locus['numcalls']  == [5, 5])
    assert np.all(locus['n_missing_only_vcf1']  == [0, 1])
    assert np.all(locus['n_missing_only_vcf2']  == [1, 0])
    assert np.all(locus['n_missing_both']  == [1, 1])

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

