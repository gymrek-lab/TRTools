import argparse
import os, sys
import numpy as np
import pytest

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
    args.vcftype1 = None
    args.vcftype2 = None
    args.verbose = False
    args.noplot = False
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
    
def test_GetBubbleLegent():
    # only 3 values
    sample_counts = [1,1,1,2,2,3,3,3,3]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 2, 3]
    assert all([a == b for a, b in zip(actual, expected)])

    # More than three, linear plot (max(val)/min(val) < 10)
    sample_counts = [1,1,1,2,2,3,3,3,3,4,4,4,5,5]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 3, 5]
    assert all([a == b for a, b in zip(actual, expected)])


    # More than three, log plot (max(val)/min(val) > 10)
    sample_counts = [1,1,1,5,5,10,10,10,10,14,14,14,100,100]
    actual = GetBubbleLegend(sample_counts)
    expected = [1, 10, 100]
    assert all([a == b for a, b in zip(actual, expected)])
