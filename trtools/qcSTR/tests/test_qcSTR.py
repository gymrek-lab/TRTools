import argparse
import os, sys
import numpy as np
import pytest

from ..qcSTR import * 


# Set up base argparser
def base_argparse(tmpdir):
    args = argparse.Namespace()
    args.vcf = None
    args.out = str(tmpdir / "test_qc")
    args.vcftype = "auto"
    args.samples = None
    args.numrecords = None
    args.period = None
    args.quality = []
    return args

# Just confirm that the method doesn't throw an error
def test_OutputDiffRefHistogram(tmpdir):
    diffs_from_ref = [0,0,0,0,1,0,-1,-2,-4,-5]
    fname = str(tmpdir / "test_hist.pdf")
    OutputDiffRefHistogram(diffs_from_ref, fname)

# Just confirm that the method doesn't throw an error
def test_OutputDiffRefBias(tmpdir):
    diffs_from_ref = [0,0,0,0,1,0,-1,-2,-4,-5]
    reflens = [1,2,3,4,5,6,7,8,9,10]
    fname = str(tmpdir / "test_bias.pdf") # will be empty because of low count <25
    OutputDiffRefBias(diffs_from_ref, reflens, fname)

# Just confirm that the method doesn't throw an error
def test_OutPutSampleCallrate(tmpdir):
    sample_calls = {'s1': 120, 's2': 10}
    fname = str(tmpdir / "test_qc1.pdf")
    OutputSampleCallrate(sample_calls, fname)

# Just confirm that the method doesn't throw an error
def test_OutPutChromCallrate(tmpdir):
    chrom_calls = {'chr1': 100, 'chr2': 200}
    fname = str(tmpdir / "test_qc2.pdf")
    OutputChromCallrate(chrom_calls, fname)

def test_output_location_is_directory(tmpdir, vcfdir, capsys):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")
    args.out = str(tmpdir) 
    retcode = main(args)
    assert "is a directory" in capsys.readouterr().err
    assert retcode == 1

def test_output_location_is_in_nonexistant_directory(tmpdir, vcfdir, capsys):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")
    args.out = str(tmpdir / 'nonexistant_dir' / 'some-file-prefix')
    retcode = main(args)
    assert "does not exist" in capsys.readouterr().err
    assert retcode == 1

def test_asking_for_qual_plot_fails_when_no_qual_info(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "test_ExpansionHunter.vcf")
    args.quality = ['per-locus']
    retcode = main(args)
    assert "doesn't have quality scores" in capsys.readouterr().err
    assert retcode == 1

def test_dont_make_qual_plot_when_no_qual_info(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "test_ExpansionHunter.vcf")
    retcode = main(args)
    assert not os.path.exists(args.out + "-quality.pdf")
    assert retcode == 0

def test_make_default_qual_plot_few_samples(tmpdir, vcfdir):
    pass

def test_make_default_qual_plot_many_samples(tmpdir, vcfdir):
    pass

def test_make_single_qual_plots_explicit(tmpdir, vcfdir):
    pass

def test_make_all_qual_plots(tmpdir, vcfdir):
    pass

def test_main(tmpdir, vcfdir):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    # correct vcf
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # vcf file with no contig
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr_nocontig.vcf")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0
    
    # Set sample list 
    args.samples = os.path.join(qcdir, "test_samplelist.txt")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # Trying to test line 186 but seems to not work. TODO update with something else
    # Set sample list with sample that doesn't exist in VCF (should skip it)
    args.samples = os.path.join(qcdir, "test_samplelist.txt")
    with pytest.warns(UserWarning, match="fabricated"):
        retcode = main(args)
    assert retcode == 0

    # Non existent vcf
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_non_exist.vcf")
    retcode = main(args)
    assert retcode == 1
