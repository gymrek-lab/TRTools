import argparse
import glob, os, shutil, sys
import numpy as np
import pytest

from ..qcSTR import *
from ..qcSTR import _QualityTypes


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
    args.quality_ignore_no_call = False
    args.refbias_metric = "mean"
    args.refbias_mingts = 100
    args.refbias_xrange_min = 0
    args.refbias_xrange_max = 100
    args.refbias_binsize = 5
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
    OutputDiffRefBias(diffs_from_ref, reflens, fname, metric="median")
    OutputDiffRefBias(diffs_from_ref, reflens, fname, metric="invalid")

# Just confirm that the method doesn't throw an error
def test_OutputSampleCallrate(tmpdir):
    sample_calls = np.array([120, 10])
    samples = ['s1', 's2']
    fname = str(tmpdir / "test_qc1.pdf")
    OutputSampleCallrate(sample_calls, samples, fname)

# Just confirm that the method doesn't throw an error
def test_OutputChromCallrate(tmpdir):
    chrom_calls = {'chr1': 100, 'chr2': 200}
    fname = str(tmpdir / "test_qc2.pdf")
    OutputChromCallrate(chrom_calls, fname)

def test_output_location_is_directory(tmpdir, vcfdir, capsys):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")
    args.out = str(tmpdir)+os.path.sep
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

def test_refbias_options(tmpdir, vcfdir, capsys):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    vcf = os.path.join(qcdir, "test_popstr.vcf")

    # Test metrics
    args = base_argparse(tmpdir)
    args.vcf = vcf
    args.refbias_metric = "median"
    retcode = main(args)
    assert retcode == 0
    args.refbias_metric = "mean"
    retcode = main(args)

    # Test mingts
    args = base_argparse(tmpdir)
    args.vcf = vcf
    args.refbias_mingts = 1
    retcode = main(args)
    assert retcode == 0
    args.refbias_mingts = -1
    retcode = main(args)
    assert retcode == 1
    assert "refbias-mingts must be" in capsys.readouterr().err

    
    # Test binsize
    args = base_argparse(tmpdir)
    args.vcf = vcf
    args.refbias_binsize = 5
    retcode = main(args)
    assert retcode == 0
    args.refbias_binsize = -1
    retcode = main(args)
    assert retcode == 1
    assert "refbias-binsize must be" in capsys.readouterr().err

    # Test xrange
    args = base_argparse(tmpdir)
    args.vcf = vcf
    args.refbias_xrange_min = 0
    args.refbias_xrange_max = 80
    retcode = main(args)
    assert retcode == 0
    args.refbias_xrange_min = 100
    retcode = main(args)
    assert retcode == 1
    assert "refbias-xrange" in capsys.readouterr().err

def test_asking_for_qual_plot_fails_when_no_qual_info(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "test_ExpansionHunter_reheadered.vcf")
    args.quality = ['per-locus']
    retcode = main(args)
    assert "doesn't have quality scores" in capsys.readouterr().err
    assert retcode == 1

def test_dont_make_qual_plot_when_no_qual_info(tmpdir, vcfdir, capsys):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "test_ExpansionHunter_reheadered.vcf")
    retcode = main(args)
    assert len(glob.glob(args.out + "*quality*")) == 0
    assert retcode == 0
    assert "skipping all quality plots" in capsys.readouterr().out

def test_make_default_qual_plot_few_samples(tmpdir, vcfdir):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    retcode = main(args)
    assert retcode == 0
    assert os.path.exists(args.out + "-quality.pdf")
    assert len(glob.glob(args.out + "-quality-*")) == 0

def test_make_default_qual_plot_many_samples(tmpdir, vcfdir):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "many_samples.vcf.gz")
    retcode = main(args)
    assert retcode == 0
    assert os.path.exists(args.out + "-quality.pdf")
    assert len(glob.glob(args.out + "-quality-*")) == 0


# From https://stackoverflow.com/a/1073382/2966505
def _clear_folder(folder):
    for root, dirs, files in os.walk(folder):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))


def test_make_single_qual_plots_explicit(tmpdir, vcfdir):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    for qual in _QualityTypes.__members__.values():
        qual = qual.value
        args.quality = [qual]
        _clear_folder(str(tmpdir))
        retcode = main(args)
        assert retcode == 0
        assert os.path.exists(args.out + "-quality-" + qual + ".pdf")
        assert len(glob.glob(args.out + "-quality*")) == 1


def test_make_all_qual_plots(tmpdir, vcfdir):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    for qual in _QualityTypes.__members__.values():
        qual = qual.value
        args.quality.append(qual)
    retcode = main(args)
    assert retcode == 0
    assert not os.path.exists(args.out + "-quality.pdf")
    for qual in _QualityTypes.__members__.values():
        qual = qual.value
        assert os.path.exists("{}-quality-{}.pdf".format(args.out, qual))

def test_all_qual_plots_with_ignore_no_call(tmpdir, vcfdir):
    # check retcode is zero and appropriate files are/are not
    # created, don't do any checks for content
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "few_loci.vcf")
    args.quality_ignore_no_call = True
    for qual in _QualityTypes.__members__.values():
        qual = qual.value
        args.quality.append(qual)
    retcode = main(args)
    assert retcode == 0
    assert not os.path.exists(args.out + "-quality.pdf")
    for qual in _QualityTypes.__members__.values():
        qual = qual.value
        assert os.path.exists("{}-quality-{}.pdf".format(args.out, qual))


def test_output_all_files(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "many_samples_multiple_chroms.vcf.gz")
    retcode = main(args)
    assert retcode == 0
    stdout = capsys.readouterr().out
    for suffix in ["-sample-callnum",
                   "-chrom-callnum",
                   "-diffref-histogram",
                   "-diffref-bias",
                   "-quality"]:
        outfile = str(tmpdir / "test_qc") + suffix + ".pdf"
        assert "Producing " + outfile in stdout
        assert os.path.exists(outfile)


def test_omit_callnum_one_chrom(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "many_samples.vcf.gz")
    retcode = main(args)
    assert retcode == 0
    stdout = capsys.readouterr().out
    for suffix in ["-sample-callnum",
                   "-diffref-histogram",
                   "-diffref-bias",
                   "-quality"]:
        outfile = str(tmpdir / "test_qc") + suffix + ".pdf"
        assert "Producing " + outfile in stdout
        assert os.path.exists(outfile)
    outfile = str(tmpdir / "test_qc") + "-chrom-callnum.pdf"
    assert not os.path.exists(outfile)
    assert "skipping " + outfile in stdout


def test_omit_callnum_one_sample(tmpdir, vcfdir, capsys):
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(vcfdir, "one_sample_multiple_chroms.vcf.gz")
    args.refbias_mingts = 1
    retcode = main(args)
    assert retcode == 0
    stdout = capsys.readouterr().out
    for suffix in ["-chrom-callnum",
                   "-diffref-histogram",
                   "-diffref-bias",
                   "-quality"]:
        outfile = str(tmpdir / "test_qc") + suffix + ".pdf"
        assert "Producing " + outfile in stdout
        assert os.path.exists(outfile)
    outfile = str(tmpdir / "test_qc") + "-sample-callnum.pdf"
    assert not os.path.exists(outfile)
    assert "skipping " + outfile in stdout


def test_main(tmpdir, vcfdir):
    qcdir = os.path.join(vcfdir, "qc_vcfs")
    # correct vcf
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")
    retcode = main(args)
    assert retcode == 0

    # vcf file with no contig
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_popstr_nocontig.vcf")
    retcode = main(args)
    assert retcode == 0
    args.vcf = os.path.join(qcdir, "test_popstr.vcf")

    # Set sample list
    args.samples = os.path.join(qcdir, "test_samplelist.txt")
    retcode = main(args)
    assert retcode == 0

    args.samples = os.path.join(qcdir, "test_fakesample.txt")
    retcode = main(args)
    assert retcode == 0

    # Non existent vcf
    args = base_argparse(tmpdir)
    args.vcf = os.path.join(qcdir, "test_non_exist.vcf")
    retcode = main(args)
    assert retcode == 1
