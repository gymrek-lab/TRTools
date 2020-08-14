# pylint: disable=W0621,C0116,C0114
import argparse
import filecmp
import glob
import os

import pytest

from ..statSTR import *
from ..statSTR import _dist_plot_args


# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.out = str(tmpdir /  "test")
    args.vcftype = "auto"
    args.samples = None
    args.sample_prefixes = None
    args.plot_afreq = False
    args.region = None
    args.thresh = False
    args.afreq = False
    args.acount = False
    args.hwep = False
    args.het = False
    args.use_length = False
    args.mean = False
    args.mode = False
    args.var = False
    args.numcalled = False
    args.entropy = False
    args.precision = 2
    args.plot_dists = None
    args.tabfiles = None
    return args

@pytest.fixture
def args_allstats(args):
    args.thresh = True
    args.het = True
    args.entropy = True
    args.hwep = True
    args.numcalled = True
    args.mode = True
    args.mean = True
    args.var = True
    args.acount = True
    args.afreq = True
    return args

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

# Test the right file or directory
def test_RightFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

# Test all the statSTR options
def test_Stats(args_allstats, vcfdir):
    args = args_allstats
    fname = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.vcf = fname
    assert main(args) == 0
    args.uselength = True
    assert main(args) == 0
    args.region = "chr1:3045469-3045470"
    assert main(args) == 0
    args.samples = os.path.join(vcfdir, "fewer_samples.txt")
    assert main(args) == 0

def test_require_a_stat(args, vcfdir, capsys):
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    assert check_args(args) is None
    assert ("Please use at least one of the flags in the Stats group" in
            capsys.readouterr().err)

def test_require_tabix_for_regions(args, vcfdir, capsys):
    # args.samples = os.path.join(vcfdir, "test_gangstr_samples.txt")
    args.region = "chr1:3045469-3045470"
    args.vcf = os.path.join(vcfdir, "test_gangstr.vcf")
    args.thresh = True
    assert main(args) == 1
    assert 'bgzipped' in capsys.readouterr().err

def test_Stats2(args, vcfdir):
    args.vcf = os.path.join(vcfdir, "mergeSTR_vcfs", "test_file_gangstr1.vcf.gz")
    args.region = "chr1:3045469-3045470"
    assert main(args)==0

def test_PlotAfreq(args, vcfdir):
    fname = os.path.join(vcfdir, "test_gangstr.vcf")
    args.vcf = fname
    args.plot_afreq = True
    assert main(args)==0

def test_output(args_allstats, vcfdir, statsdir):
    """
    Run statSTR on a file which statSTR has been run on in the past
    and confirm the results haven't changed
    """
    args = args_allstats
    fname = os.path.join(vcfdir, "many_samples.vcf.gz")
    args.vcf = fname
    # exclude an allele which doesn't have
    # reproducible stats even up to two decimal places
    args.region = "1:1-3683401"
    assert main(args) == 0
    assert filecmp.cmp(
        args.out + ".tab",
        os.path.join(statsdir, "many_samples_all.tab")
    )
    assert len(glob.glob(args.out + "*.png")) == 0


def test_output_samplestrat(args_allstats, vcfdir, statsdir):
    """
    Run statSTR on a file which statSTR has been run on in the past
    and confirm the results haven't changed
    """
    args = args_allstats
    fname = os.path.join(vcfdir, "many_samples.vcf.gz")
    args.vcf = fname
    args.samples = (
        os.path.join(vcfdir, "many_samples_subsample1.txt") + "," +
        os.path.join(vcfdir, "many_samples_subsample2.txt")
    )
    # exclude an allele which doesn't have
    # reproducible stats even up to two decimal places
    args.region = "1:1-1833757"
    assert main(args) == 0
    assert filecmp.cmp(
        args.out + ".tab",
        os.path.join(statsdir, "many_samples_all_strat.tab")
    )


# Test plot_dists
# Tests can only ensure plotting methods don't crash.
# Run these locally to view plots to make sure they look sane

def test_nocrash_plot_dists_all(args_allstats, vcfdir):
    args = args_allstats
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_all_smooth(args_allstats, vcfdir):
    args = args_allstats
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.plot_dists = 'smooth'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_subset(args, vcfdir):
    args.thresh = True
    args.mode = True
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        fname = "{}-{}.png".format(args.out, stat)
        if stat in {'thresh', 'mode'}:
            assert os.path.exists(fname)
        else:
            assert not os.path.exists(fname)


def test_crash_plot_dists_no_plottable_stats(args, vcfdir, capsys):
    args.afreq = True
    args.plot_dists = 'histo'
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    assert check_args(args) is None
    assert "corresponding statistic to plot" in capsys.readouterr().err


def test_crash_plot_dists_bad_arg_combos(args, vcfdir, statsdir, capsys):
    args.thresh = True
    args.vcf = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.plot_dists = 'histo'
    args.tabfiles = os.path.join(statsdir, "many_samples_all.tab")
    assert check_args(args) is None
    assert "not both" in capsys.readouterr().err
    args.vcf = None

    args.plot_dists = None
    assert check_args(args) is None
    assert "both --tabfiles and --plot-dists" in capsys.readouterr().err
    args.plot_dists = 'histo'

    args.out = 'stdout'
    assert check_args(args) is None
    assert "Cannot plot distributions to stdout" in capsys.readouterr().err
    args.out = "foo"

    args.samples = os.path.join(vcfdir, "fewer_samples.txt")
    assert check_args(args) is None
    assert ("Cannot specify --tabfiles and --samples" in
        capsys.readouterr().err)
    args.samples = None

    args.region = "1:1-1000000"
    assert check_args(args) is None
    assert ("Cannot specify --tabfiles and --region" in
        capsys.readouterr().err)
    args.region = None

def test_crash_plot_dists_bad_tabfile(args_allstats, capsys):
    args = args_allstats
    args.tabfiles = "foo"
    args.plot_dists = 'histo'
    assert main(args) == 1
    assert "foo does not exist" in capsys.readouterr().err


def test_nocrash_plot_dists_all_tabfile(args_allstats, statsdir):
    args = args_allstats
    args.tabfiles = os.path.join(statsdir, "many_samples_all.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_all_many_tabfiles(args_allstats, statsdir):
    args = args_allstats
    args.tabfiles = os.path.join(statsdir, "many_samples_all.tab")
    args.tabfiles += "," + args.tabfiles
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_subset_tabfile(args, statsdir):
    args.thresh = True
    args.mode = True
    args.tabfiles = os.path.join(statsdir, "many_samples_all.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        fname = "{}-{}.png".format(args.out, stat)
        if stat in {'thresh', 'mode'}:
            assert os.path.exists(fname)
        else:
            assert not os.path.exists(fname)


def test_crash_plot_dists_stat_not_in_tabfile(args, statsdir, capsys):
    args.mode = True
    args.tabfiles = os.path.join(statsdir, "afreq_mean_hwep.tab")
    args.plot_dists = 'histo'
    assert main(args) == 1
    assert "Expected to find column mode in" in capsys.readouterr().err


def test_nocrash_plot_dists_subset_many_tabfiles(args, statsdir):
    args.thresh = True
    args.tabfiles = os.path.join(statsdir, "many_samples_all.tab")
    args.tabfiles += "," + args.tabfiles
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        fname = "{}-{}.png".format(args.out, stat)
        if stat in {'thresh'}:
            assert os.path.exists(fname)
        else:
            assert not os.path.exists(fname)


def test_nocrash_plot_dists_diff_tabfiles_same_subset(args, statsdir):
    args.mean = True
    args.tabfiles = os.path.join(statsdir, "afreq_mean_hwep.tab")
    args.tabfiles += "," + os.path.join(statsdir, "acount_mean_het.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        fname = "{}-{}.png".format(args.out, stat)
        if stat in {'mean'}:
            assert os.path.exists(fname)
        else:
            assert not os.path.exists(fname)


def test_nocrash_plot_dists_prefixes(args_allstats, statsdir):
    args = args_allstats
    args.tabfiles = os.path.join(statsdir, "many_samples_all_strat.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_prefixes_tabfile(args_allstats, statsdir):
    args = args_allstats
    args.tabfiles = os.path.join(statsdir, "many_samples_all_strat.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_nocrash_plot_dists_prefixes_many_tabfiles(args_allstats, statsdir):
    args = args_allstats
    args.tabfiles = os.path.join(statsdir, "many_samples_all_strat.tab")
    args.tabfiles += ',' +  args.tabfiles
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        assert os.path.exists("{}-{}.png".format(
            args.out, stat))


def test_crash_plot_dists_tabfiles_different_prefixes(args, statsdir, capsys):
    args.mean = True
    args.tabfiles = os.path.join(statsdir, "many_samples_all_strat.tab")
    args.tabfiles += "," + os.path.join(statsdir,
                                        "many_samples_all_strat_diffnames.tab")
    args.plot_dists = 'histo'
    assert main(args) == 1
    assert "mean-2" in capsys.readouterr().err


def test_nocrash_plot_dists_diff_tabfiles_same_subset_same_prefixes(args,
                                                                    statsdir):
    args.mean = True
    args.tabfiles = os.path.join(statsdir, "afreq_mean_hwep.tab")
    args.tabfiles += ',' + os.path.join(statsdir, "acount_mean_het.tab")
    args.plot_dists = 'histo'
    assert main(args) == 0
    for stat in _dist_plot_args:
        fname = "{}-{}.png".format(args.out, stat)
        if stat == 'mean':
            assert os.path.exists(fname)
        else:
            assert not os.path.exists(fname)


def test_crash_plot_dists_prefixes_stat_not_in_tabfile(args, statsdir, capsys):
    args.thresh = True
    args.tabfiles = os.path.join(statsdir, "afreq_mean_hwep.tab")
    args.tabfiles += ',' + os.path.join(statsdir, "acount_mean_het.tab")
    args.plot_dists = 'histo'
    assert main(args) == 1
    assert 'thresh' in capsys.readouterr().err


def test_crash_plot_dists_prefixes_stat_not_in_tabfiles(args, statsdir, capsys):
    args.hwep = True
    args.tabfiles = os.path.join(statsdir, "afreq_mean_hwep.tab")
    args.tabfiles += ',' + os.path.join(statsdir, "acount_mean_het.tab")
    args.plot_dists = 'histo'
    assert main(args) == 1
    assert 'hwep' in capsys.readouterr().err

