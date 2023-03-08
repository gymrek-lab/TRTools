import argparse
import os

import pytest

from ..statSTR import *


# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcfs = None
    args.out = str(tmpdir /  "test")
    args.vcftype = "auto"
    args.samples = None
    args.sample_prefixes = None
    args.plot_afreq = False
    args.region = None
    args.thresh = False
    args.afreq = False
    args.acount = False
    args.nalleles = False
    args.nalleles_thresh = 0.01
    args.hwep = False
    args.het = False
    args.use_length = False
    args.mean = False
    args.mode = False
    args.var = False
    args.numcalled = False
    args.entropy = False
    args.precision = 4
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
def test_Stats(args, vcfdir, capsys):
    fname = os.path.join(vcfdir, "few_samples_few_loci.vcf.gz")
    args.vcf = fname
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
    args.nalleles = True
    assert main(args) == 0
    args.uselength = True
    assert main(args) == 0
    args.region = "chr1:3045469-3045470"
    assert main(args) == 0
    args.samples = os.path.join(vcfdir, "fewer_samples.txt")
    assert main(args) == 0
    # test that using only nonexistant samples errors out
    args.samples = os.path.join(vcfdir, "missing_samples.txt")
    assert main(args) == 1
    assert 'no samples' in capsys.readouterr().err.lower()

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

# True if equal
def _imprecise_compare(token1, token2):
    token1 = float(token1)
    token2 = float(token2)
    if np.isnan(token1) and np.isnan(token2):
        return True
    elif token1 == 0 and token2 == 0:
        return True
    elif token1 < 1e-3:
        token1 = _precision(token1, 2)
        token2 = _precision(token2, 2)
    elif token1 < 1e-4:
        # only use one bit of precision
        # so that calculations with numeric instability can
        # still be compared at very low values
        token1 = _precision(token1, 1)
        token2 = _precision(token2, 1)
    else:
        token1 = _precision(token1, 3)
        token2 = _precision(token2, 3)

    return abs(token1 - token2) < token1*1e-6


# Converts 5.461e-6 to 5e-6
def _precision(val, prec):
    exp = 10**(np.floor(np.log10(val)) - prec + 1)
    return np.floor(val/exp)*exp


def _make_allele_dict(allele_dict_string):
    d = {}
    for pair in allele_dict_string.split(','):
        key, val = pair.split(':')
        d[key] = val
    return d


def _comp_output_files(fname1, fname2):
    with open(fname1) as file1, open(fname2) as file2:
        i1 = iter(file1)
        i2 = iter(file2)
        header1 = next(i1)
        header2 = next(i2)
        assert header1 == header2
        header = header1.split('\t')

        linenum = 1
        while True:
            try:
                l1 = next(i1)
            except StopIteration:
                try:
                    l2 = next(i2)
                    raise ValueError(
                        "file1 is a truncated version of file2, equal up till truncation"
                    )
                except StopIteration:
                    # files are the same and ended at the same time
                    return
            try:
                l2 = next(i2)
            except StopIteration:
                raise ValueError(
                    "file2 is a truncated version of file1, equal up till truncation"
                )
            linenum += 1
            l1 = l1.split('\t')
            l2 = l2.split('\t')
            if len(l1) != len(l2):
                raise ValueError(
                    ("Locus {} (line num {}) is of different "
                     "length between the two files").format(
                         l1[0:3], linenum)
                )
            for idx, (token1, token2) in enumerate(zip(l1, l2)):
                column = header[idx]
                # check some columns for a less than exact match
                if ('hwep' in column or
                        'het' in column or
                        'entropy' in column or
                        'var' in column or
                        'mean' in column or
                        'mode' in column or
                        'thresh' in column):
                    if not _imprecise_compare(token1, token2):
                        raise ValueError(
                            ("Locus {} (line num {}) is different "
                             "between the two files at column {} "
                             "with vals file1:{} file2:{}").format(
                                 l1[0:3], linenum, header[idx], token1, token2)
                        )
                elif 'afreq' in column or 'acount' in column:
                    if token1 == '.' or token2 == '.':
                        if token1 == '.' and token2 == '.':
                            continue
                        else:
                            raise ValueError(
                                ("Locus {} (line num {}) is different "
                                 "between the two files at column {}"
                                 "file1:{} file2:{}").format(
                                     l1[0:3], linenum, header[idx], token1,
                                     token2)
                            )
                    dict1 = _make_allele_dict(token1)
                    dict2 = _make_allele_dict(token2)
                    if len(dict1) != len(dict2):
                        raise ValueError(
                            ("Locus {} (line num {}) is different "
                             "between the two files at column {}"
                             "with different numbers of alleles! file1:{} "
                             "file2:{}").format(
                                 l1[0:3], linenum, header[idx], token1, token2)
                        )
                    for key in dict1:
                        if not key in dict2:
                            raise ValueError(
                                ("Locus {} (line num {}) is different "
                                 "between the two files at column {}"
                                 "allele {} is in file1 but not "
                                 "file2").format(
                                     l1[0:3], linenum, header[idx], key)
                            )
                        if not _imprecise_compare(dict1[key],
                                                  dict2[key]):
                            raise ValueError(
                                ("Locus {} (line num {}) is different "
                                 "between the two files at column {} "
                                 "with vals for allele {} file1:{} "
                                 "file2:{}").format(
                                     l1[0:3], linenum, header[idx], key,
                                     token1, token2)
                            )
                elif token1 != token2:
                    raise ValueError(
                        ("Locus {} (line num {}) is different "
                         "between the two files at column {} "
                         "with vals file1:{} file2:{}").format(
                             l1[0:3], linenum, header[idx], token1, token2
                         )
                    )


def test_output(args, vcfdir, statsdir):
    """
    Run statSTR on a file which statSTR has been run on in the past
    and confirm the results haven't changed

    If you need to change the output file, generate it with precision 4
    """
    fname = os.path.join(vcfdir, "many_samples.vcf.gz")
    args.vcf = fname
    args.thresh = True
    args.afreq = True
    args.acount = True
    args.nalleles = True
    args.nalleles_thresh = 0.1
    args.hwep = True
    args.het = True
    args.entropy = True
    args.mean = True
    args.mode = True
    args.var = True
    args.numcalled = True
    # exclude an allele which doesn't have
    # reproducible stats even up to two decimal places
    assert main(args) == 0
    _comp_output_files(
        args.out + ".tab",
        os.path.join(statsdir, "many_samples_all.tab")
    )


def test_output_samplestrat(args, vcfdir, statsdir):
    """
    Run statSTR on a file which statSTR has been run on in the past
    and confirm the results haven't changed

    If you need to change the output file, generate it with precision 4
    """
    fname = os.path.join(vcfdir, "many_samples.vcf.gz")
    args.vcf = fname
    args.samples = (
        os.path.join(vcfdir, "many_samples_subsample1.txt") + "," +
        os.path.join(vcfdir, "many_samples_subsample2.txt")
    )
    args.thresh = True
    args.afreq = True
    args.acount = True
    args.nalleles = True
    args.nalleles_thresh = 0.1
    args.hwep = True
    args.het = True
    args.entropy = True
    args.mean = True
    args.mode = True
    args.var = True
    args.numcalled = True
    # exclude an allele which doesn't have
    # reproducible stats even up to two decimal places
    assert main(args) == 0
    _comp_output_files(
        args.out + ".tab",
        os.path.join(statsdir, "many_samples_all_strat.tab")
    )


