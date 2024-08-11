import argparse
import gzip
import os

import pytest

from ..annotaTR import *

# Set up base argparser
@pytest.fixture(name='args')
def args_fixture(tmpdir):
    return args(tmpdir)

def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.vcftype = "auto"
    args.out = str(tmpdir / "test")
    args.outtype = ["vcf"]
    args.dosages = None
    args.ref_panel = None
    return args

@pytest.fixture
def testBeagleDir(vcfdir):
    return vcfdir + "/beagle"

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1
    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.ref_panel = fname
    retcode = main(args)
    assert retcode==1
    args.ref_panel = None
    args.vcf = os.path.join(vcfdir, "missing_samples.txt")
    retcode = main(args)
    assert retcode==1

def test_WrongOutDir(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.out = os.path.join("fakedir","test")
    retcode = main(args)
    assert retcode==1
    args.out = vcfdir + "/"
    retcode = main(args)
    assert retcode==1
   
def test_OutTypes(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.outtype = ["vcf", "pgen"]
    args.dosages = "bestguess_norm"
    retcode = main(args)
    assert retcode==0
    args.outtype = ["pgen", "vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["pgen"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["vcf", "vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["dummy"]
    retcode = main(args)
    assert retcode==1

def test_VCFType(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "badtype"
    retcode = main(args)
    assert retcode==1
    args.vcftype = "auto"
    retcode = main(args)
    assert retcode==0

def test_DosageType(args, vcfdir):
    # Non-beagle VCF
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.dosages = "bestguess"
    retcode = main(args)
    assert retcode==0
    args.dosages = "badtype"
    retcode = main(args)
    assert retcode==1
    args.dosages = "beagleap"
    retcode = main(args)
    assert retcode==1
    args.dosages = "beagleap_norm"
    retcode = main(args)
    assert retcode==1
    # Beagle VCF
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "hipstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.dosages = "beagleap_norm"
    retcode = main(args)
    assert retcode==0
    # No dosage specified should cause error for pgen output
    args.dosages = None
    args.outtype = ["pgen"]
    retcode = main(args)
    assert retcode==1
    # Non-norm dosage should also cause error for pgen output
    args.dosages = "beagleap"
    retcode = main(args)
    assert retcode==1
    # norm dosages should be fine for pgen
    args.dosages = "beagleap_norm"
    retcode = main(args)
    assert retcode==0
    
def test_LoadRefpanel(args, vcfdir):
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    retcode = main(args)
    assert retcode == 0
    args.vcftype = "auto"
    retcode = main(args)
    assert retcode == 0
    # Bad refpanel
    args.ref_panel = os.path.join(vcfdir, "missing_samples.txt")
    retcode = main(args)
    assert retcode == 1
    # Don't support refpanel if popstr
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.vcftype = "popstr"
    retcode = main(args)
    assert retcode == 1
    # Error if no TRs in refpanel
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel_nostrs.vcf.gz")
    args.vcftype = "hipstr"
    retcode = main(args)
    assert retcode == 1
    # Error if duplicate TR in ref panel
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel_duplocus.vcf.gz")
    args.vcftype = "hipstr"
    retcode = main(args)
    assert retcode == 1
    # Fail if missing a required info header
    args.vcf = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel_missinginfoheader.vcf.gz")
    args.vcftype = "hipstr"
    retcode = main(args)
    assert retcode == 1
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.vcftype = "eh"
    retcode = main(args)
    assert retcode == 1
    # Fail if wrong ref panel
    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_hipstr.sorted.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.vcftype = "hipstr"
    retcode = main(args)
    assert retcode == 1
   