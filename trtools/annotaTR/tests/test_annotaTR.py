import argparse
import cyvcf2
import gzip
import os

import pytest

from ..annotaTR import *
from trtools.testsupport.utils import assert_same_vcf

# Set up base argparser
@pytest.fixture(name='args')
def args_fixture(tmpdir):
    return args(tmpdir)

def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.vcftype = "auto"
    args.region = None
    args.out = str(tmpdir / "test")
    args.outtype = ["vcf"]
    args.dosages = None
    args.ref_panel = None
    args.match_refpanel_on = "rawalleles"
    args.ignore_duplicates = False
    args.debug = False
    args.chunk_size = 1000
    args.warn_on_AP_error = False
    return args

@pytest.fixture
def testBeagleDir(vcfdir):
    return vcfdir + "/beagle"

@pytest.fixture
def antrvcfdir(vcfdir):
   return os.path.join(vcfdir, "annotaTR_vcfs")

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
    args.outtype = ["vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["vcf", "vcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["dummy"]
    retcode = main(args)
    assert retcode==1
    args.outtype = ["pgen"]
    retcode = main(args)
    assert retcode==0
    # Check the pvar file can be ready by cyvcf2?
    #pvarfile = args.out + ".pvar"
    #test_reader = cyvcf2.VCF(pvarfile)

def test_NoOperation(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.outtype = ["vcf"]
    retcode = main(args)
    assert retcode==1

def test_VCFType(args, vcfdir):
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "badtype"
    args.dosages = "bestguess"
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
    fname = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_second_50_STRs_imputed.vcf.gz")
    args.vcf = fname
    args.vcftype = "hipstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_first_50_annotated.vcf.gz")
    args.warn_on_AP_error = True # should't fail with missing AP
    args.outtype = ["pgen","vcf"]
    retcode = main(args)
    assert retcode==0
    args.warn_on_AP_error = False # set back
    with pytest.raises(ValueError):
        main(args)
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
    args.match_refpanel_on = "trimmedalleles"
    retcode = main(args)
    assert retcode == 0
    args.match_refpanel_on = "locid"
    with pytest.raises(ValueError):
        main(args)
    args.match_refpanel_on = "rawalleles" # set back for future tests
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
    with pytest.raises(ValueError):
        main(args)
    args.ignore_duplicates = True
    retcode = main(args)
    assert retcode == 0
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

def test_TrimAlleles():
    ref_allele = "AAAT"
    alt_alleles = ["AAATA","AAATAA","AAATAAA"]
    new_ref, new_alt = TrimAlleles(ref_allele, alt_alleles)
    assert(new_ref == ".")
    assert(new_alt == ["A","AA","AAA"])
    ref_allele = "AAATGAC"
    alt_alleles = ["AAATAGAC","AAATAAGAC","AAATAAAGAC"]
    new_ref, new_alt = TrimAlleles(ref_allele, alt_alleles)
    assert(new_ref == ".")
    assert(new_alt == ["A","AA","AAA"])
    ref_allele = "GAC"
    alt_alleles = ["AGAC","AAGAC","AAAGAC"]
    new_ref, new_alt = TrimAlleles(ref_allele, alt_alleles)
    assert(new_ref == ".")
    assert(new_alt == ["A","AA","AAA"])
    ref_allele = "AATAAGTAAATAAATAAATAA"
    alt_alleles = ["AATAAGTAAATAAATAA"]
    new_ref, new_alt = TrimAlleles(ref_allele, alt_alleles)
    assert(new_ref == "TAAA")
    assert(new_alt[0] == ".")

"""
These tests run annotaTR and compare its output
to output that has been generated by a previous version of
annotaTR and saved in the repo. The results are expected
to be identical.

These tests are too strict and will often break because
annotaTR output has been intentionally changed
However, the presence of these tests is important because
it should prevent any unexpected changes in output.
If you've reviewed the change in output and find it acceptable,
use trtools/testsupport/sample_vcfs/annotaTR_vcfs/create_test_files.sh
to regenerate the test files with the new version of mergeSTR.
"""

def test_OutputFilesSame(args, vcfdir, antrvcfdir):
    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.dosages = "bestguess"
    assert main(args) == 0
    assert_same_vcf(args.out + ".vcf", antrvcfdir + "/gangstr_bestguess.vcf", max_lines_to_compare=200)

    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.dosages = "bestguess_norm"
    assert main(args) == 0
    assert_same_vcf(args.out + ".vcf", antrvcfdir + "/gangstr_bestguess_norm.vcf", max_lines_to_compare=200)

    args.vcf = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_hipstr.sorted.vcf.gz")
    args.vcftype = "hipstr"
    args.dosages = "bestguess_norm"
    assert main(args) == 0
    assert_same_vcf(args.out + ".vcf", antrvcfdir + "/hipstr_bestguess_norm.vcf", max_lines_to_compare=200)

    args.vcf = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_second_50_STRs_imputed.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_first_50_annotated.vcf.gz")
    args.vcftype = "hipstr"
    args.dosages = "bestguess_norm"
    assert main(args) == 0
    assert_same_vcf(args.out + ".vcf", antrvcfdir + "/hipstr_beagle.vcf", max_lines_to_compare=200)

    args.vcf = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.vcftype = "hipstr"
    args.dosages = "beagleap"
    args.match_refpanel_on = "trimmedalleles"
    assert main(args) == 0
    assert_same_vcf(args.out + ".vcf", antrvcfdir + "/beagleap_trimmed.vcf", max_lines_to_compare=200)
