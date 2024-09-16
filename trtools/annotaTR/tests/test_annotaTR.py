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
    args.update_ref_alt = False
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
    args.outtype = ["gzvcf"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["gzvcf", "pgen"]
    retcode = main(args)
    assert retcode==0
    args.outtype = ["gzvcf", "vcf"]
    retcode = main(args)
    assert retcode==1
    # Should get pgen error if input VCF has fewer
    # TRs than ref panel
    args.vcf = os.path.join(vcfdir, "beagle", "beagle_imputed_noTRs.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_tinyrefpanel.vcf.gz")
    args.match_refpanel_on = "locid"
    args.dosages = "bestguess_norm"
    args.outtype = ["pgen"]
    args.out = "test"
    retcode = main(args)
    assert retcode==1
    for fname_ext in ("pgen", "pvar", "psam"):
        if os.path.exists("test."+fname_ext):
            os.remove("test."+fname_ext)
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

def test_UpdateRefAlt(args, vcfdir):
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "hipstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.dosages = "beagleap"
    args.update_ref_alt = True
    # Won't work if matching on anything besides locid
    args.match_refpanel_on = "rawalleles"
    retcode = main(args)
    assert retcode==1
    # Won't work with no ref panel
    args.match_refpanel_on = "locid"
    args.ref_panel = None
    retcode = main(args)
    assert retcode==1

    # Try on good file with alleles that do match refpanel
    args.match_refpanel_on = "locid"
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_goodalleles.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_tinyrefpanel.vcf.gz")
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

    # Try on dummy file with bad alleles that don't match refpanel
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_badalleles.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_tinyrefpanel.vcf.gz")
    args.vcf = fname
    with pytest.raises(ValueError):
        main(args)


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

def test_LoadRegion(args, vcfdir):
    # Test good region
    fname = os.path.join(vcfdir, "dumpSTR_vcfs", "trio_chr21_gangstr.sorted.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.dosages = "bestguess"
    args.region = "chr21:9489666-9546720"
    retcode = main(args)
    assert retcode==0
    
    # Test region with refpanel
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "gangstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    retcode = main(args)
    args.region = "chr21:14282813-14303433"
    retcode = main(args)
    assert retcode==0

    # Region not in ref panel
    # Fails because we find no TRs
    args.region = "chr19:14282813-14303433"
    retcode = main(args)
    assert retcode==1

    # Malformatted region
    # Fails because we find no TRs
    args.region = "XXXXX"
    retcode = main(args)
    assert retcode==1

    # Note, cyvcf2 doesn't complain about malformatted regions
    # and just will not return any intervals
    # TODO: We might want to check regions here and in other tools
    # like statstr and prancstr where users can set a region

def test_LoadRefpanel(args, vcfdir):
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "hipstr"
    args.ref_panel = os.path.join(vcfdir, "beagle", "beagle_refpanel.vcf.gz")
    args.match_refpanel_on = "rawalleles"    
    args.debug = True
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
    # Test on example where locid should work
    args.vcf = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_second_50_STRs_imputed.vcf.gz")
    args.ref_panel = os.path.join(vcfdir, "beagle", "1kg_snpstr_21_first_100k_first_50_annotated.vcf.gz")
    args.match_refpanel_on = "locid"
    args.vcftype = "hipstr"
    retcode = main(args)
    assert retcode == 0
    # Invalid match option
    args.match_refpanel_on = "badoption"
    with pytest.raises(ValueError):
        GetLocusKey(None, match_on="bad match") 
    retcode = main(args)
    assert retcode == 1
    args.match_refpanel_on = "rawalleles" # set back for future tests
    # Load mix of SNPs/TRs but no ref panel
    fname = os.path.join(vcfdir, "beagle", "beagle_imputed_withap.vcf.gz")
    args.vcf = fname
    args.vcftype = "hipstr"
    args.dosages = "bestguess"
    args.ref_panel = None
    args.debug = False
    retcode = main(args)
    assert retcode == 1
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

def test_CheckAlleleCompatibility():
    # Alleles identical
    panel_ref = "AAT"
    panel_alt = ["AATAAT","AATAATAAT"]
    record_ref = "AAT"
    record_alt = ["AATAAT","AATAATAAT"]
    assert CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    # Alleles in target subset of those in panel
    panel_ref = "AATAAT"
    panel_alt = ["AATAATAAT","AATAATAATAAT"]
    record_ref = "AAT"
    record_alt = ["AATAAT","AATAATAAT"]
    assert CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    panel_ref = "AATAAG"
    panel_alt = ["AATAATAAG","AATAATAATAAG"]
    record_ref = "AAG"
    record_alt = ["AATAAG","AATAATAAG"]
    assert CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    panel_ref = "AAGAAT"
    panel_alt = ["AAGAATAAT","AATAATAATAAT"]
    record_ref = "AAG"
    record_alt = ["AAGAAT","AATAATAAT"]
    assert CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    # Different numbers of alleles
    panel_ref = "AAGAAT"
    panel_alt = ["AAGAATAAT"]
    record_ref = "AAG"
    record_alt = ["AAGAAT","AATAATAAT"]
    assert not CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    panel_ref = "AAGAAT"
    panel_alt = ["AAGAATAAT","AATAATAATAAT"]
    record_ref = "AAG"
    record_alt = ["AAGAAT"]
    assert not CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    # Alleles in target not a subset of refpanel
    panel_ref = "AATAAT"
    panel_alt = ["AATAAGAAT","AATAATAATAAT"]
    record_ref = "AAT"
    record_alt = ["AATAAT","AATAATAAT"]
    assert not CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

    # Offsets not the same
    panel_ref = "AATAAT"
    panel_alt = ["AATAATAAT","AATAATAATAATAAT"]
    record_ref = "AAT"
    record_alt = ["AATAAT","AATAATAAT"]
    assert not CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt)

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