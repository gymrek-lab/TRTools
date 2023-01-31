import argparse

import cyvcf2
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm

from .. import associaTR

# necessary to make comparisons actually equal
associaTR.load_and_filter_genotypes.allele_len_precision = 10
associaTR.pval_precision = 5

@pytest.fixture
def test_associaTR_dir(vcfdir):
    return vcfdir + "/associaTR"

@pytest.fixture(name='args')
def args(test_associaTR_dir):
    args = argparse.Namespace()
    args.outfile = 'test_association_results.tsv'
    # this has dosages in it, but will be ignored unless dosages are specified
    args.tr_vcf = test_associaTR_dir + "/many_samples_biallelic_dosages.vcf.gz"
    args.phenotype_name = 'test_pheno'
    args.traits = [test_associaTR_dir + "/traits_0.npy"]
    args.vcftype = 'auto'
    args.same_samples = False
    args.sample_list = None
    args.region = None
    args.non_major_cutoff = 0
    args.beagle_dosages = False
    args.plotting_phenotype = None
    args.paired_genotype_plot = False
    args.plot_phenotype_residuals = False
    args.plotting_ci_alphas = []
    args.imputed_ukb_strs_paper_period_check = False
    return args

format_precision = 2
def my_format(float_):
    return np.format_float_scientific(float_, precision=format_precision, unique=False)

diff_size = 2 # allow for a diff of ~2% in percision slop
def comp_floats(float1, float2):
    if not np.sign(float1) == np.sign(float2):
        assert False, (float1, float2)
    f1 = my_format(abs(float1))
    f2 = my_format(abs(float2))
    if not f1[:2] == f2[:2]:
        assert False, (float1, float2)
    if not abs(int(f1[2:4]) - int(f2[2:4])) <= diff_size:
        assert False, (float1, float2)
    if not f1[5:] == f2[5:]:
        assert False, (float1, float2)

def compare_my_gwas_to_plink(my_gwas_file, plink_file, phenotype_name, skip_filtered=False):
    out_df = pd.read_csv(my_gwas_file, sep='\t')
    plink_df = pd.read_csv(plink_file, sep='\t')

    if skip_filtered:
        out_df = out_df.loc[out_df['locus_filtered'] == 'False', :].reset_index()

        plink_df = plink_df.loc[plink_df['ERRCODE'] == '.', :]
        # the above would strip our loci with only one allele, but plink
        # dos not error those, so manually do that here
        plink_df = plink_df.loc[plink_df['REF'].str.len() != plink_df['ALT'].str.len(), :].reset_index()

    assert out_df.shape[0] == plink_df.shape[0]
    for line in range(out_df.shape[0]):
        out_p = out_df.loc[line, 'p_' + phenotype_name]
        # This was filtered, plink does something different with filtered loci, that's okayh
        if not skip_filtered and np.isnan(out_p):
            if ',' in out_df.loc[line, 'alleles']:
                assert plink_df.loc[line, 'ERRCODE'] != '.'
            # otherwise the only difference is sequence, not length, so plink
            # will test but my code won't, which is okay
            continue
        comp_floats(out_p, plink_df.loc[line, 'P'])
        ref_len = out_df.loc[line, 'ref_len']
        alleles = [float(x) for x in out_df.loc[line, 'alleles'].split(',')]
        assert len(alleles) == 2
        copy_count_diff = abs(alleles[0] - alleles[1])
        sign = 1 if ref_len == min(alleles) else -1
        comp_floats(out_df.loc[line, 'coeff_' + phenotype_name]*copy_count_diff * sign, plink_df.loc[line, 'BETA'])
        comp_floats(out_df.loc[line, 'se_' + phenotype_name]*copy_count_diff, plink_df.loc[line, 'SE'])

def test_one_trait_file(args, test_associaTR_dir):
    args.same_samples = True
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single.plink2.trait_0.glm.linear", args.phenotype_name)

def test_two_trait_files(args, test_associaTR_dir):
    args.same_samples = True
    args.traits = [test_associaTR_dir + "/traits_0.npy", test_associaTR_dir + "/traits_1.npy"]
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/combined.plink2.trait_0.glm.linear", args.phenotype_name)

def test_one_trait_file_sample_merge(args, test_associaTR_dir):
    args.traits = [test_associaTR_dir + "/traits_0_40_samples.npy"]
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_40.plink2.trait_0.glm.linear", args.phenotype_name)

def test_two_trait_files_sample_merge(args, test_associaTR_dir):
    args.traits = [test_associaTR_dir + "/traits_0_40_samples.npy", test_associaTR_dir + "/traits_1_45_samples.npy"]
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/combined_35.plink2.trait_0.glm.linear", args.phenotype_name)

def test_one_trait_file_sample_subset(args, test_associaTR_dir):
    args.same_samples = True
    args.sample_list = test_associaTR_dir + "/samples_6_to_45.txt"
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_40.plink2.trait_0.glm.linear", args.phenotype_name)

def test_one_trait_file_sample_merge_and_subset(args, test_associaTR_dir):
    args.traits = [test_associaTR_dir + "/traits_0_40_samples.npy"]
    args.sample_list = test_associaTR_dir + "/45_samples.txt"
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_35.plink2.trait_0.glm.linear", args.phenotype_name)

def test_region(args, test_associaTR_dir):
    args.same_samples = True
    associaTR.main(args)
    with open(args.outfile) as outfile:
        lines = outfile.readlines()
    args.region = '1:993134-3781638'
    associaTR.main(args)
    with open(args.outfile) as outfile:
        assert next(outfile) == lines[0]
        for idx, line in enumerate(outfile.readlines()):
            assert line == lines[idx + 77]
        assert idx ==  366 - 77  # happens to be the number of variants in that range

    args.region = '2:993134-3781638'
    associaTR.main(args)
    with open(args.outfile) as outfile:
        # no output
        assert len(outfile.readlines()) == 1

def test_non_major_count_cutoff(args, test_associaTR_dir):
    args.same_samples = True
    args.non_major_cutoff = 5
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_cutoff_5.plink2.trait_0.glm.linear", args.phenotype_name, skip_filtered=True)

def test_dosages(args, test_associaTR_dir):
    args.same_samples = True
    args.beagle_dosages = True
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_dosages.plink2.trait_0.glm.linear", args.phenotype_name)

def test_dosage_sample_subset(args, test_associaTR_dir):
    args.same_samples = True
    args.beagle_dosages = True
    args.sample_list = test_associaTR_dir + "/samples_6_to_45.txt"
    associaTR.main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single_40_dosages.plink2.trait_0.glm.linear", args.phenotype_name)

# test multiallelic
# first multiallelic allele has 3 separate alleles
# second has 3 alleles, but lens of 0 and 2 are the same
# first multiallelic allele has allelec counts
#  95, 3, 2
# and dosage allele totals
# [62.58339285850525, 19.654988, 17.76162], non major = 37.416608
# second has allele counts (only 49 called samples)
#  89, 9
# and dosage allele totals
# [80.34407, 17.65593]

# specifically tests if recoding and coalescing alleles is working properly
def test_multiallelic(args, test_associaTR_dir):
    args.same_samples = True
    args.tr_vcf = test_associaTR_dir + "/many_samples_multiallelic_dosages.vcf.gz"
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    covars = np.load(args.traits[0])
    covars = np.hstack((covars, np.ones((covars.shape[0], 1))))
    outcome = covars[:, 0].copy()
    vcf = cyvcf2.VCF(args.tr_vcf)

    # test var 1
    var = next(vcf)
    gts = var.genotype.array()[:, :-1]
    new_gts = np.full(gts.shape, np.nan)
    # recode based on lengths compared to ref
    new_gts[gts == 0] = 0
    new_gts[gts == 1] = -1
    new_gts[gts == 2] = 1
    summed_gts = np.sum(new_gts, axis=1)
    covars[:, 0] = summed_gts
    result = sm.OLS(outcome, covars).fit()
    comp_floats(out_df.loc[0, 'p_' + args.phenotype_name], result.pvalues[0])
    comp_floats(out_df.loc[0, 'coeff_' + args.phenotype_name], result.params[0])
    comp_floats(out_df.loc[0, 'se_' + args.phenotype_name], result.bse[0])
   
    # test var 2
    # remove the missing sample
    covars = covars[1:, :]
    outcome = outcome[1:]
    var = next(vcf)
    gts = var.genotype.array()[1:, :-1]
    new_gts = np.full(gts.shape, np.nan)
    # recode based on lengths compared to ref
    new_gts[gts == 0] = 0
    new_gts[gts == 1] = -2
    new_gts[gts == 2] = 0
    summed_gts = np.sum(new_gts, axis=1)
    covars[:, 0] = summed_gts
    result = sm.OLS(outcome, covars).fit()
    comp_floats(out_df.loc[1, 'p_' + args.phenotype_name], result.pvalues[0])
    comp_floats(out_df.loc[1, 'coeff_' + args.phenotype_name], result.params[0])
    comp_floats(out_df.loc[1, 'se_' + args.phenotype_name], result.bse[0])

# specifically tests if recoding and coalescing alleles is working properly
def test_multiallelic_dosages(args, test_associaTR_dir):
    args.same_samples = True
    args.beagle_dosages = True
    args.tr_vcf = test_associaTR_dir + "/many_samples_multiallelic_dosages.vcf.gz"
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    covars = np.load(args.traits[0])
    covars = np.hstack((covars, np.ones((covars.shape[0], 1))))
    outcome = covars[:, 0].copy()
    vcf = cyvcf2.VCF(args.tr_vcf)

    # test var 1
    var = next(vcf)
    summed_dosages = np.zeros(len(vcf.samples))
    # recode based on lengths compared to ref
    summed_dosages -= (var.format('AP1') + var.format('AP2'))[:, 0]
    summed_dosages += (var.format('AP1') + var.format('AP2'))[:, 1]
    covars[:, 0] = summed_dosages
    result = sm.OLS(outcome, covars).fit()
    comp_floats(out_df.loc[0, 'p_' + args.phenotype_name], result.pvalues[0])
    comp_floats(out_df.loc[0, 'coeff_' + args.phenotype_name], result.params[0])
    comp_floats(out_df.loc[0, 'se_' + args.phenotype_name], result.bse[0])
   
    # test var 2
    # remove the missing sample
    covars = covars[1:, :]
    outcome = outcome[1:]
    var = next(vcf)
    # both other alleles are length zero relative to ref
    summed_dosages = -2*(var.format('AP1') + var.format('AP2'))[1:, 0]
    covars[:, 0] = summed_dosages
    result = sm.OLS(outcome, covars).fit()
    comp_floats(out_df.loc[1, 'p_' + args.phenotype_name], result.pvalues[0])
    comp_floats(out_df.loc[1, 'coeff_' + args.phenotype_name], result.params[0])
    comp_floats(out_df.loc[1, 'se_' + args.phenotype_name], result.bse[0])

def test_multiallelic_cutoff(args, test_associaTR_dir):
    args.same_samples = True
    args.tr_vcf = test_associaTR_dir + "/many_samples_multiallelic_dosages.vcf.gz"
    args.non_major_cutoff = 3
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert np.all(~out_df.loc[:, 'coeff_' + args.phenotype_name].isna())
    args.non_major_cutoff = 8
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert np.isnan(out_df.loc[0, 'coeff_' + args.phenotype_name])
    assert ~np.isnan(out_df.loc[1, 'coeff_' + args.phenotype_name])
    args.non_major_cutoff = 10
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert np.all(out_df.loc[:, 'coeff_' + args.phenotype_name].isna())

def test_dosage_multiallelic_cutoff(args, test_associaTR_dir):
    args.same_samples = True
    args.beagle_dosages = True
    args.tr_vcf = test_associaTR_dir + "/many_samples_multiallelic_dosages.vcf.gz"
    args.non_major_cutoff = 10
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert np.all(~out_df.loc[:, 'coeff_' + args.phenotype_name].isna())
    args.non_major_cutoff = 20
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert ~np.isnan(out_df.loc[0, 'coeff_' + args.phenotype_name])
    assert np.isnan(out_df.loc[1, 'coeff_' + args.phenotype_name])
    args.non_major_cutoff = 38
    associaTR.main(args)
    out_df = pd.read_csv(args.outfile, sep='\t')
    assert np.all(out_df.loc[:, 'coeff_' + args.phenotype_name].isna())

# TODO test fields other than p coeff and se

# TODO test binary

# TODO test plotting phenotype, in addition to paired genotype plot and residuals and cis
