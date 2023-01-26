import argparse

import numpy as np
import pandas as pd
import pytest

from ..associaTR import *

# necessary to make comparisons actually equal
load_and_filter_genotypes.allele_len_precision = 10

@pytest.fixture
def test_associaTR_dir(vcfdir):
    return vcfdir + "/associaTR"

@pytest.fixture(name='args')
def args(test_associaTR_dir):
    args = argparse.Namespace()
    args.outfile = 'test_association_results.tsv'
    args.str_vcf = test_associaTR_dir + "/many_samples_biallelic.vcf.gz"
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

def my_format(float_):
    return np.format_float_scientific(float_, precision=pval_precision, unique=False)

def compare_my_gwas_to_plink(my_gwas_file, plink_file, phenotype_name):
    out_df = pd.read_csv(my_gwas_file, sep='\t')
    plink_df = pd.read_csv(plink_file, sep='\t')
    assert out_df.shape[0] == plink_df.shape[0]
    for line in range(out_df.shape[0]):
        #print(out_df.loc[line, :])
        out_p = out_df.loc[line, 'p_' + phenotype_name]
        # This was filtered, plink does something different with filtered loci, that's okayh
        if out_p == 1:
            continue
        assert my_format(out_p) == my_format(plink_df.loc[line, 'P'])
        ref_len = out_df.loc[line, 'ref_len']
        alleles = [float(x) for x in out_df.loc[line, 'alleles'].split(',')]
        assert len(alleles) == 2
        copy_count_diff = abs(alleles[0] - alleles[1])
        sign = 1 if ref_len == min(alleles) else -1
        assert my_format(out_df.loc[line, 'coeff_' + phenotype_name]*copy_count_diff * sign) == my_format(plink_df.loc[line, 'BETA'])
        assert my_format(out_df.loc[line, 'coeff_' + phenotype_name]*copy_count_diff * sign) == my_format(plink_df.loc[line, 'BETA'])
        assert my_format(out_df.loc[line, 'se_' + phenotype_name]*copy_count_diff) == my_format(plink_df.loc[line, 'SE'])

def test_one_trait_file(args, test_associaTR_dir):
    args.same_samples = True
    main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/single.plink2.trait_0.glm.linear", args.phenotype_name)
        
def test_two_trait_files(args, test_associaTR_dir):
    args.same_samples = True
    args.traits = [test_associaTR_dir + "/traits_0.npy", test_associaTR_dir + "/traits_1.npy"]
    main(args)
    compare_my_gwas_to_plink(args.outfile, test_associaTR_dir + "/combined.plink2.trait_0.glm.linear", args.phenotype_name)

