#!/usr/bin/env python
"""
Perform STR association tests
"""

import argparse

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input VCF file", type=str)
    inout_group.add_argument("--out", help="Output prefix", type=str)
    pheno_group = parser.add_argument_group("Phenotypes")
    pheno_group.add_argument("--pheno", help="Phenotypes file", type=str)
    pheno_group.add_argument("--mpheno", help="Use (n+2)th column from --pheno", type=int, default=1)
    pheno_group.add_argument("--missing-phenotype", help="Missing phenotype code", type=str, default="-9")
    covar_group = parser.add_argument_group("Covariates")
    covar_group.add_argument("--covar", help="Covariates file", type=str)
    covar_group.add_argument("--covar-name", help="Names of covariates to load. Comma-separated", type=str)
    covar_group.add_argument("--covar-number", help="Column number of covariates to load. Comma-separated", type=str)
    assoc_group = parser.add_argument_group("Association testing")
    assoc_group.add_argument("--linear", help="Perform linear regression", action="store_true")
    assoc_group.add_argument("--logistic", help="Perform logistic regression", action="store_true")
    assoc_group.add_argument("--region", help="Only process this region (chrom:start-end)", type=str)
    assoc_group.add_argument("--samples", help="File with list of samples to include", type=str)
    assoc_group.add_argument("--exclude-samples", help="File with list of samples to exclude", type=str)
    assoc_group.add_argument("--infer-snpstr", help="Infer which positions are SNPs vs. STRs", action="store_true")
    assoc_group.add_argument("--allele-tests", help="Also perform allele-based tests", action="store_true")
    fm_group = parser.add_argument_group("Fine mapping")
    fm_group.add_argument("--condition", help="Comma-separated list of positions (chrom:start) to condition on", type=str)
    args = parser.parse_args()

if __name__ == "__main__":
    main()
