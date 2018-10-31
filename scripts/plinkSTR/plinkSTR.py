#!/usr/bin/env python
"""
Perform STR association tests

Example:
./plinkSTR.py \
--vcf /storage/mgymrek/KD/imputed/KD.chr22.imputed.vcf.gz \
--out /storage/mgymrek/KD/assoc/KD.chr22.assoc.tab \
--fam /storage/mgymrek/KD/pheno/KD.fam \
--sex \
--logistic \
--infer-snpstr \
--allele-tests

"""

# Constants
MIN_STR_LENGTH = 8 # min ref length for an STR

# Imports
import sys
sys.path.append("../utils")

import warnings
warnings.filterwarnings("ignore")

import argparse
import common
import numpy as np
import pandas as pd
import vcf

def RestrictSamples(data, samplefile, include=True):
    samples = pd.read_csv(samplefile, names=["FID", "IID"])
    if include:
        data = pd.merge(data, samples, on=["FID", "IID"])
    else:
        data = pd.merge(data, samples, on=["FID", "IID"], how="left", indicator=True)
        data = data[data["_merge"]=="left_only"]
        data = data.drop("_merge", 1)
    return data

def AddCovars(data, fname, covar_name, covar_number):
    """
    Add covariates to phenotype data frame and return names of covar columns
    """
    default_cols = ["FID", "IID", "Father_ID", "Mother_ID"]
    if covar_name:
        colnames = default_cols+covar_name.split(",")
        cov = pd.read_csv(fname, delim_whitespace=True, \
                          names=colnames, usecols=colnames)
    elif covar_number:
        colnames = default_cols+["C"+item for item in covar_number.split(",")]
        cov = pd.read_csv(fname, delim_whitespace=True, \
                          names=colnames, \
                          usecols = list(range(4))+[int(item)-1 for item in covar_number.split(",")])
    else:
        cov = pd.read_csv(fname, delim_whitespace=True)
        if "FID" not in cov.columns: cov.columns = default_cols+cov.columns[len(default_cols):]
    data = pd.merge(data, cov, on=["FID","IID"])
    covarcols = [item for item in data.columns if item not in default_cols]
    return data, covarcols

def LoadPhenoData(fname, fam=True, missing=-9, mpheno=1, sex=False):
    """
    Load phenotype data from fam or pheno file
    Only return samples with non-missing phenotype
    If using sex as a covariate, only return samples with sex specified

    Input:
    - fname (str): input filename
    - fam (bool): True if the file is .fam. Else assume .pheno
    - missing (str): Value for missing phenotypes
    - mpheno (int): If using .pheno file, take phenotype from column 2+mpheno
    - sex (bool): Using sex as a covariate, so remove samples with no sex specified
    """
    if fam:
        data = pd.read_csv(fname, delim_whitespace=True, names=["FID", "IID", "Father_ID", "Mother_ID", "sex", "phenotype"])
        if sex:
            data = data[data["sex"]!=0] # (1=male, 2=female, 0=unknown)
    else:
        data = pd.read_csv(fname, delim_whitespace=True, names=["FID", "IID","phenotype"], usecols=[0,1,1+mpheno])
    data = data[data["phenotype"].apply(str) != missing]
    data["phenotype"] = data["phenotype"].apply(int)
    return data

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input VCF file", type=str)
    inout_group.add_argument("--out", help="Output prefix", type=str)
    inout_group.add_argument("--fam", help="FAM file with phenotype info", type=str)
    inout_group.add_argument("--samples", help="File with list of samples to include", type=str)
    inout_group.add_argument("--exclude-samples", help="File with list of samples to exclude", type=str)
    pheno_group = parser.add_argument_group("Phenotypes")
    pheno_group.add_argument("--pheno", help="Phenotypes file (to use instead of --fam)", type=str)
    pheno_group.add_argument("--mpheno", help="Use (n+2)th column from --pheno", type=int, default=1)
    pheno_group.add_argument("--missing-phenotype", help="Missing phenotype code", type=str, default="-9")
    covar_group = parser.add_argument_group("Covariates")
    covar_group.add_argument("--covar", help="Covariates file", type=str)
    covar_group.add_argument("--covar-name", help="Names of covariates to load. Comma-separated", type=str)
    covar_group.add_argument("--covar-number", help="Column number of covariates to load. Comma-separated", type=str)
    covar_group.add_argument("--sex", help="Include sex from fam file as covariate", action="store_true")
    assoc_group = parser.add_argument_group("Association testing")
    assoc_group.add_argument("--linear", help="Perform linear regression", action="store_true")
    assoc_group.add_argument("--logistic", help="Perform logistic regression", action="store_true")
    assoc_group.add_argument("--region", help="Only process this region (chrom:start-end)", type=str)
    assoc_group.add_argument("--infer-snpstr", help="Infer which positions are SNPs vs. STRs", action="store_true")
    assoc_group.add_argument("--allele-tests", help="Also perform allele-based tests using each separate allele", action="store_true")
    assoc_group.add_argument("--allele-tests-length", help="Also perform allele-based tests using allele length", action="store_true")
    fm_group = parser.add_argument_group("Fine mapping")
    fm_group.add_argument("--condition", help="Comma-separated list of positions (chrom:start) to condition on", type=str)
    args = parser.parse_args()

    # Load phenotype information
    if args.fam is not None:
        pdata = LoadPhenoData(args.fam, fam=True, missing=args.missing_phenotype, sex=args.sex)
    elif args.pheno is not None:
        if args.sex: ERROR("--sex only works when using --fam (not --pheno)")
        pdata = LoadPhenoData(args.pheno, fam=False, missing=args.missing_phenotype, mpheno=args.mpheno)
    else:
        common.ERROR("Must specify phenotype using either --fam or --pheno")

    # Load covariate information
    if args.covar is not None:
        pdata, covarcols = AddCovars(pdata, args.covar, args.covar_name, args.covar_number)
    elif args.sex is not None: covarcols = ["sex"]
    else: covarcols = []

    # Include/exclude samples
    if args.samples is not None:
        pdata = RestrictSamples(pdata, args.samples, include=True)
    if args.exclude_samples is not None:
        pdata = RestrictSamples(pdata, args.exclude_samples, include=False)

    # Set sample ID to FID_IID to match vcf
    pdata["sample"] = pdata.apply(lambda x: x["FID"]+"_"+x["IID"], 1)
    sample_order = list(pdata["sample"])

    # Perform association test for each record
    reader = vcf.Reader(open(args.vcf, "rb"))
    if args.region: reader = reader.fetch(args.region)
    for record in reader:
        # Infer whether we should treat as a SNP or STR
        is_str = True # by default, assume all data is STRs
        if args.infer_snpstr:
            if len(record.REF)==1 and len(record.ALT)==1 and len(record.ALT[0])==1:
                is_str = False
            if len(record.REF) < MIN_STR_LENGTH: continue # probably an indel
        # Extract genotypes in sample order - TODO
        # Perform regression - TODO
        # Output results - TODO
        print(record)
        
if __name__ == "__main__":
    main()
