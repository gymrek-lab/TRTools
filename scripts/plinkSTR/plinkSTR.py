#!/usr/bin/env python
"""
Perform STR association tests

Example - KD:
./plinkSTR.py --vcf /storage/mgymrek/KD/imputed/KD.chr22.imputed.vcf.gz --out stdout --fam /storage/mgymrek/KD/pheno/KD.fam --sex --logistic --infer-snpstr --allele-tests --allele-tests-length  --remove-rare-str-alleles 0.05 --region 22:17655257-17655258

Example - PGC:
./plinkSTR.py \
--vcf /home/gymrek/ssaini/per_locus/11.57523883/PGC.11.57523883.imputed.noAF.vcf.gz \
--covar /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov \
--covar-name C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 \
--out stdout \
--fam /home/gymrek/pgc_imputation/PGC_eur.fam \
--sex --logistic --infer-snpstr --allele-tests --allele-tests-length \
--remove-rare-str-alleles 0.05
"""

# TODO
# - categorical variables (e.g. PGC cohort)
# - linear regression
# - add standard errors
# - make it faster...
# - check for VIF

# Constants
MIN_STR_LENGTH = 8 # min ref length for an STR

# Imports
import sys
sys.path.append("/home/gymrek/workspace/STRTools/scripts/utils")

import warnings
warnings.filterwarnings("ignore")

import argparse
import common
import numpy as np
import pandas as pd
from statsmodels.formula.api import logit
import vcf

def LoadCondition(vcffile, condition, sample_order):
    reader2 = vcf.Reader(open(vcffile, "rb"))
    chrom, start = condition.split(":")
    region = "%s:%s-%s"%(chrom, start, int(start)+1)
    reader2.fetch(region)
    for record in reader2:
        print record.start, int(start), record.ID
        if record.start==int(start):
            return LoadGT(record, sample_order, is_str=False)
    common.ERROR("Could not find SNP to condition on")

def GetAssocType(is_str, alt=-1, alt_len=-1, name=None):
    """
    Return string describing association type
    """
    if not name is None: return name
    if not is_str: return "SNP"
    else:
        if alt >= 0: return "STR-alt-%s"%alt
        elif alt_len >= 0: return "STR-length-%s"%alt_len
        else: return "STR"

def PrintHeader(outf, case_control=False, quant=True, comment_lines=[]):
    """
    Print header info for association output

    Input:
    - outf (file handle): output file handle
    - case_control (bool): specific logistic regression output
    - quant (bool): specify linear regression output
    """
    if case_control:
        header = ["chrom", "start", "type", "p-val", "OR", "se(OR)", "CI95", "maf", "N"]
    else: header = ["chrom", "start", "type", "p-val", "coef", "stderr", "CI95", "maf", "N"]
    for cl in comment_lines:
        outf.write("#"+cl.strip()+"\n")
    outf.write("\t".join(header)+"\n")
    outf.flush()

def OutputAssoc(chrom, start, assoc, outf, assoc_type="STR"):
    """
    Write association output

    Input:
    - chrom (str)
    - start (int)
    - assoc (dict): contains association results
    - outf (file handle): output file handle
    - assoc_type (str): type of association
    """
    if assoc is None: return
    items = [chrom, start, assoc_type, assoc["pval"], assoc["coef"], assoc["stderr"], assoc["CI"], assoc["maf"], assoc["N"]]
    outf.write("\t".join([str(item) for item in items])+"\n")
    outf.flush()

def PerformAssociation(data, covarcols, case_control=False, quant=True, minmaf=0.05, exclude_samples=[], maxiter=100):
    """
    Perform association tests

    Input:
    - data (pd.DataFrame): has columns GT, phenotype, and covarcols
    - covarcols (list<str>): names of columns to use as covars
    - case_control (bool): indicate to perform logistic regression
    - quant (bool): indicate to perform linear regression
    - minmaf (float): don't attempt regression below this MAF

    Output:
    - assoc (dict). Returns association results. Return None on error
    """
    assoc = {}
    formula = "phenotype ~ GT+"+"+".join(covarcols)
    maf = sum(data["GT"])*1.0/(2*data.shape[0])
    assoc["maf"] = "%.3f"%maf
    assoc["N"] = data.shape[0]
    if minmaf != 1 and (maf <= minmaf or (maf >= 1-minmaf)):
        return None # don't attempt regression
    if case_control:
        try:
            pgclogit = logit(formula=formula, data=data[data["sample"].apply(lambda x: x not in exclude_samples)][["phenotype", "GT"]+covarcols]).fit(disp=0, maxiter=maxiter)
            #print pgclogit.summary() # TODO remove after debug
        except:
            assoc["coef"] = "NA"
            assoc["pval"] = "NA"
            assoc["stderr"] = "NA"
            assoc["CI"] = "NA"
            return assoc
        assoc["coef"] = "%.3f"%np.exp(pgclogit.params["GT"])
#        try:
        assoc["pval"] = "%.2E"%pgclogit.pvalues["GT"]
        assoc["stderr"] = "%.3f"%(np.exp(pgclogit.params["GT"])*pgclogit.bse["GT"])
        assoc["CI"] = "-".join(["%.3f"%(np.exp(pgclogit.conf_int().loc["GT",cind])) for cind in [0, 1]])
#        except:
#            assoc["pval"] = "NA"
#            assoc["stderr"] = "NA"
#            assoc["CI"] = "NA"
    else:
        return None # TODO implement linear
    return assoc

def LoadGT(record, sample_order, is_str=True, use_alt_num=-1, use_alt_length=-1, rmrare=0):
    """
    Load genotypes from a record and return values in the sample order

    Input:
    - record (vcf._Record): input record
    - sample_order (list<str>): list of sample ids. Return genotypes in this order
    - is_str (bool): If false, treat as a SNP and use GT field. 
                     If true, treat as STR and by default use length
    - use_alt_num (int): If >=0, treat as bi-allelic using this allele number as the reference
    - use_alt_length (int): If >=0, treat as bi-allelic using this allele length as the reference

    Output:
    - genotypes (list<int>): list of genotype values using given sample order
    - exclude_samples: list of sample to exclude later
    """
    gtdata = {}
    exclude_samples = []
    alleles = [record.REF]+record.ALT
    afreqs = [(1-sum(record.aaf))]+record.aaf
    for sample in record:
        if not is_str: # Note, code opposite to match plink results
            gtdata[sample.sample] = sum([1-int(item) for item in sample.gt_alleles])
        else:
            if use_alt_num > -1:
                gtdata[sample.sample] = sum([int(int(item)==use_alt_num) for item in sample.gt_alleles])
            elif use_alt_length >-1:
                gtdata[sample.sample] = sum([int(len(alleles[int(item)])==use_alt_length) for item in sample.gt_alleles])
            else:
                f1, f2 = [afreqs[int(item)] for item in sample.gt_alleles]
                if f1 < rmrare or f2 < rmrare: 
                    exclude_samples.append(sample.sample)
                gtdata[sample.sample] = sum([len(alleles[int(item)]) for item in sample.gt_alleles])
    d = [gtdata[s] for s in sample_order]
    return d, exclude_samples

def RestrictSamples(data, samplefile, include=True):
    """
    Include or exclude specific samples

    Input:
    - data (pd.DataFrame): data frame. must have columns FID, IID
    - samplefile (str): filename of samples. Must have two columns (FID, IID)
    - include (bool): If true, include these samples. Else exclude them

    Output:
    - data (pd.DataFrame): modified dataframe
    """
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
    default_cols = ["FID", "IID"]
    other_defaults = ["Father_ID", "Mother_ID", "sex","phenotype"]
    if covar_name:
        colnames = default_cols+covar_name.split(",")
        cov = pd.read_csv(fname, delim_whitespace=True, \
                          usecols=colnames)
    elif covar_number:
        colnames = default_cols+["C"+item for item in covar_number.split(",")]
        cov = pd.read_csv(fname, delim_whitespace=True, \
                          names=colnames, \
                          usecols = list(range(2))+[1+int(item) for item in covar_number.split(",")])
    else:
        cov = pd.read_csv(fname, delim_whitespace=True)
        if "FID" not in cov.columns: cov.columns = default_cols+cov.columns[len(default_cols):]
    data = pd.merge(data, cov, on=["FID","IID"])
    covarcols = [item for item in data.columns if item not in default_cols and item not in other_defaults]
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
    data["phenotype"] = data["phenotype"].apply(int)-1 # convert to 0/1
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
    covar_group.add_argument("--cohort-pgc", help="Use cohort from PGC FIDs as a covariate", action="store_true")
    assoc_group = parser.add_argument_group("Association testing")
    assoc_group.add_argument("--linear", help="Perform linear regression", action="store_true")
    assoc_group.add_argument("--logistic", help="Perform logistic regression", action="store_true")
    assoc_group.add_argument("--region", help="Only process this region (chrom:start-end)", type=str)
    assoc_group.add_argument("--infer-snpstr", help="Infer which positions are SNPs vs. STRs", action="store_true")
    assoc_group.add_argument("--allele-tests", help="Also perform allele-based tests using each separate allele", action="store_true")
    assoc_group.add_argument("--allele-tests-length", help="Also perform allele-based tests using allele length", action="store_true")
    assoc_group.add_argument("--minmaf", help="Ignore bi-allelic sites with low MAF", type=float, default=0.05)
    assoc_group.add_argument("--str-only", help="Used with --infer-snptr, only analyze STRs", action="store_true")
    assoc_group.add_argument("--remove-rare-str-alleles", help="Remove genotypes with alleles less than this freq", default=0.0, type=float)
    assoc_group.add_argument("--max-iter", help="Maximum number of iterations for logistic regression", default=100, type=int)
    fm_group = parser.add_argument_group("Fine mapping")
    fm_group.add_argument("--condition", help="Condition on this position chrom:start", type=str)
    args = parser.parse_args()
    # Some initial checks
    if int(args.linear) + int(args.logistic) != 1: ERROR("Must choose one of --linear or --logistic")

    # Load phenotype information
    common.MSG("Loading phenotype information...")
    if args.fam is not None:
        pdata = LoadPhenoData(args.fam, fam=True, missing=args.missing_phenotype, sex=args.sex)
    elif args.pheno is not None:
        if args.sex: ERROR("--sex only works when using --fam (not --pheno)")
        pdata = LoadPhenoData(args.pheno, fam=False, missing=args.missing_phenotype, mpheno=args.mpheno)
    else:
        common.ERROR("Must specify phenotype using either --fam or --pheno")
    common.MSG("Loaded %s samples..."%pdata.shape[0])

    # Load covariate information
    common.MSG("Loading covariate information...")
    covarcols = []
    if args.covar is not None:
        pdata, covarcols = AddCovars(pdata, args.covar, args.covar_name, args.covar_number)
    if args.sex: covarcols.append("sex")
    if args.cohort_pgc:
        pdata["cohort"] = pdata["FID"].apply(lambda x: x.split("*")[0])
        covarcols.append("cohort")
    common.MSG("Loaded %s samples..."%pdata.shape[0])

    # Include/exclude samples
    common.MSG("Loading sample information...")
    if args.samples is not None:
        pdata = RestrictSamples(pdata, args.samples, include=True)
    if args.exclude_samples is not None:
        pdata = RestrictSamples(pdata, args.exclude_samples, include=False)
    common.MSG("Left with %s samples..."%pdata.shape[0])

    # Setup VCF reader
    common.MSG("Set up VCF reader")
    reader = vcf.Reader(open(args.vcf, "rb"))

    # Set sample ID to FID_IID to match vcf
    common.MSG("Set up sample info")
    pdata["sample"] = pdata.apply(lambda x: x["FID"]+"_"+x["IID"], 1)
    reader_samples = set(reader.samples)
    pdata = pdata[pdata["sample"].apply(lambda x: x in reader_samples)]
    sample_order = list(pdata["sample"])
    pdata = pdata[["phenotype","sample"]+covarcols]
    common.MSG("Left with %s samples..."%pdata.shape[0])

    # Get data to condition on
    if args.condition is not None:
        cond_gt = LoadCondition(args.vcf, args.condition, sample_order)
        pdata["condition"] = cond_gt[0]
        covarcols.append("condition")

    # Prepare output file
    if args.out == "stdout":
        outf = sys.stdout
    else: outf = open(args.out, "w")
    PrintHeader(outf, case_control=args.logistic, quant=args.linear, comment_lines=[" ".join(sys.argv)])

    # Perform association test for each record
    common.MSG("Perform associations... with covars %s"%str(covarcols))
    if args.region: reader = reader.fetch(args.region)
    for record in reader:
        # Check MAF 
        aaf = sum(record.aaf)
        aaf = min([aaf, 1-aaf])
        if aaf < args.minmaf: continue
        # Infer whether we should treat as a SNP or STR
        is_str = True # by default, assume all data is STRs
        if args.infer_snpstr:
            if len(record.REF)==1 and len(record.ALT)==1 and len(record.ALT[0])==1:
                is_str = False
            if is_str and len(record.REF) < MIN_STR_LENGTH: continue # probably an indel
            if not is_str and args.str_only: continue
        # Extract genotypes in sample order, perform regression, and output
        common.MSG("   Load genotypes...")
        gts, exclude_samples = LoadGT(record, sample_order, is_str=is_str, rmrare=args.remove_rare_str_alleles)
        pdata["GT"] = gts
        if is_str: minmaf = 1
        else: minmaf = args.minmaf
        common.MSG("   Perform association...")
        assoc = PerformAssociation(pdata, covarcols, case_control=args.logistic, quant=args.linear, minmaf=minmaf, exclude_samples=exclude_samples, maxiter=args.max_iter)
        common.MSG("   Output association...")
        OutputAssoc(record.CHROM, record.POS, assoc, outf, assoc_type=GetAssocType(is_str, name=record.ID))
        # Allele based tests
        common.MSG("   Allele based tests...")
        if is_str and args.allele_tests:
            alleles = [record.REF]+record.ALT
            for i in range(len(record.ALT)+1):
                gts, exclude_samples = LoadGT(record, sample_order, is_str=True, use_alt_num=i)
                pdata["GT"] = gts
                assoc = PerformAssociation(pdata, covarcols, case_control=args.logistic, quant=args.linear, exclude_samples=exclude_samples, maxiter=args.max_iter)
                OutputAssoc(record.CHROM, record.POS, assoc, outf, assoc_type=GetAssocType(is_str, alt=alleles[i]))
        if is_str and args.allele_tests_length:
            for length in set([len(record.REF)] + [len(alt) for alt in record.ALT]):
                gts, exclude_samples = LoadGT(record, sample_order, is_str=True, use_alt_length=length)
                pdata["GT"] = gts
                assoc = PerformAssociation(pdata, covarcols, case_control=args.logistic, quant=args.linear, exclude_samples=exclude_samples, maxiter=args.max_iter)
                OutputAssoc(record.CHROM, record.POS, assoc, outf, assoc_type=GetAssocType(is_str, alt_len=length))
        
if __name__ == "__main__":
    main()
