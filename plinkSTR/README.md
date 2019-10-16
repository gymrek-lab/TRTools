# PlinkSTR 

Perform STR association tests 

## Usage 
To run plinkSTR-KD use the following command: 
```
./plinkSTR.py \
  --vcf /storage/mgymrek/KD/imputed/KD.chr22.imputed.vcf.gz \
  --out stdout --fam /storage/mgymrek/KD/pheno/KD.fam \
  --sex --logistic --infer-snpstr --allele-tests --allele-tests-length \
  --remove-rare-str-alleles 0.05 --region 22:17655257-17655258
```

To run plinkSTR-PGC use the following command: 
```
./plinkSTR.py \
  --vcf /home/gymrek/ssaini/per_locus/11.57523883/PGC.11.57523883.imputed.noAF.vcf.gz \
  --covar /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov \
  --covar-name C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 \
  --out stdout \
  --fam /home/gymrek/pgc_imputation/PGC_eur.fam \
  --sex --logistic --infer-snpstr --allele-tests --allele-tests-length \
  --remove-rare-str-alleles 0.05
```

### Parameters Options

Input/Output parameters: 
* **`--vcf`**: input vcf files 
* **`--out`**: prefix to name output files 
* **`--fam`**: include the FAM file with phenotype info 
* **`--samples`**: include the file with the list of samples to include 
* **`--exclude-samples`**: include the file with the list of samples to exclude 

Phenotype parameters: 
* **`--pheno`**: include the phenotypes file instead of the FAM files 
* **`--mpheno`**: uses the (n+2)th column of the phenotypes file 
* **`--missing-phenotype`**: input the missing phenotype code 

Covariate parameters: 
* **`--covar`**: include the covariates file 
* **`--covar-name`**: input the names of covariates to load, make sure it's comma-separated 
* **`--covar-number`**: input the column number of covariates to load, make sure it's comma-separated
* **`--sex`**: include the sex from the FAM file as a covariate 
* **`--cohort-pgc`**: use the cohort from PGC FIDs as a covariate 

Association Testing parameters: 
* **`--linear`**: performs linear regression 
* **`--logistic`**: performs logistic regression 
* **`--region`**: processes only a specific region, input as (chrom:start-end) 
* **`--infer-snpstr`**: infers which positions are SNPs vs. STRs
* **`--allele-tests`**: perofrms allele-based tests using each separate allele
* **`--allele-tests-length`**: performs allele-based tests using allele length 
* **`--minmaf`**: ignores bi-allelic sites with low MAF 
* **`--str-only`**: used with parameter --infer-snpstr to only analyze STRs
* **`--remove-rare-str-alleles`**: removes genotypes with alleles less than the frequency provided 
* **`--max-iter`**: set the maximum number of iterations for logistic regression 

Fine Mapping parameters: 
* **`--condition`**: condition on a specific position (chrom:start) 
