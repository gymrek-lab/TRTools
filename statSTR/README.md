# StatSTR 

This tool is used for computing stats on VCF files 

## Usage 
To run statSTR use the following command: 
```
statSTR 
  --vcf /storage/mgymrek/gangstr-analysis/round2/trio/trio_gangstr_filtered_level1.vcf.gz \
  --out test.tab \
  --thresh \
  [filter options]
```

Required Parameters: 
* **`--vcf`**: input the STR VCF file 
* **`--out`**: prefix to name output files

Filtering Group parameters: 
* **`--samples`**: this is the file containing the list of samples to include 
* **`--region`**: restrict to specific regions (chrom:start-end) 

Stats group parameters: 
* **`--thres`**: an output threshold field; threshold is set to the max observed allele length
* **`--afreq`**: to output allele frequencies 
* **`--acount`**: to output allele counts 
* **`--hwep`**: output HWE p-values per loci 
* **`--het`**: output observed heterozygote counts used for HWE per loci 
* **`--use-length`**: calculate per-locus stats (het, HWE) collapsing alleles by length 


