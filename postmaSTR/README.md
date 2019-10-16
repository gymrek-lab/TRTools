# PostmaSTR 

PostmaSTR is a tool for post-processing and priotization of STR genotypes

## Usage 
To run postmaSTR use the following command: 
```
```

### Parameters 

Required parameters: 
* **`--vcf`**: input the STR vcf file 
* **`--fam`**: input the fam file 
* **`--out`**: input the prefix for the output files 

Locus-level filters: 
* **`--use-length`**: calculates the per-locus stats (het, HWE) collapsing alleles by length 

Call-level filters specific to affected samples: 
* **`--affec-max-expansion-prob-het`**: filters calls with probability of heterozygous expansion less than this 
* **`--affec-min-expansion-prob-het`**: filters calls with probability of heterozygous expansion more than this 
* **`--affec-max-expansion-prob-hom`**: filters calls with probability of homozygous expansion less than this
* **`--affec-min-expansion-prob-hom`**: filters calls with probability of homozygous expansion more than this
* **`--affec-max-expansion-prob-total`**: filters calls with probability of total expansion less than this
* **`--affec-min-expansion-prob-total`**: filters calls with probability of total expansion more than this
* **`--affec-min-call-count`**: this is the minimum number of affected calls per TR region

Call-level filters specific to unaffected samples:
* **`--unaff-max-expansion-prob-het`**: filters calls with probability of heterozygous expansion less than this
* **`--unaff-min-expansion-prob-het`**: filters calls with probability of heterozygous expansion more than this
* **`--unaff-max-expansion-prob-hom`**: filters calls with probability of homozygous expansion less than this
* **`--unaff-min-expansion-prob-hom`**: filters calls with probability of homozygous expansion more than this
* **`--unaff-max-expansion-prob-total`**: filters calls with probability of total expansion less than this
* **`--unaff-min-expansion-prob-total`**: filters calls with probability of total expansion more than this
* **`--unaff-min-call-count`**: this is the minimum number of unaffected calls per TR region

Debugging parameters: 
* **`--num-records`**: this only processes this many records
* **`--die-on-warning`**: quits if a record can't be parsed 
* **`--verbose`**: prints out extra info 


