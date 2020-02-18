# MergeSTR 
<a href="#usage">Usage</a> | <a href="#example">Example</a> | <a href="#supported">Supported Fields</a> 

MergeSTR is a tool for merging VCF files from supported TR genotyping tools.

If TR genotyping was performed separately on different samples or batches of samples, mergeSTR can be used to combine the resulting VCFs into one file. This is often necessary for downstream steps such as: computing per-locus statistics, performing per-locus filtering, and association testing.

While other VCF libraries have capabilities to merge VCF files, they do not always handle multi-allelic STRs properly, especially if the allele definitions are different across files. MergeSTR is TR-aware and currently handles VCF files obtained by: GangSTR, HipSTR, ExpansionHunter, popSTR, or adVNTR. See below for specific VCF fields supported for each genotyper.

<a name="usage"></a>
## Usage 
MergeSTR takes as input two or more VCF files with TR genotypes and outputs a combined VCF file. Note, input VCF files must be sorted, indexed, and have the appropriate `##contig` header lines.

To run mergeSTR use the following command: 
```
mergeSTR \
  --vcfs <vcf file1, vcf file2> \
  --out test \
  [additional options]
```

Required Parameters: 
* `--vcf <VCF>`: Comma-separated list of VCF files to merge. All must have been created by the same TR genotyper. Must be bgzipped, sorted, and indexed. 
* `--vcftype <string>`: Type of VCF files being merged. Default='auto'. Must be one of: 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'.
* `--out <string>`: prefix to name output files

Special Merge Options: 
* `--update-sample-from-file`: Append file names to sample names. Useful if sample names are repeated across VCF files.

Optional Additional Parameters: 
* `--verbose`: Prints out extra information 
* `--quiet`: Doesn't print out anything 

<a name="example"></a>
## Example MergeSTR command 

If you have installed the STRTools package, you can run an example mergeSTR command using files in this repository.

```
FILE1=${REPODIR}/test/common/sample_vcfs/mergeSTR_vcfs/test_file_gangstr1.vcf.gz
FILE2=${REPODIR}/test/common/sample_vcfs/mergeSTR_vcfs/test_file_gangstr2.vcf.gz
mergeSTR \
   --vcfs ${FILE1},${FILE2} \
   --out test_run
```
where `$REPODIR` points to the root path of this repository.

If you are testing from source, you can run:
```
python mergeSTR.py \
   --vcfs ${FILE1},${FILE2} \
   --out test_run
```

This command should create a file `test_run.vcf` with the merged genotypes.

<a name="supported"></a>

## Supported VCF fields

In addition to proper merging of alleles at multi-allelic sites, MergeSTR supports the following VCF fields for each tool. Fields not listed are currently ignored when merging. INFO fields below are expected to be constant across loci being merged.

| Tool | Supported INFO fields | Supported FORMAT fields |
|-----|---------|--------|
| GangSTR | END, RU, PERIOD, REF, EXPTHRESH | DP,Q,REPCN,REPCI,RC,ENCLREADS,FLNKREADS,ML,INS,STDERR,QEXP |
| HipSTR | START, END, PERIOD| GB,Q,PQ,DP,DSNP,PSNP,PDP,GLDIFF,DSTUTTER,DFLANKINDEL,AB,FS,DAB,ALLREADS,MALLREADS|
| ExpansionHunter | END, REF, REPID, RL, RU, SVTYPE |ADFL,ADIR,ADSP,LC,REPCI,REPCN,SO |
| PopSTR | Motif| AD, DP, PL|
| AdVNTR |END, RU, RC |DP, SR, FL, ML |
