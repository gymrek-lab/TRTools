# CompareSTR: a tool for comparing tandem repeat callsets
<a href="#usage">Usage</a> | <a href="#outputs">Outputs</a> | <a href="#examples">Example</a> 

CompareSTR is a tool for comparing TR callsets. Example use cases include:
* Comparing calls to a "ground truth" set, e.g. from capillary electrophoresis data
* Comparing calls for the same tool using different parameter settings

It outputs overall, per-locus, and per-sample concordance information. It optionally will stratify results based on a user-specified FORMAT field (e.g. depth, or quality score) and by repeat motif length.

<a name="usage"></a>
## Usage
CompareSTR takes as input two VCF files with overlapping TRs and samples and outputs metrics and plots based on comparing calls across the two VCFs. The input VCF files must be sorted, indexed, and have the appropriate `##contig` header lines. The VCF files must also have alleles specified in the same way. e.g. GangSTR-style using full allele sequences, or ExpansionHunter style with SV tags.

To run compareSTR use the following command:
```
compareSTR \
  --vcf1 <vcf file> --vcf2 <vcf file> \
  --out test \
  [additional options]
```

Required Parameters:
* `--vcf1 <VCF>`: First VCF file to compare (must be sorted, bgzipped, and indexed).
* `--vcf2 <VCF>`: Second VCF file to compare (must be sorted, bgzipped, and indexed).
* `--out <string>`: Prefix to name output files

Recommended Parameters:
* `--vcftype <string>`: Type of VCF files being merged. Default='auto'. Must be one of: 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'. If using `auto`, compareSTR will try to figure out the VCF type. This might not work for VCF files not from one of the supported types. e.g. if you are comparing to an experimental ground truth dataset. Therefore you must make sure `--vcf1` and `--vcf2` have alleles specified in the same way, then set `--vcftype` to whatever tool was used to generate `--vcf1` or `--vcf2`.

Filtering Options:
* `--samples <string>`: File containing list of samples to include. If not specified, all samples are used.
* `--region <string>`: Restrict to this region chrom:start-end

Metrics to stratify results:
* `--stratify-fields`: Comma-separated list of FORMAT fields to stratify by. e.g. DP,Q
* `--stratify-binsizes`: Comma-separated list of min:max:binsize to stratify each field on. Must be same length as `--stratify-fields`. e.g. 0:50:5,0:1:0.1
* `--stratify-file`: Specify which file to look at the `--stratify-fields` in. If set to 0, apply to both files. If set to 1, apply only to `--vcf1`. If set to 2, apply only to `--vcf2`.
* `--period`: Report results overall and also stratified by repeat unit length (period).

Plotting options:
* `--bubble-min`: Minimum x/y axis value to display on bubble plots.
* `--bubble-max`: Maximum x/y axis value to display on bubble plots.

Other options:
* `--verbose`: Print helpful debugging info
* `--noplot`: Don't output any plots. Only produce text output.

<a name="outputs"></a>
## Outputs

In output files, compareSTR reports the following metrics:
* Length concordance: % of genotypes concordant between the two VCF files when only considering TR allele lengths
* Sequence concordance: % of genotypes concordant between the two VCF files when considering TR allele sequence. Currently only relevant for HipSTR. Otherwise will be identical to length concordance
* R2: Pearson r2 between the sum of allele lengths at each call compared between the two VCF files.

compareSTR outputs the following text files and plots:

* outprefix+"-overall.tab": Has columns period, concordance-seq, concordance-len, r2, numcalls. Plus additional columns for any FORMAT fields to stratify results on. This file has one line for all results (period="ALL") and a different line for each period analyzed separately. If stratifying by format fields, it will have additional lines for each range of values for each of those fields.
* outprefix+"-bubble-period$period.pdf": "Bubble" plot, which plots the sum of allele lengths for each call in `--vcf1` vs. `--vcf2`. Allele lengths are given in terms of bp difference from the reference genome. The size of each bubble gives the number of calls at each cooordinate. A seperate plot is output for all TRs (period="ALL") and for each period.
* outprefix+"-locuscompare.tab": Has columns chrom, start, metric-conc-seq, metric-conc-len, sample. last column gives the number of samples considered at each locus. There is one line for each TR.
* outprefix+"-locuscompare.pdf": Plots the length concordance metric for each TR locus considered.
* outprefix+"-samplecompare.tab": Has columns sample, metric-conc-seq, metric-conc-len, numcalls. One line per sample
* outprefix+"-samplecompare.pdf": Plots the length concordance metric for each sample considered.

<a name="examples">Examples</a>

Compare two callsets:
```
FILE1=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
FILE2=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf2.vcf.gz
compareSTR \
  --vcf1 ${FILE1} --vcf2 ${FILE2} \
  --out test-compare
```
where `$REPODIR` points to the root path of this repository.

Similarly, to compare two callsets, but stratify by the DP and Q format fields in the first VCF file and output metrics separately by period (and also modify the bubble plot dimensions):
```
FILE1=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
FILE2=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf2.vcf.gz
compareSTR \
  --vcf1 ${FILE1} --vcf2 ${FILE2} \
  --stratify-fields DP,Q \
  --stratify-binsizes 0:50:10,0:1:0.1 \
  --stratify-file 1 \
  --period \
  --bubble-min -50 --bubble-max 50 \
  --out test-compare
```
