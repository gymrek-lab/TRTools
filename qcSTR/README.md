# qcSTR: a tool for generatng quality control plots for TR genotypes
<a href="#usage">Usage</a> | <a href="#outputs">Outputs</a> | <a href="#examples">Example</a> 

qcSTR is a tool for generating various plots that are useful for diagnosing issues in TR calling.

<a name="usage"></a>
## Usage
qcSTR takes as input a VCF file and outputs several plots in pdf format. To run qcSTR, use the following command:

```
qcSTR \
   --vcf <vcf file> \
   --out test \
   [additional options]
```

Required Parameters:
* `--vcf <string>`: Input VCF file
* `--out <string>`: Prefix to name output files

Recommended Parameters:
* `--vcftype <string>`: Type of VCF files being merged. Default='auto'. Must be one of: 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'.

Optional Parameters:
* `--samples <string>`: File containing list of samples to include. If not specified, all samples are used.
* `--period <int>`: Restrict to TRs with this motif length. e.g. to restrict to dinucleotide repeats, use `--period 2`.

<a name="outputs"></a>
## Outputs

qcSTR outputs the following plots:
* outprefix+"xx-sample-callnum.pdf": a barplot giving the number of calls for each sample. Can be used to determine failed or outlier samples.
* outprefix+"xx-chrom-callnum.pdf": a barplot giving the number of calls for each chromosome. Can be useful to determine if the expected number of calls per chromosome are present.
* outprefix+"-diffref-histogram.pdf": a histogram of the difference from the reference allele (in number of repeat units) for each allele called. Can be used to visualize if there is a strong bias toward calling deletions vs. insertions compared to the reference, which might indicate a problem.
* outprefix+"xx-diffref-bias.pdf": plots reference length (bp) vs. the mean difference in length of each allele called compared to the reference allele. It is expected that the mean difference should be around 0 for most settings. When this value starts to deviate from 0, e.g. for very long repeats, it could indicate a drop in call quality.

<a name="examples"></a>
## Example qcSTR command

```
FILE=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
qcSTR \
  --vcf ${FILE1} \
  --out test-qc
```
where `$REPODIR` points to the root path of this repository.
