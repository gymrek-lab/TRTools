# DumpSTR

DumpSTR is a tool for filtering VCF file with STR genotypes produced by HipSTR or GangSTR.

[Usage](#usage) | [Recommended HipSTR filters](#hipstr) | [Recommended GangSTR filters](#gangstr)


<a name="usage"></a>
## Usage
To run dumpSTR use the following command:
```
dumpSTR \
  --vcf <HipSTR or GangSTR VCF> \
  --out <string> \
  [filter options]
```

Required parameters:
* **--vcf <GangSTR or HipSTR VCF>** VCF file output by HipSTR or GangSTR. Can be unzipped or bgzipped.
* **--out <string>** prefix to name output files

Parameters for locus-level filtering:
* **`--min-locus-callrate <float>`**: Filters loci with too few calls
* **`--min-locus-hwep <float>`**: Filters loci departing from Hardy Weinberg equilibrium at some p-value threshold
* **`--min-locus-het <float>`**: Filters loci with low heteroyzgosity
* **`--max-locus-het <float>`**: Filters loci with high heterozygosity
* **`--use-length`**: Use allele lengths, rather than sequences, to compute heterozygosity and HWE (only relevant for HipSTR)
* **`--filter-regions <BEDFILE[,BEDFILE12,...]>`**: Filter loci overlapping the specified set of regions. Must be used with `--filter-regions-names`. Can supply a comma-separated list to each to apply multiple region filters.
* **`--filter-regions-names <string[,string2,...]>`**: Filter names for each BED file specified in `--filter-regions`.
* **`--filter-hrun`**: Filter repeats with long homopolymer runs
* **`--drop-filtered`**: Do not output loci that were filtered, or loci with no calls remaining after filtering.

A file is provided in `dumpSTR/filter_files/hg19_segmentalduplications.bed.gz` for filtering segmental duplications in hg19.

General parameters for call-level filtering (apply to both HipSTR or GangSTR):
* **`--min-call-DP <int>`**: Minimum call coverage
* **`--max-call-DP <int>`**: Maximum call coverage
* **`--min-call-Q <int>`**: Minimum Q value

HipSTR-specific call-level filters:
* **`--max-call-flank-indel <float>`**: Maximum rate of indels in flanking regions
* **`--max-call-stutter <float>**`: Maximum call stutter rate
* **`--min-supp-reads <int>**`: Require this many supporting reads for each allele length called.

GangSTR-specific call-level filters:
* **`--expansion-prob-het <float>`**: Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this
* **`--expansion-prob-hom <float>`**: Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this
* **`--expansion-prob-total <float>`**: Expansion prob-value threshold. Filters calls with probability of total expansion less than this
* **`--filter-span-only`**: Filter loci where only spanning reads were identified
* **`--filter-spanbound-only`**: Filter loci where only spanning or bounding reads were identified
* **`--filter-badCI`**: Filter regions where the ML estimate is not in the CI

<a name="hipstr"></a>
## Recommended HipSTR filters

The following command gives recommended filters for HipSTR calls:

```
dumpSTR \
    --vcf ${HIPVCF} \
    --out ${OUTPREFIX} \
    --filter-regions dumpSTR/filter_files/hg19_segmentalduplications.bed.gz \
    --filter-regions-names SEGDUP \
    --min-supp-reads 1 \
    --min-call-DP 10 \
    --max-call-DP 1000 \
    --min-call-Q 0.9 \
    --max-call-flank-indel 0.15 \
    --max-call-stutter 0.15
```

For population-level datasets (e.g. 50+ samples), the Hardy-Weinberg filter is recommended:
```
--min-locus-hwep 0.01
```

<a name="gangstr"></a>
## Recommended GanSTR filters

The following command gives recommended filters for dumpSTR calls:

```
dumpSTR \
    --vcf ${GANGSTRVCF} \
    --out ${OUTPREFIX} \
    --max-call-DP 1000 \
    --filter-spanbound-only \
    --filter-badCI \
    --filter-regions dumpSTR/filter_files/hg19_segmentalduplications.bed.gz \
    --filter-regions-names SEGDUP
```

Additional filters are recommended based on the application of interest. For more info, see (Filtering GangSTR output)[https://github.com/gymreklab/GangSTR/wiki/Filtering-GangSTR-output].