# DumpSTR

[Usage](#usage) | [Filter options](#filters) | [Output files](#outputs) | [Recommended filters](#recommended)

DumpSTR is a tool for filtering VCF file with STR genotypes produced by supported TR genotyping tools. It can perform both call-level and locus-level filtering and outputs a filtered VCF file.

<a name="usage"></a>
## Usage
To run dumpSTR use the following command:
```
dumpSTR \
  --vcf <VCF> \
  --out <string> \
  [filter options]
```

Required parameters:
* `--vcf <VCF>` VCF file output by a supported genotyping tool.
* `--vcftype <string>`: Type of VCF files being merged. Default='auto'. Must be one of: 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'.
* `--out <string>` prefix to name output files.

DumpSTR will output a new VCF file named `$out.vcf`, a sample log file `$out.samplog.tab`, and a locus log file `$out.loclog.tab`. See a description of output files below.

Specific filters available are described below.

<a name="filters"></a>
## Filter options

DumpSTR offers the following types of filters:

* [Locus-level filters](#locus)
* [Call-level filters](#call)
  * [GangSTR](#gangstr)
  * [HipSTR](#hipstr)
  * [ExpansionHunter](#eh)
  * [PopSTR](#popstr)
  * [AdVNTR](#advntr)



<a name="locus"></a>
### Locus-level filters

These filters are not specific to any tool and can be applied to any VCF file:

| DumpSTR option | Filter Description |
| ----| ------|
| `--min-locus-callrate <float>` | Filters loci with too few calls |
| `--min-locus-hwep <float>` | Filters loci departing from Hardy Weinberg equilibrium at some p-value threshold. Based on a two-sided binomial test comparing the observed vs. expected percent of calls that are homozygous |
| `--min-locus-het <float>` | Filters loci with low heteroyzgosity. Where heterozygosity is equal to: 1-sum_i p_i^2, where p_i is the frequency of allele i |
| `--max-locus-het <float>` | Filters loci with high heterozygosity |
| `--use-length` | Use allele lengths, rather than sequences, to compute heterozygosity and HWE (only relevant for HipSTR, which reports sequence level differences in TR alleles) |
| `--filter-regions <BEDFILE[,BEDFILE12,...]>` | Filter TRs overlapping the specified set of regions. Must be used with `--filter-regions-names`. Can supply a comma-separated list to each to apply multiple region filters. Bed files must be sorted and tabix-indexed. |
| `--filter-regions-names <string[,string2,...]>` | Filter names for each BED file specified in `--filter-regions`. |
| `--filter-hrun` | Filter repeats with long homopolymer runs. Only used for HipSTR VCF files otherwise ignored. Ignores pentanucleotides with homopolymer runs of 5bp or longer, or hexanucleotides with homopolymer runs of 6bp or longer. |
|`--drop-filtered` | Do not output loci that were filtered, or loci with no calls remaining after filtering. |

TRs passing all locus-level filters will be marked as "PASS" in the FILTER field. Those failing will have a list of failing filters in the FILTER field. If `drop-filtered` is specified, only loci passing all filters will be output.

<a name="call"></a>
### Call-level filters

Different call-level filters are available for each supported TR genotyping tool:

<a name="gangstr"></a>
#### GangSTR call-level filters

| DumpSTR option | Filter Description |
|------| -----|
| `--gangstr-min-call-DP <int>` | Minimum call coverage. Based on DP field. |
| `--gangstr-max-call-DP <int>` | Maximum call coverage. Based on DP field. |
| `--gangstr-min-call-Q <float>` | Minimum call quality score. Based on Q field. |
| `--gangstr-expansion-prob-het <float>` | Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this. Based on QEXP field. |
| `--gangstr-expansion-prob-hom <float>` | Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this. Based on QEXP field. |
| `--gangstr-expansion-prob-total <float>` | Expansion prob-value threshold. Filters calls with probability of homozygous  or heterozygous expansion less than this. Based on QEXP field. |
| `--gangstr-filter-span-only` | Filter out all calls that only have spanning read support. Based on RC field. |
| `--gangstr-filter-spanbound-only` | Filter out all reads except spanning and bounding. Based on RC field. | 
| `--gangstr-filter-badCI` | Filter regions where the ML estimate is not in the CI. Based on REPCN and REPCI fields. |
| `--gangstr-require-support <int>` | Require each allele call to have at least this many supporting reads. Based on ENCLREADS, RC, and FLNKREADS fields.|
| `--gangstr-readlen <int>` | Read length used (bp). Required if using --require-support. |

<a name="hipstr"></a>
#### HipSTR call-level filters

| DumpSTR option | Filter Description |
| ----| ------|
| `--hipstr-max-call-flank-indel <float>` | Maximum call flank indel rate. Computed as DFLANKINDEL/DP |
| `--hipstr-max-call-stutter <float>` | Maximum call stutter rate. Computed as DSTUTTER/DP |
| `--hipstr-min-supp-reads <int>` | Minimum supporting reads for each allele. Based on ALLREADS and GB fields |
| `--hipstr-min-call-DP <int>` | Minimum call coverage. Based on DP field. |
| `--hipstr-max-call-DP <int>` | Maximum call coverage. Based on DP field. |
| `--hipstr-min-call-Q <float>` | Minimum call quality score. Based on Q field. |

<a name="popstr"></a>
#### PopSTR call-level filters
| DumpSTR option | Filter Description |
| ----| ------|
| `--popstr-min-call-DP <int>` | Minimum call coverage. Based on DP field. |
| `--popstr-max-call-DP <int>` | Maximum call coverage. Based on DP field. |
| `--popstr-require-support <int>` | Require each allele call to have at least n supporting reads. Based on AD field.|

<a name="eh"></a>
#### ExpansionHunter call-level filters
| DumpSTR option | Filter Description |
| ----| ------|
| `--eh-min-ADFL <int>` | Minimum number of flanking reads consistent with the allele. Based on ADFL field. |
| `--eh-min-ADIR <int>` | Minimum number of in-repeat reads consistent with the allele. Based on ADIR field. |
| `--eh-min-ADSP <int>` | Minimum number of spanning reads consistent with the allele. Based on ADSP field. |
| `--eh-min-call-LC <int>` | Minimum call coverage. Based on LC field. |
| `--eh-max-call-LC <int>` | Maximum call coverage. Based on LC field. |

<a name="advntr"></a>
#### AdVNTR call-level filters
| DumpSTR option | Filter Description |
| ----| ------|
| `--advntr-min-call-DP <int>` | Minimum call coverage. Based on DP field. |
| `--advntr-max-call-DP <int>` | Maximum call coverage. Based on DP field. |
| `--advntr-min-spanning <int>` | Minimum spanning read count (SR field) |
| `--advntr-min-flanking <int>` | Minimum flanking read count (FR field) | 
| `--advntr-min-ML <float>` | Minimum value of maximum likelihood (ML field) |


<a name="recommended"></a>
## Recommended filters

<a name="outputs"></a>
## Output files


## Basic dumpSTR command 

```
dumpSTR \
        --vcf ./test_files/test_files.vcf \
        --out test_run
```

Check test_run_original in Example_files folder for example output. 

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
## Recommended GangSTR filters

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
