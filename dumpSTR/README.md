# DumpSTR : A tool for filtering tandem repeat VCF files

[Usage](#usage) | [Filter options](#filters) | [Output files](#outputs) 

DumpSTR is a tool for filtering VCF file with TR genotypes produced by supported genotyping tools. It can perform both call-level and locus-level filtering and outputs a filtered VCF file.

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


<a name="outputs"></a>
## Output files

DumpSTR outputs the following files:

* `$out.vcf`: Filtered VCF file. Filtered loci have a list of failing filters in the FILTER column. An additional FORMAT:FILTER field is added to each call. This is set to PASS for passing calls. For failing calls, this is set to a list of filter reasons and the genotype is set to missing.
* `$out.samplog.tab`: Output sample-level log info. This is a tab-delimited file with columns: sample, number of calls, and mean coverage at that sample.
* `$out.loclog.tab`: Output locus-level log info. It contains the mean call rate at passing TR loci. It also contains a separate line for each filter with the number of TR loci failing that filter.
