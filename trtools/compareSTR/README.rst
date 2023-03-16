.. overview_directive
.. |compareSTR overview| replace:: CompareSTR compares different TR callsets generated on the same samples against the same reference panel. CompareSTR outputs overall, per-locus, and per-sample concordance information.
.. overview_directive_done

CompareSTR
==========

|compareSTR overview|

Example use cases include:

* Comparing calls to a "ground truth" set, e.g. from capillary electrophoresis data, to find call errors
* Comparing calls for the same tool using different parameter settings to identify differences due to bioinformatic processing
* Comparing calls for different tools. This only works if they used the same set of reference TRs. Please note that compareSTR only compares TRs with matching chromosomes and allele-positions (after ignoring flanking base pairs)

CompareSTR optionally will stratify results based on a user-specified FORMAT field (e.g. depth, or quality score) and/or by repeat motif length.

Note: CompareSTR is designed to be used as a QC tool. While it may be able to pick up certain biological differences in some applications (e.g. identifying de novo mutations by comparing parent and child callsets or somatic mutations by comparing callsets from different tissues), use-case specific analyses may be better performed by more specialized tools.

Note: CompareSTR has the ability to stratify comparisons based on quality scores. However, beware that quality scores output by different genotypers may not be directly comparable. You can use `qcSTR <https://trtools.readthedocs.io/en/latest/source/qcSTR.html>`_ to visualize the distribution of quality scores in each VCF file seprately.

Outputs
-------

CompareSTR calculates the following metrics in the overall summary, per-locus and per-sample:

* fraction_concordant_len: fraction of genotypes with the same length in the two VCF files. We report this metric because it is intuitive and we expect many users would want it.
  However, as a measurement of raw accuracy this statistic is not very useful, as it is highly influenced by allele frequencies. For example, if an STR has two possible alleles,
  one common (99% frequency) and one rare (1% frequency), then a caller which always calls the common allele would have 99% accuracy, which taken out of context would seem to be a good accuracy,
  even though it is actually guessing at random.
* fraction_concordant_seq: fraction of genotypes with the same sequence in the two VCF files. Currently HipSTR is the only caller which can report alleles
  with the same length but different sequences. Otherwise will be identical to fraction_concordant_len. This can be misleading similarly to fraction_concordant_len.
* r: The Pearson correlation between the sum of allele lengths at each call compared between the two VCF files, where allele lengths are measured as number of repeat copies different from the reference.

CompareSTR calculates the following metric only at a per-locus level:

* Mean absolute difference: for each sample at a locus, for each VCF take the sum of repeat lengths across both chromosomes in that VCF. Take the absolute difference
  of those sums between VCFs, and then report the average of those values. Like fraction_concordant_len, this is a measurement of raw accuracy and so can be misleading.
  This metric is useful for highly multiallelic STRs where all alleles are individually relatively infrequent (so as to not be particularly misleading), and
  when it is useful to know how close the two callers are to one another, as opposed to how often they are match exactly.

These metrics and numcalls only reflect the (sample, locus) pairs that were called by both callers. Pairs which were not called by both callers are omitted instead of being counted as 
disconcordant.

Note that phased genotypes are comapred as a whole instead of per allele for generating the concordance metrics,
so 1|2 compared to 1|2 only generates one successful 'match' instead of 2 and 1|2 compared to 3|2
only generates one unsuccessful 'match' instead of one successful one and one unsuccessful one.

compareSTR outputs the following text files and plots:

* :code:`<outprefix>-overall.tab`: Has columns period, fraction_concordant_seq, fraction_concordant_len, r, numcalls, n_missing_only_vcf1, n_missing_only_vcf2, n_missing_both. Plus additional columns for any FORMAT fields to stratify results on. This file has one line for all results (period="ALL") and a different line for each period analyzed separately if request by :code:`--period`. If stratifying by format fields, it will have additional lines for each range of values for each of those fields.
* :code:`<outprefix>-bubble-period$period.pdf`: "Bubble" plot, which plots the sum of allele lengths for each call in :code:`--vcf1` vs. :code:`--vcf2`. Allele lengths are given in terms of difference in number of repeat units from the reference. The size of each bubble gives the number of calls at each cooordinate. A seperate plot is output for all TRs (period="ALL") and for each period if requested by :code:`--period`.
* :code:`<outprefix>-locuscompare.tab`: Has columns chrom, start, fraction_concordant_seq, fraction_concordant_len, mean_absolute_difference, r, numcalls, n_missing_only_vcf1, n_missing_only_vcf2, n_missing_both. There is one line for each TR.
* :code:`<outprefix>-locuscompare.pdf`: Plots the length concordance metric for each TR locus considered.
* :code:`<outprefix>-samplecompare.tab`: Has columns sample, fraction_concordant_seq, fraction_concordant_len, numcalls, n_missing_only_vcf1, n_missing_only_vcf2, n_missing_both. One line per sample
* :code:`<outprefix>-samplecompare.pdf`: Plots the length concordance metric for each sample considered.

See `Example Commands`_ below for example compareSTR commands for different supported TR genotypers based on example data files in this repository. More detailed use cases are also given in the vignettes https://trtools.readthedocs.io/en/develop/VIGNETTES.html.


Usage
-----
CompareSTR takes as input two VCF files with overlapping TRs and samples and outputs metrics and plots based on comparing calls across the two VCFs. The input VCF files must be sorted, indexed, and have the appropriate `##contig` header lines. CompareSTR only considers the subset of samples shared across two VCF files being compared, based on sample ID in the VCF headers.

Note: if comparing two VCFs with **phased** calls, ensure the chromosome ordering matches. A 'motherAllele|fatherAllele' representation will not match with a 'fatherAllele|motherAllele' representation.

To run compareSTR use the following command::

  compareSTR \
    --vcf1 <file1.vcf.gz> --vcf2 <file2.vcf.gz> \
    --out <string> \
    [additional options]

Required Parameters:

* :code:`--vcf1 <VCF>`: First VCF file to compare (must be sorted, bgzipped, and indexed. See `Instructions on Compressing and Indexing VCF files`_ below)
* :code:`--vcf2 <VCF>`: Second VCF file to compare (must be sorted, bgzipped, and indexed. See `Instructions on Compressing and Indexing VCF files`_ below)
* :code:`--out <string>`: Prefix to name output files

Filtering Options:

* :code:`--samples <string>`: File containing list of samples to include. If not specified, all samples are used.
* :code:`--region <string>`: Restrict to this region chrom:start-end

Options to stratify the overall summary:

* :code:`--stratify-fields`: Comma-separated list of FORMAT fields to stratify by. e.g. DP,Q
* :code:`--stratify-binsizes`: Comma-separated list of min:max:binsize to stratify each field on. Must be same length as :code:`--stratify-fields`. e.g. 0:50:5,0:1:0.1 . The range [min, max] is inclusive.
* :code:`--stratify-file`: Specify which file to look at the :code:`--stratify-fields` in. If set to 0, apply to both files. If set to 1, apply only to :code:`--vcf1`. If set to 2, apply only to :code:`--vcf2`.
* :code:`--period`: Report results overall and also stratified by repeat unit length (period).

Options influencing which statistics are reported or how they are computed:

* :code:`--ignore-phasing`: Treat all calls as if they are unphased, ignoring any phasing information in the VCF. Particularly useful with HipSTR derived calls which are reported as phased even though they are only
  phased locally, not globally, and so the phasing is not directly usable in comparisons.
* :code:`--vcf2-beagle-probabilities`: For use when vcf2 was imputed by Beagle. Use the imputed probabilities in the AP fields instead of the 
* :code:`--balanced-accuracy` : output the balanced accuracy assuming that the first vcf file contains the true genotypes and the second vcf file is being measured against that standard
  (possibly because the second genotype file is imputed). For phased genotypes, the accuracy for each possible allele length is the percent of alleles with that true length (as specified by the first vcf)
  which were assigned that length in the second vcf. The balanced accuracy is the average of the accuracies of each possible allele length. This is different than standard measures of accuracy or concordance
  because each possible allele length contributes equally to this score, regardless of allele frequency imbalances. For nonphased TODO Mathematically TODO

Plotting options:

* :code:`--bubble-min`: Minimum x/y axis value to display on bubble plots.
* :code:`--bubble-max`: Maximum x/y axis value to display on bubble plots.

Other options:

* :code:`--verbose`: Print helpful debugging info
* :code:`--noplot`: Don't output any plots. Only produce text output.
* :code:`--vcftype1 <string>`: Type of VCF file 1.
* :code:`--vcftype2 <string>`: Type of VCF file 2.

Instructions on Compressing and Indexing VCF files
--------------------------------------------------
CompareSTR requires input files to be compressed and indexed. Use the following commands to create compressed and indexed vcf files::

  bgzip file.vcf
  tabix -p vcf file.vcf.gz

Example Commands
----------------

Below are :code:`compareSTR` examples using VCFs from supported TR genotypers. Data files can be found at https://github.com/gymreklab/TRTools/tree/master/example-files::

  # AdVNTR (comparing a file against itself. Not very interesting. Just for demonstration)
  # Note, you first need to reheader files to add required contig lines to VCF headers
  bcftools reheader -f hg19.fa.fai -o NA12878_advntr_reheader.vcf.gz NA12878_chr21_advntr.sorted.vcf.gz
  tabix -p vcf NA12878_advntr_reheader.vcf.gz 
  FILE1=NA12878_advntr_reheader.vcf.gz
  compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out advntr_vs_advntr --noplot

  # HipSTR vs. ExpansionHunter
  compareSTR \
      --vcf1 NA12878_chr21_hipstr.sorted.vcf.gz \
      --vcf2 NA12878_chr21_eh.sorted.vcf.gz \
      --vcftype1 hipstr --vcftype2 eh --out hipstr_vs_eh

  # HipSTR vs. GangSTR
  compareSTR \
      --vcf1 NA12878_chr21_hipstr.sorted.vcf.gz \
      --vcf2 NA12878_chr21_gangstr.sorted.vcf.gz \
      --vcftype1 hipstr --vcftype2 gangstr --out hipstr_vs_gangstr

  # PopSTR (comparing a file against itself. Not very interesting. Just for demonstration)
  FILE1=trio_chr21_popstr.sorted.vcf.gz
  compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out popstr_vs_popstr


