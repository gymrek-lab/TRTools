.. overview_directive
.. |qcSTR overview| replace:: qcSTR generates plots that are useful for diagnosing issues in TR calling.
.. overview_directive_done

qcSTR
=====

|qcSTR overview|

Usage
-----
qcSTR takes as input a VCF file and outputs several plots in pdf format. To run qcSTR, use the following command::

    qcSTR \
  	--vcf <vcf file> \
   	--out test \
   	[additional options]


Required Parameters:

* :code:`--vcf <string>`: Input VCF file
* :code:`--out <string>`: Prefix to name output files


Optional Input Filters:

* :code:`--samples <string>`: File containing list of samples to include. If not specified, all samples are used.
* :code:`--period <int>`: Restrict to TRs with this motif length. e.g. to restrict to dinucleotide repeats, use :code:`--period 2`.
* :code:`--vcftype <string>`: Type of VCF files being merged. Default = :code:`auto`. Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.

Quality Plot Options:

* :code:`--quality`:  These options determine if the plot is stratified, and what 
  the x-axis represents. The y-axis always measures percentages,
  The x-axis is always cumulative decreasing. Multiple options can be specified
  by separating them with commas, no spaces, which will produce a different
  plot for each option e.g. :code:`per-locus,sample-strat`

  * :code:`per-locus`
    Compute the call quality at each locus averaged across all samples.
    Plot the distribution of those loci qualities.
    (This is the default for > 5 samples)
  * :code:`sample-strat` 
    Plot a separate line for each sample of the distribution of loci qualities
    for that sample.
    (This is the default for <= 5 samples)
  * :code:`per-sample`
    Compute the call quality for each sample averaged across all loci.
    Plot the distribution of those sample qualities.
  * :code:`locus-strat` 
    Plot a separate line for each locus of the distribution of sample qualities
    at that locus.
  * :code:`per-call`
    Plot the distribution of the quality of all calls.

* :code:`--quality-log-scale` 
  Make the quality plot x-axis logarithmic.


Outputs
-------

qcSTR outputs the following plots:

* outprefix+"xx-sample-callnum.pdf": a barplot giving the number of calls for each sample. Can be used to determine failed or outlier samples.
* outprefix+"xx-chrom-callnum.pdf": a barplot giving the number of calls for each chromosome. Can be useful to determine if the expected number of calls per chromosome are present.
* outprefix+"-diffref-histogram.pdf": a histogram of the difference from the reference allele (in number of repeat units) for each allele called. Can be used to visualize if there is a strong bias toward calling deletions vs. insertions compared to the reference, which might indicate a problem. The red line gives the cumulative fraction of TRs below each reference length.
* outprefix+"xx-diffref-bias.pdf": plots reference length (bp) vs. the mean difference in length of each allele called compared to the reference allele. It is expected that the mean difference should be around 0 for most settings. When this value starts to deviate from 0, e.g. for very long repeats, it could indicate a drop in call quality.

Example qcSTR command
---------------------

Example::

	FILE=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
	qcSTR \
  	--vcf ${FILE} \
  	--out test-qc

where :code:`$REPODIR` points to the root path of this repository.

