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
  Samples in the list that are not included in the input vcf or
  are misspelled are silently ignored.
* :code:`--period <int>`: Restrict to TRs with this motif length. e.g. to restrict to dinucleotide repeats, use :code:`--period 2`.
* :code:`--vcftype <string>`: Type of VCF files being merged. Default = :code:`auto`. Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.

If you wish to run qcSTR on a more complicated subset of the input VCF, we suggest you use
:code:`dumpSTR` or `bcftools view <http://samtools.github.io/bcftools/bcftools.html#view>`_ to
filter the input VCF first, and then run qcSTR on the vcf those commands
outputed.

Quality Plot Options:

* :code:`--quality`:  This option determines if the plot is stratified, and what 
  distribution the y-axis represents. The x-axis always is the quality score, from one to
  zero, and the y-axis is always an increasing CDF. This can be specified multiple
  times to produce multiple plots (e.g. :code:`--quality per-locus --quality
  per-sample`) - each produced plot will have the type appended to its name.

  * :code:`per-locus`
    Compute the call quality at each locus averaged across all samples.
    Plot the distribution of those loci qualities.
    (This is the default for > 5 samples)
  * :code:`sample-stratified` 
    Plot a separate line for each sample of the distribution of loci qualities
    for that sample.
    (This is the default for <= 5 samples)
    Note: If you specify this for a vcf with many samples,
    the code may slow to a halt and/or the plot may be cluttered.
  * :code:`per-sample`
    Compute the call quality for each sample averaged across all loci.
    Plot the distribution of those sample qualities.
  * :code:`locus-stratified` 
    Plot a separate line for each locus of the distribution of sample qualities
    at that locus.
    Note: If you specify this for a vcf with many loci,
    the code may slow to a halt and/or the plot may be cluttered.
  * :code:`per-call`
    Plot the distribution of the quality of all calls.

* :code:`--quality-ignore-no-call` - by default, (sample, locus) pairs which
  were not called are treated as calls with zero quality for the the quality plot.
  With this option enabled, instead they are ignored. This may cause the
  plotting to crash if it causes some samples/loci to have <= 1 valid call.


Outputs
-------

qcSTR outputs the following plots:

* :code:`<outprefix>-sample-callnum.pdf`: a barplot giving the number of calls for each sample. Can be used to determine failed or outlier samples.
* :code:`<outprefix>-chrom-callnum.pdf`: a barplot giving the number of calls for each chromosome. Can be useful to determine if the expected number of calls per chromosome are present.
* :code:`<outprefix>-diffref-histogram.pdf`: a histogram of the difference from the reference allele (in number of repeat units) for each allele called. Can be used to visualize if there is a strong bias toward calling deletions vs. insertions compared to the reference, which might indicate a problem. The red line gives the cumulative fraction of TRs below each reference length.
* :code:`<outprefix>-diffref-bias.pdf`: plots reference length (bp) vs. the mean difference in length of each allele called compared to the reference allele. It is expected that the mean difference should be around 0 for most settings. When this value starts to deviate from 0, e.g. for very long repeats, it could indicate a drop in call quality.
* :code:`<outprefix>-quality.pdf`: plots the distribution of the quality of
  calls for this vcf. Will not be produced for vcfs which do not have quality
  metrics. If you specify the type of quality plot you wish to see with
  the :code:`--quality` option, then instead you will get a file named 
  :code:`<outprefix>-quality-<type>.pdf` for each type of plot you requested.
  The following are example quality plots:


.. image:: images/quality-per-locus.png
.. image:: images/quality-sample-stratified.png
.. image:: images/quality-per-sample.png
.. image:: images/quality-locus-stratified.png
.. image:: images/quality-per-call.png

Example qcSTR command
---------------------

Example::

	FILE=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
	qcSTR \
  	--vcf ${FILE} \
  	--out test-qc

where :code:`$REPODIR` points to the root path of this repository. See "Additional Examples" below for additional example qcSTR commands for different supported TR genotypers based on example data files in this repository.

Wishlist
--------
A :code:`--quality-log-scale` option to expand the level of differentiation of qualities near 1.
A :code:`--quality-smooth` option to smooth the quality plots using :code:`sklearn.neighbors.KernelDensity(kernel='gaussian')`


Additional Examples
-------------------

Below are additional :code:`qcSTR` examples using VCFs from supported TR genotypers. Data files can be found in the :code:`example-files` directory of this repository::

  # GangSTR
  qcSTR --vcf ${REPODIR}/example-files/trio_chr21_gangstr.sorted.vcf.gz --out test_qc_gangstr --period 4 --quality per-locus

  # HipSTR
  qcSTR --vcf ${REPODIR}/example-files/trio_chr21_hipstr.sorted.vcf.gz --out test_qc_hipstr --vcftype hipstr --samples example-files/ex-samples.txt

  # ExpansionHunter
  qcSTR --vcf ${REPODIR}/example-files/NA12878_chr21_eh.sorted.vcf.gz --out test_qc_eh

  # AdVNTR
  qcSTR --vcf ${REPODIR}/example-files/NA12878_chr21_advntr.sorted.vcf.gz --out test_qc_advntr

  # PopSTR
  qcSTR --vcf ${REPODIR}/example-files/trio_chr21_popstr.sorted.vcf.gz --out test_qc_popstr

