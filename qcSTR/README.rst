qcSTR
=====

qcSTR is a tool for generating various plots that are useful for diagnosing issues in TR calling.

Usage
-----
qcSTR takes as input a VCF file and outputs several plots in pdf format. To run qcSTR, use the following command:

.. code-block::

	qcSTR \
  	--vcf <vcf file> \
   	--out test \
   	[additional options]


Required Parameters:

* :code:`--vcf <string>`: Input VCF file
* :code:`--out <string>`: Prefix to name output files

Recommended Parameters:

* :code:`--vcftype <string>`: Type of VCF files being merged. Default=:code:'auto'. Must be one of: :code:'gangstr', :code:'advntr', :code:'hipstr', :code:'eh', :code:'popstr'.

Optional Parameters:

* :code:`--samples <string>`: File containing list of samples to include. If not specified, all samples are used.
* :code:`--period <int>`: Restrict to TRs with this motif length. e.g. to restrict to dinucleotide repeats, use :code:`--period 2`.

Outputs
-------

qcSTR outputs the following plots:

* outprefix+"xx-sample-callnum.pdf": a barplot giving the number of calls for each sample. Can be used to determine failed or outlier samples.
* outprefix+"xx-chrom-callnum.pdf": a barplot giving the number of calls for each chromosome. Can be useful to determine if the expected number of calls per chromosome are present.
* outprefix+"-diffref-histogram.pdf": a histogram of the difference from the reference allele (in number of repeat units) for each allele called. Can be used to visualize if there is a strong bias toward calling deletions vs. insertions compared to the reference, which might indicate a problem. The red line gives the cumulative fraction of TRs below each reference length.
* outprefix+"xx-diffref-bias.pdf": plots reference length (bp) vs. the mean difference in length of each allele called compared to the reference allele. It is expected that the mean difference should be around 0 for most settings. When this value starts to deviate from 0, e.g. for very long repeats, it could indicate a drop in call quality.

Example qcSTR command
---------------------

.. code-block::

	FILE=${REPODIR}/test/common/sample_vcfs/compareSTR_vcfs/compare_vcf1.vcf.gz
	qcSTR \
  	--vcf ${FILE1} \
  	--out test-qc

where :code:`$REPODIR` points to the root path of this repository.

