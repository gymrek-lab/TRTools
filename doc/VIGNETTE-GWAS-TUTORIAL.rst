Performing a TR-based GWAS using imputed genotypes
==================================================

Overview
--------

* TRTools tools used: annotaTR
* Other tools used: Beagle, bcftools, plink2

The tutorial walks through the following steps:

* :ref:`Step 1: Imputing TRs from genotype data <step1>`

  * :ref:`Step 1.1: Genotype preprocessing <step1_1>`
  * :ref:`Step 1.2: Imputing TRs in each batch <step1_2>`
  * :ref:`Step 1.3: Extracting TRs from each batch <step1_3>`
  * :ref:`Step 1.4: Merging results from each batch <step1_4>`

* :ref:`Step 2: Computing dosages with annotaTR <step2>`
* :ref:`Step 3: Running GWAS with TR genotypes <step3>`

Background
----------

To run a GWAS, you will need to have genotypes and phenotype data for your phenotype of interest available for a large cohort. Whereas SNP-based GWAS typically tests each variant for linear association between the number of alternate alleles an individual has vs. the phenotype value, the TR-based GWAS below tests for linear association between "TR dosages" and the phenotype at each TR. There are multiple ways to define TR dosage defined below, but all are related to the mean number of copies of the repeat each person has across their two alleles at a particular TR. This overall framework is similar to that described in `Margoliash et al. 2023 <https://pubmed.ncbi.nlm.nih.gov/38116119/>`_.

.. image:: source/images/tr-vs-snp-gwas.png
   :width: 600

This tutorial shows how to perform a GWAS using TR genotypes obtained from imputation. Alternatively, genotypes can be obtained directly from WGS data using one of the supported callers_. If you already have quality filtered TR calls available and do not need to perform imputation, you can jump directly to step 2 below.

A WDL workflow performing steps 1 and 2 is available on our `CAST workflows <https://github.com/CAST-genomics/cast-workflows/blob/main/tr-imputation/wdl/batch_imputation.wdl>`_ page.

.. _step1:

Step 1: Imputing TRs from genotype data
---------------------------------------

Imputation overview
~~~~~~~~~~~~~~~~~~~

The steps below assume you will be imputing, rather than directly genotyping, TRs for your cohort. There are multiple scenarios when imputing TRs might be preferred over direct genotyping:

* If you only have genotype data (typically SNPs and indels) available for your GWAS cohort and therefore cannot directly genotype TRs
* Even when WGS data is available, it can be expensive and time-consuming to perform genome-wide TR genotyping. For example, the pipeline below is based on what our group has applied to the `All of Us dataset <https://workbench.researchallofus.org/>`_ which at the time of writing has WGS for around 246,000 individuals. Running HipSTR on a single sample takes around 16 hours/sample, and would cost tens to hundreds of thousands of dollars at this scale. Instead, we can rely on imputation, which can give accurate genotypes at the majority of TRs and is much faster.

The steps below guide you through performing TR imputation on cohorts of up to hundreds of thousands of samples.

Prerequisites
~~~~~~~~~~~~~

You will need the following files and tools to run the imputation step:

* **Individual-level SNP/indel genotypes (VCF format)**. You will need SNP/indel genotypes for your GWAS cohort. Either phased or unphased genotypes is fine. Imputation with Beagle will require these to be in VCF format. More on preprocessing genotypes below. 

* **SNP-TR reference panel (VCF or BREF3 format)**. You will need a precomputed referene panel that contains phased SNP and TR genotypes from an orthogonal cohort. The panel needs to be in either VCF or BREF3 format to be compatible with Beagle. We have generated two such panels:

  * `Saini et al. <https://gymreklab.com/2018/03/05/snpstr_imputation.html>`_ This panel is in the GRCh37 reference build and contains 445,725 STRs and 27M SNPs/indels for 2,504 samples from the 1000 Genomes Project. Genotypes had been imputed into 1000 Genomes samples based on calls in the Simons Simplex Collection, which is predominantly European. It is also restricted to STRs called by HipSTR with repeat unit lengths <= 6bp. This panel is available in VCF format with one file/chromosome.
  * `Ziaei-Jam et al. <https://github.com/gymrek-lab/EnsembleTR/blob/fix-ref/README.md#version-iii-of-reference-snptr-haplotype-panel-for-imputation-of-tr-variants>`_ This panel is in the GRCh38 reference build and contains 1,070,698 TRs and 70M Snps/indels from 3,202 samples from the 1000 Genomes Project. TR genotypes are based on `EnsembleTR <https://github.com/gymrek-lab/ensembleTR>_` and contain both STRs (repeat unit 1-6bp) and VNTRs (repeat unit 7+bp). The steps below were specifically tested with this panel but should also be mostly relevant to imputation with the Saini reference panel. This panel is available in VCF and BREF3 formats with one file/chromosome.

* **Genetic map (optional)**. You can optionally provide a genetic map to Beagle with cM coordinates, which are more accurate for modeling recombination than bp coordinates. We used the `GRCh38 map files available from Beagle <https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/>`_.

* **Tools**

  * `Beagle <https://faculty.washington.edu/browning/beagle/beagle.27May24.118.jar>`_: Beagle is used for imputing TRs from SNPs/indels. The steps below were tested with :code:`beagle.27May24.118.jar` which the link points to.
  * `bcftools <https://samtools.github.io/bcftools/bcftools.html>`_: Bcftools is used for multiple steps including merging files and extracting TRs from the Beagle output.
  * `tabix <https://anaconda.org/bioconda/tabix>`_: Used for indexing VCF files


Imputation steps
~~~~~~~~~~~~~~~~

.. _step1_1:


Step 1.1: Genotype preprocessing
________________________________


The inputs to Step 1.2 below are one VCF file per chromosome per batch of samples.

Beagle requires genotypes of the target samples to be input in VCF format. If your files are in another format (e.g. PGEN or Plink BED) you will need to first convert them to VCF. Further, the SNP-TR reference panel files and map files are split by chromosome, so you will similarly want your input genotype files to be split by chromosome.

Finally, if you have a very large cohort (more than several thousand samples) imputation with Beagle can be very memory intensive. To avoid memory errors, we recommend splitting your genotypes into batches of 1,000 samples each. We have found the :code:`bcftools plugin split` command helpful for creating all the batches at once. An example command is below::

	bcftools plugin split full_vcf_chr1.vcf.gz -G sample_groups.txt -Oz -o .

    	for f in *.vcf.gz; do tabix -p vcf $f; done

In this command:

* :code:`full_vcf_chr21.vcf.gz` is the VCF file with genotypes for the whole cohort on chr21.
* :code:`sample_groups.txt` is a file defining how the samples should be split into batches. An example is below::

	sample1 -       batch1_chr21
	sample2 -       batch1_chr21
	sample3 -       batch1_chr21
	sample4 -       batch2_chr21
	sample5 -       batch2_chr21
	sample6 -       batch2_chr21

The :code:`split` command above will result in files :code:`batch1_chr21.vcf.gz` and :code:`batch2_chr21.vcf.gz`.

An full workflow for generating these VCF subsets in All of Us is available from our `CAST workflows <https://github.com/CAST-genomics/cast-workflows/tree/main/subset_vcf>`_ page. (Note in that workflow given the cohort size, we further broke up the split step by genomic region and then concatenate the results for all the regions from each chromosome in a final step).

.. _step1_2:


Step 1.2: Imputing TRs in each batch
____________________________________

TODO

.. _step1_3:


Step 1.3: Extracting TRs from each batch
________________________________________

TODO


.. _step1_4:


Step 1.4: Merging results from each batch
_________________________________________

TODO



.. _step2:

Step 2: Computing dosages with annotaTR
---------------------------------------

TODO

.. _step3:

Step 3: Running GWAS with TR genotypes
--------------------------------------

TODO

.. _callers: https://trtools.readthedocs.io/en/stable/CALLERS.html