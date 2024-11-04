Performing a TR-based GWAS using imputed genotypes
==================================================

Overview
--------

TRTools tools used: annotaTR

Other tools used: Beagle, bcftools, plink2

The tutorial walks through the following steps:

* :ref:`Step 1: Imputing TRs from genotype data <step1>`
* :ref:`Step 2: Computing dosages with annotaTR <step2>`
* :ref:`Step 3: Running GWAS with TR genotypes <step3>`

Background
----------

To run a GWAS, you will need to have genotypes and phenotype data for your phenotype of interest available for a large cohort. Whereas SNP-based GWAS typically tests each variant for linear association between the number of alternate alleles an individual has vs. the phenotype value, the TR-based GWAS below tests for linear association between "TR dosages" and the phenotype at each TR. There are multiple ways to define TR dosage defined below, but all are related to the mean number of copies of the repeat each person has across their two alleles at a particular TR. This overall framework is similar to that described in `Margoliash et al. 2023 <https://pubmed.ncbi.nlm.nih.gov/38116119/>`_.

.. image:: source/images/tr-vs-snp-gwas.png
   :width: 600

This tutorial shows how to perform a GWAS using TR genotypes obtained from imputation. Alternatively, genotypes can be obtained directly from WGS data using one of the supported callers_. If you already have quality filtered TR calls available and do not need to perform imputation, you can jump directly to step 2 below.

.. _step1:

Step 1: Imputing TRs from genotype data
---------------------------------------

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