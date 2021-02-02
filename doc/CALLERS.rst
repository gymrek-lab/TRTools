.. _callers:

Supported TR Genotypers
=======================

TRTools currently supports 5 tandem repeat genotypers.
We summarize them in the first table and provide some basic parameters of their functionality in the second.
For more information on a genotyper, please see its website linked below.

+----------------------------+--------------------------------------+
| Genotyper (version tested) |     Use case notes                   |
+============================+======================================+
|      AdVNTR_ (v1.3.3)      | Infers allele lengths. May           |
|                            | alternatively identify putative      |
|                            | frameshift mutations within          |
|                            | VNTRs (6+bp repeat units).           |
|                            | Designed for targeted genotyping of  |
|                            | VNTRs.                               |
|                            | May be run on large panels of        |
|                            | TRs but is compute-intenstive.       |
+----------------------------+--------------------------------------+
| ExpansionHunter_ (v3.2.2)  | Handles repeats with                 |
|                            | structures such as interruptions or  |
|                            | nearby repeats.                      |
|                            | Designed for targeted genotyping of  |
|                            | expansions at                        |
|                            | known pathogenic TRs but may be run  |
|                            | genome-wide on short and             |
|                            | expanded TRs using a custom TR panel.|
+----------------------------+--------------------------------------+
|    GangSTR_ (2.4.4)        | Designed for genome-wide genotyping  |
|                            | of short or expanded TRs.            |
+----------------------------+--------------------------------------+
|    HipSTR_ (v0.6.2)        | Designed for genome-wide genotyping  |
|                            | of STR (1-6bp repeat units) alleles  |
|                            | shorter than the read length.        |
|                            | Can phase repeats with SNPs.         |
+----------------------------+--------------------------------------+
|    PopSTR_ (v2.0)          | Designed for genome-wide genotyping  |
|                            | of short or expanded TRs.            |
+----------------------------+--------------------------------------+

|

+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+
| Genotyper (version tested) |  Repeat unit lengths     | Alleles longer than reads? | Allele type inferred   |  # TRs in reference      | Sequencing technology   | # Samples at a time    |
+============================+==========================+============================+========================+==========================+=========================+========================+
|      AdVNTR_ (v1.3.3)      |  6-100bp                 | No                         | Length, frameshifts    |   158,522 (genic hg19)   | Illumina, PacBio        | Single                 |
+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+
| ExpansionHunter_ (v3.2.2)  | 1-6bp. Can handle        | Yes                        | Length                 |   25 (hg19)              | PCR-free Illumina       | Single                 |
|                            | complex repeat structures|                            |                        |                          |                         |                        |
|                            | specified by regular     |                            |                        |                          |                         |                        |
|                            | expressions              |                            |                        |                          |                         |                        |
+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+
|    GangSTR_ (2.4.4)        | 1-20bp                   | Yes                        | Length                 |  829,233 (hg19)          | Paired-end Illumina     | Many                   |
+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+
|    HipSTR_ (v0.6.2)        | 1-9bp                    | No                         | Length, sequence       | 1,620,030 (hg19)         | Illumina                | Many                   |
+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+
|    PopSTR_ (v2.0)          | 1-6bp                    | Yes                        | Length                 | 540,1401 (hg38)          | Illumina                | Many                   |
+----------------------------+--------------------------+----------------------------+------------------------+--------------------------+-------------------------+------------------------+

Since each of these tools take as input a list of TRs to genotype, they could also be used on custom panels of TR loci.
Tool information and reference panel numbers shown above are based on downloads from the github repository of each tool as of July 2, 2020.

In addition, TRTools has preliminary support for TR genotypes imputed by Beagle_
from callsets generated by any of the above tools. See :ref:`our description of that workflow
<working_with_beagle>` for more information.

TRTools can be extended to support other genotypers that generate VCF files.
We welcome community contributions to help support them. If that interests you, please
see :ref:`Contributing` for more information.

.. toctree::

   BEAGLE

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _Beagle: https://faculty.washington.edu/browning/beagle/beagle.html
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR



