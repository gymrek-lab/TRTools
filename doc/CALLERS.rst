Supported TR Genotypers
=======================

TRTools currently supports 5 tandem repeat genotypers.
Here we summarize these genotypers and provide some basic specification of their functionality.
For more information on a genotyper, please see its website linked below.

+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
| Method (version tested) |  Supported TR classes    |  Num. TRs in reference  | Supported seq. tech.    |     Use case notes                   |     
+=========================+==========================+=========================+=========================+======================================+
|   AdVNTR_ (v1.3.3)      | Repeat units lengths     |   158,522 (genic hg19)\*|    Illumina/PacBio      | Infers allele lengths by default. May|
|                         | 6-100bp                  |                         |                         | alternatively identify putative      |
|                         |                          |                         |                         | frameshift mutations within VNTRs    |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
|ExpansionHunter_ (v3.2.2)| Designed for STRs (rep.) |   25 (hg19)\*           |    PCR-free Illumina    | Designed for inferring expanded.     |
|                         | unit <=6bp. Can handle   |                         |                         | repeats. Can handle repeats with     |
|                         | complex repeat structures|                         |                         | structures such as interruptions or  |
|                         | specified by regular     |                         |                         | regions containing multiple nearby   |
|                         | expressions              |                         |                         | repeats. Designed to run on a single |
|                         | e.g. (CAG)*(CCG)*        |                         |                         | sample.                              |     
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
| GangSTR_ (2.4.4)        | Designed for repeat unit |  829,233                | Paired-end Illumina     | Infers allele lengths only. Handles  |
|                         | lengths 1-20bp.          |                         |                         | both short TRs and TR expansions.    |
|                         |                          |                         |                         | Designed to run genome-wide on one or|
|                         |                          |                         |                         | more samples.                        |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
| HipSTR_ (v0.6.2)        | Repeat unit lengths 1-9bp| 1,620,030 (hg19)\*      | Illumina                | Infers both allele lengths and       |
|                         |                          |                         |                         | sequence. Can phase with SNPs. Only  |
|                         |                          |                         |                         | handles TRs spanned by short reads.  |
|                         |                          |                         |                         | Designed to run genome-wide on one or|
|                         |                          |                         |                         | more samples.                        |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
| PopSTR_ (v2.0)          | Repeat unit length 1-6bp.| 540,1401 (hg38)         | Illumina                | Infers allele lengths only. Handles  |
|                         |                          |                         |                         | both short TRs and TR expansions.    |
|                         |                          |                         |                         | Designed to run genome-wide on one or|
|                         |                          |                         |                         | more samples.                        |
+-------------------------+--------------------------+-------------------------+-------------------------+--------------------------------------+
\* Indicates that users may also supply custom reference panels of TR loci. Reference panel numbers shown above are based on downloads from the github repository of each tool as of July 2, 2020.

TRTools can be extended to support other genotypers that generate files meeting the VCF standard.
We welcome community contributions to help support them. If that interests you, please
see :ref:`Contributing` for more information.

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR

