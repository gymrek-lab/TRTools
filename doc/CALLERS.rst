Supported TR Genotypers
=======================

TRTools currently supports 5 tandem repeat genotypers. It also supports the Beagle imputation software (see :ref:`below <Beagle_section>`).
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

TRTools can be extended to support other genotypers that generate VCF files.
We welcome community contributions to help support them. If that interests you, please
see :ref:`Contributing` for more information.

..
    please ensure this list of links remains the same as the one in the main README

.. _AdVNTR: https://advntr.readthedocs.io/en/latest/
.. _ExpansionHunter: https://github.com/Illumina/ExpansionHunter
.. _GangSTR: https://github.com/gymreklab/gangstr
.. _HipSTR: https://hipstr-tool.github.io/HipSTR/
.. _PopSTR: https://github.com/DecodeGenetics/popSTR

.. _Beagle_section:

Beagle
------

The Beagle_ software can take genotypes called by a TR genotyper in a set of reference samples and impute them into other samples that do not have directly genotyped TRs.
TRTools supports TR genotypes produced by any of the above genotypers and then imputed into other samples with Beagle except for PopSTR genotypes. For each tool
in this tool suite, unless it's docs specifically say otherwise, that tool can be used on Beagle VCFs as if those VCFs were produced directly by the underlying TR genotyper,
with no additional flags or arguments needed, as long as the steps below were followed to make sure the Beagle VCF is properly formatted.

Caveats:

* Beagle provides phased best-guess genotypes for each imputed sample at each TR locus. When run with the :code:`ap` or :code:`gp` flags Beagle will also output
  probabilities for each possible haplotype/genotype, respectively. These probabilities are also called dosages. While dosages are often more informative for downstream
  analyses than the best-guess genotypes located in the :code:`GT` format field (for instance, for association testing), TRTools currently does *not* support dosage
  based analyses and instead will only look at the :code:`GT` field. Feel free to submit PRs with features that handle dosages (see the :ref:`Contributing` docs).
* At each locus Beagle returns the most probable phased genotype. This will often but not always correspond to the most probable unphased genotype. For instance,
  it is possible that :code:`P(A|A) > P(A|B)` and :code:`P(A|A) > P(B|A)`, but :code:`P(A/A) = P(A|A) < P(A|B) + P(B|A) = P(A/B)`. Similarly, it is possible that
  :code:`P(A|B) > P(C|D)` and :code:`P(A|B) > P(D|C)`, but :code:`P(A/B) = P(A|B) + P(B|A) < P(C|D) + P(D|C) = P(C/D)`.
  TRTools currently does not take this into
  account and just uses the phased genotypes returned by Beagle. If you deem this to be an issue, feel free to submit PRs to help TRTools take this into account
  (see the :ref:`Contributing` docs).
* For callers which return sequences, not just lengths (e.g. HipSTR), if there are loci with multiple plausible sequences of the same length, then its possible
  that the most probable genotype returned by Beagle does not have the most probable length. For example, the following could be true of a single haplotype:
  :code:`Len(S_1) = L_1, Len(S_2) = L_1, Len(S_3) = L_2` and :code:`P(S_1) < P(S_3), P(S_2) < P(S_3)` but :code:`P(S_3) < P(S_1) + P(S_2)`.

An overview of steps to perform before Beagle imputation:

* The samples being imputed into must have directly genotyped loci that are also genotyped in the reference samples. This allows those samples to be 'matched' with samples in the reference.
* The genotypes of both the reference samples and samples of interest must be phased. That can be done by statistically phasing the genotypes prior to running Beagle imputation.
* The referece samples must also not contain any missing genotypes. Possible methods for dealing with that include removing loci with missing genotypes or using imputation to impute
  the missing genotypes prior to imputing the TRs.

The VCFs that Beagle outputs need to be preprocessed before use by TRTools. We have provided a tool :code:`trtools_prep_beagle_vcf.sh` to run on those VCFs.
After running this script, the files should be usable by any of the tools in TRTools.

In case of error, it may be useful to know what steps the script attempts to perform:

* It copies over source and command meta header lines from the reference panel to the imputed VCF so
  that it is clear which genotyper's syntax is being used to represent the STRs in the VCF.
* It copies over contig and ALT lines which is required for downstream tools including mergeSTR and is good practice to include in the VCF header.
* It annotates each STR with the necessary INFO fields from the reference panel that Beagle dropped from the imputed VCF.
* The imputed VCF contains both TR loci and the shared loci (commonly SNPs) that were used for the imputation.
  This script removes the non-STR loci (identified as those loci not having STR-specific INFO fields).

.. _Beagle: http://faculty.washington.edu/browning/beagle/beagle.html
.. _bcftools: https://samtools.github.io/bcftools/bcftools.html
