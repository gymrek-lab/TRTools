The TRHarmonizer Library
========================

The TRHarmonizer Python library provides a uniform interface for accessing VCFs created by 
different tandem repeat (TR) genotypers. This library is the shared basis for all the 
command-line tools in the TRTools package. It is designed to cleanly handle differences in 
how different genotypers represent alleles, quality-scores and other metadata describing TR 
genotypes. This allows coding against a uniform interface while analyzing genetic variation 
at TRs regardless of which genotyper was used. The TRHarmonizer library also allows third 
parties to leverage the harmonization functionality outside of the command-line tools provided 
in TRTools.

A major challenge in analyzing TR genotypes is that alleles are represented differently in VCF 
outputs of different genotypers. The example below for chr21:47251618 (hg19) genotyped in 
Platinum Genomes sample NA12878 shows the different ways reference and alternate alleles are 
specified in VCFs by the genotypers which TRTools currently supports.

* adVNTR*, GangSTR, HipSTR:

  * REF: AGTTAGTTAGTTAGTT
  * ALT: AGTTAGTTAGTTAGTTAGTT

* ExpansionHunter:

  * REF: A
  * ALT: <STR5>
  * INFO: REF=4;RU=AGTT

* PopSTR:

  * REF: AGTTAGTTAGTTAGTT
  * ALT: <5>
  * INFO: Motif=AGTT

\* Note that while this TR was not called by AdVNTR because its motif is too short, 
AdVNTR output represents alleles in the same format as HipSTR and GangSTR.

Furthermore, consider the example at chr21:16402147

* adVNTR, GangSTR, HipSTR:

  * REF: AAATAAATAAATAAATAAAT
  * ALT: AAATAAATAAATAAAT

* ExpansionHunter:

  * REF: A
  * ALT: <STR4>
  * INFO: REF=5;RU=AAAT

* PopSTR:

  * REF: AAATAAATAAATAAATAAATAATAAA
  * ALT: <5.5>
  * INFO: Motif=AAAT

Here we see that popSTR’s representation of alleles changes to specify impurities and partial 
repeats.

The key function of the TRHarmonizer module, HarmonizeRecord, takes as input a cyvcf2 [1]_ 
record (a cyvcf2.Variant object) and a VCF type (one of:  “advntr”, “eh”, “gangstr”, 
“hipstr” or “popstr”, corresponding to the supported genotypers) and outputs a TRRecord 
object storing alleles and other metadata in a 
standardized format. This allows downstream analyses to proceed agnostic of the genotyper 
which created the record. The TRRecord stores allele length genotypes as the number of 
copies of the motif corresponding to that length. This number is a float to allow for 
impurities and partial repeats. For genotypers which infer sequence alleles, the record 
additionally stores the sequence of the allele in all uppercase. In addition to alleles,
a TRRecord also provides a uniform method for accessing the TR motif, per-sample quality 
scores and other metadata supplied by the underlying genotyper.

TRHarmonizer is designed to be lightweight, and as such there are similar yet more complex 
use-cases that TRHarmonizer intentionally does not support. It does not have any insight 
into sequencing technologies which produce data that is later processed by TR genotypers 
into VCFs. As such it relies on the alleles, calls and associated quality scores output 
by the genotypers, each of which use their own models to compute quality scores. TRTools 
makes no attempt to modify those scores based on sequencing errors or other sources of error. 

TRHarmonizer also does not handle differences in variant coordinates, whether due to differences 
in choice of variant reference set or differences between calling algorithms. Note that this is 
only relevant to compareSTR, as that is the only one of our tools designed to process TRs from 
multiple VCFs produced by different genotypers simultaneously. The types of differences related 
to variant coordinates that TRHarmonizer does not handle includes:

* Repeat regions which some callers choose to represent as a single variant and other callers 
  represent as multiple variants
* Overlapping variants of different lengths due to decisions about whether to phase the 
  repeat variant with other nearby variants
* Overlapping variants of different lengths due to different choices as to which parts of a 
  locus constitute impure repeats and which constitute flanking regions

Rather, TRHarmonizer restricts itself to comparing variants called by different callers 
whose reference alleles start and end at the same base pairs. Handling different variant 
representations is a complex problem that has been the subject of significant work [2]_ [3]_ and 
is best handled by haplotype comparison tools which have been tailored to the specific 
use-case at hand. If there are compelling use-cases specific to TRs that are not handled by 
prior work then TRHarmonizer could be extended to tackle them, but that is currently out of scope.

Finally, TRHarmonizer can be readily extended to support any TR genotyping tool built on top of 
any sequencing or genotyping technology as long as the tool produces a valid VCF file representing 
each TR as a distinct record in the VCF. Supporting additional tools simply requires adding a
short function to the TRHarmonizer module converting records to the standardized format described 
above.

.. [1] Brent S Pedersen, Aaron R Quinlan, cyvcf2: fast, flexible variant analysis with Python, Bioinformatics, Volume 33, Issue 12, 15 June 2017, Pages 1867-1869, https://doi.org/10.1093/bioinformatics/btx057 
.. [2] Cleary, John G., et al. "Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines." BioRxiv (2015): 023754.
.. [3] https://github.com/Illumina/hap.py

