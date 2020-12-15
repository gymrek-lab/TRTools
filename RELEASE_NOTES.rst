4.0.0
-----

Features:

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  This makes both the underlying VCF reading code and the TRTools code
  significantly faster and more memory efficient. For instance, the loading of
  VCFs into memory is now > 15x faster for VCFs with many samples.

Command line interface changes:

* StatSTR's --region option now requires the input VCF to be bgzipped and tabix indexed.

Python interface changes:

* The trtools.utils.tr_harmonizer module has been reworked to use cyvcf2,
  and in doing so a large portion of its interface has changed in small ways.
* The big conceptual change is that instead of repeatedly calling a method
  on a TRRecord object like GetStringGenotype for each sample in the VCF,
  instead you call the new corresponding method GetStringGenotypes once,
  and it returns a numpy array of values where the first axis of the array 
  ranges over the samples.
* The way missing calls and samples with lower than maximal
  ploidy are handled is now tested and documented. These representations
  of these genotypes have been aligned with cyvcf2's standards.
  For more info, see the docs of the index, length and 
  string genotype getter methods.

Quality of life improvements:

* StatSTR, when printing output to a file, now prints timing diagnostics to stdout.
