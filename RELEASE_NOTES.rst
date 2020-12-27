4.0.0
-----

Features:

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  This makes both the underlying VCF reading code and the TRTools code
  significantly faster and memory efficient. For instance, the loading of 
  VCFs into memory is now > 15x faster for VCFs with many samples.
* DumpSTR has a new flag --zip to produce a bgzipped and tabix-indexed output VCF

Command line interface changes:

* StatSTR's --region option now requires the input VCF to be bgzipped and tabix indexed.
* If DumpSTR is used on an input VCF with unexpectedly typed
  INFO fields 'AC', 'REFAC', 'HET', 'HWEP', 'HRUN' or FORMAT field 'FILTER',
  it now errors out and asks you to rename those fields before rerunning 
  DumpSTR. (If they already exist but have the correct number and type DumpSTR
  will overwrite them and issue a warning in case that was not intended)
* CompareSTR's docs used to claim that when comparing alleles from different callers
  those callers must use the same allele notation (e.g. implying that ExpansionHunter's
  '<STR12>' and GangSTR's 'ACACACACACAC' notation would always mismatch). That statement
  was never true for length based comparisons - CompareSTR has always been able to
  do length based comparisons regardless of notation. The incorrect claim has been
  removed from CompareSTR's docs.
* CompareSTR's docs now explicitly state how to order phased and unphased calls to
  prevent calls from spuriously mismatching

Output changes:

* DumpSTR call level filters now have the value of the call which triggered
  the filter appended to the filter name in the FILTER format field. (e.g.
  GangSTRCallMinDepth12 because the field had a depth of 12 and that's lower
  than the required min depth)
* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in the 
  samplog output file
* DumpSTR now adds ##FILTER=<ID=PASS,Description="All filters passed">
  to the header line
* When DumpSTR filters a call and replaces each of its format fields with the no call
  '.', fields with number>1 are now replaced with '.,.' e.g. for a number=2 field, rather
  than just a single '.'
* MergeSTR header lines are now copied over from the input VCFs instead of
  only copying over a few recognized fields (e.g. ID and Length
  were the only contig fields that were previously retained, but URL wouldn't be)
* MergeSTR output alt alleles for eh and popstr are now ordered by length.
  MergeSTR output alt alleles for advntr, gangstr and hipstr, when there are multiple
  alt alleles of the same length, those alleles are now ordered alphabetically instead
  of arbitrarily.
* CompareSTR no longer outputs the file <prefix>-callcompare.tab - the existence
  of that file was never documented, and besides, all its information could
  be seen more easily simply by looking at the input VCFs

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

Bug fixes:

* The AC, REFAC fields that DumpSTR output used to be incorrect, are now correct
* If you specify --drop-filtered DumpSTR will no longer set all values in the 
  output .loclog.tab file to zero and instead set them to their proper values
  (which are the same as if you had not specified --drop-filtered)
* MergeSTR now outputs the same phase as the input files instead of always outputting
  unphased data
* MergeSTR now correctly outputs Number=A, G or R correctly in FORMAT fields instead
  of outputing Number=-1, -2 or -3

Quality of life improvements:

* StatSTR, when printing output to a file, now prints timing diagnostics to stdout.
* DumpSTR will fail faster if output directory does not exist
* When encountering issues with identifying the caller type for each input VCF,
  MergeSTR now prints an error and gracefully returns instead of dying to
  an uncaught exception
* MergeSTR incompatible INFO field warnings now specify which locus has an
  incompatible field

Regressions:

* The --gangstr-require-support filter has been disabled.

Outstanding bugs:

* The dumpSTR ExpansionHunter ADFL ADIR ADSP filters have never worked
* DumpSTR remains untested on ExpansionHunter filters and files
* DumpSTR remains untested on loci with variable ploidy and/or partially
  genotyped samples (e.g. .|2)
* When running CompareSTR with the --stratify options where --stratify-file
  is either not specified or is explicitly set to zero, for each format field
  all calls where the value of that field in vcf1 does not fall into the same
  bin as the value of that field in vcf2 are silently not compared for that format field.
  The correct behavior here is probably to create paired bins based on a range
  of values from vcf1 and a range from vcf2. Regardless, the behavior here should
  be documented.
