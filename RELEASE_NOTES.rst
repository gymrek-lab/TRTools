4.0.0
-----

Features:

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  This makes both the underlying VCF reading code and the TRTools code
  significantly faster and more memory efficient. For instance, the loading of
  VCFs into memory is now > 15x faster for VCFs with many samples.
* DumpSTR has a new flag --zip to produce a bgzipped and tabix-indexed output VCF

Command line interface changes:

* StatSTR's --region option now requires the input VCF to be bgzipped and tabix indexed.
* If DumpSTR is used on an input VCF with unexpectedly typed
  INFO fields 'AC', 'REFAC', 'HET', 'HWEP', 'HRUN' or FORMAT field 'FILTER',
  it now errors out and asks you to rename those fields before rerunning 
  DumpSTR. (If they already exist but have the correct number and type DumpSTR
  will overwrite them and issue a warning in case that was not intended)

Output changes:

* DumpSTR call level filters now have the value of the call which triggered
  the filter appended to the filter name in the FILTER format field. (e.g.
  GangSTRCallMinDepth12 because the field had a depth of 12 and that's lower
  than the required min depth)
* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in the 
  samplog output file
* When running DumpSTR, loci where all the calls were either already nocalls
  or were filtered by call-level filters before the locus-level filters were run are now
  marked as 'NO_CALLS_REMAINING' instead of 'PASS'.
* When DumpSTR filters a call and replaces each of its format fields with the no call
  '.', fields with more than one value are now represented correctly. For example,
  for 2 values '.,.' is used rather than just a single '.'

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
* DumpSTR now correctly adds ##FILTER=<ID=PASS,Description="All filters passed">
  to the header line
* DumpSTR now no longer says HipSTRCallFlankIndels is applied to nocalls

Quality of life improvements:

* StatSTR, when printing output to a file, now prints timing diagnostics to stdout.
* DumpSTR will fail faster if output directory does not exist

Regressions:

* The --gangstr-require-support filter has been disabled.

Outstanding bugs:

* The dumpSTR ExpansionHunter ADFL ADIR ADSP filters have never worked
* DumpSTR remains untested on ExpansionHunter filters and files
* DumpSTR remains untested on loci with variable ploidy and/or partially
  genotyped samples (e.g. .|2)

