Release 4.0.0

Features:
* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  This makes both the underlying VCF reading code and the TRTools code
  significantly faster and memory efficient. For instance, the loading of 
  VCFs into memory is now > 15x faster for VCFs with many samples.
* DumpSTR has a new flag --zip to produce a bgzipped and tabix-indexed output VCF

Command Line Interface changes:
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
* DumpSTR now adds ##FILTER=<ID=PASS,Description="All filters passed">
  to the header line
* When DumpSTR filters a call and replaces each of its format fields with the no call
  '.', fields with number>1 are now replaced with '.,.' e.g. for a number=2 field, rather
  than just a single '.'

Bug fixes:
* The AC, REFAC fields that DumpSTR output used to be incorrect, are now correct
* If you specify --drop-filtered DumpSTR will no longer set all values in the 
  output .loclog.tab file to zero and instead set them to their proper values
  (which are the same as if you had not specified --drop-filtered)

Quality of life improvements:
* DumpSTR will fail faster if output directory does not exist

Regressions:
* The --gangstr-require-support filter has been disabled.

Outstanding bugs:
* The dumpSTR ExpansionHunter ADFL ADIR ADSP filters have never worked
* DumpSTR remains untested on ExpansionHunter filters and files
* DumpSTR remains untested on loci with variable ploidy and/or partially
  genotyped samples (e.g. .|2)
