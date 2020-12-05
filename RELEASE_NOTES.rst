Release 3.0.0

Features:
* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  Reading the underlying VCFs should be > 15x faster for VCFs with many samples,
  and all the trtools code should also be significantly faster.
* DumpSTR has a new flag --zip to produce a bgzipped and tabix-indexed output VCF

Command Line Interface changes:
* If DumpSTR is used on an input VCF with unexpectedly typed
  INFO fields 'AC', 'REFAC', 'HET', 'HWEP', 'HRUN' or FORMAT field 'FILTER',
  it now errors out and asks you to rename those fields before rerunning 
  DumpSTR. (If they already exist but have the correct number and type DumpSTR
  will overwrite them and issue a warning in case that was not intended)

Output changes:
* DumpSTR call level filters now have the value of the call which triggered
  the filter appended to the filter name in the FILTER format field.
* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in the 
  samplog output file
* DumpSTR now adds ##FILTER=<ID=PASS,Description="All filters passed">
  to the header line (TODO even when locus filters are not being run?)
* When DumpSTR filters a call and replaces the format fields with '.', fields
  with number>1 are now replaced with '.,.' e.g. for a number=2 field, rather
  than just a single '.'

Bug fixes:
* The AC, REFAC fields that DumpSTR output used to be incorrect, are now correct
* If you specify --drop-filtered DumpSTR will no longer set all the values in the 
  output .loclog.tab file to zero and instead set them to their proper values
  (which are the same as if you had not specified --drop-filtered)

Quality of life improvements:
* DumpSTR will fail faster if output directory does not exist

Qs:
* DumpSTR filters? What do I do with half called samples?
* Currently --gangstr-require--support is disabled
* Need to undo
    not putting numeric reasons at end of filter lines
    not updating format FILTER description
    update region filter description
* Need to update tests with current dumpSTR output
* Need to add bcftools annotate to bad preexisting fields help
