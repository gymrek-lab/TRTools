Release 3.0.0

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  Should be > 15x faster for VCFs with many samples

* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in the 
  samplog output file
* DumpSTR call level filters now have the value of the call which triggered
  the filter appended to the filter name in the FILTER format field.
* When DumpSTR filters a call and replaces the format fields with '.', fields
  with number>1 are now replaced with '.,.' e.g. for a number=2 field, rather
  than just a single '.'
* DumpSTR now adds ##FILTER=<ID=PASS,Description="All filters passed">
  to the header line (TODO even when locus filters are not being run?)

* DumpSTR filters? What do I do with half called samples?
* Add write gzipped file option
