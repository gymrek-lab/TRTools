4.0.0
-----

Features:

* StatSTR now can calculate the entropy of the allele distribution at each locus with the
  --entropy flag


3.0.3
-----

Bug fixes:

* Fixed a spot where qcSTR would crash because we passed Pandas a set instead of a list
* MergeSTR now writes out the header for the GT FORMAT field
