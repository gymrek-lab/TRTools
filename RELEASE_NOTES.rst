5.1.0
-----

New features:

* Added prancSTR for mosaicism detection
* Added simTR for simulating NGS reads with stutter errors at TRs

5.0.2
-----

Bug fixes:

* MergeSTR now will no longer sometimes emit an alternate allele identical to the ref allele when
  dealing with flanking base pairs.

5.0.1
-----

Bug fixes:

* Remove outdated call in qcSTR to np.float()

5.0.0
-----

New features:

* associaTR has been released!

Current limitations:

* Does not support binary phenotypes yet
* Does not support producing data for plotting individual loci yet
* Values in the output file aside from the p-values, coefficients and 
  standard errors have not been fully tested

4.2.1
-----

Bug fixes:

* Fix bioconda build

4.2.0
-----

New features:

* TRTools can now read VCFs produced by Beagle imputation.

Bug fixes:

* MergeSTR now successfully merges files containing multiple chromosomes instead of emitting
  a 'stuck in infinite loop' message and crashing
* MergeSTR no longer crashes if run with the --verbose flag and the last position in each of the
  VCFs being merged isn't identical.
* StatSTR now errors out if any of the files listed in --samples contain no samples that are present
  in the input VCF.
* DumpSTR now reads call depth from the LC format field when the DP format field is not present.
  This was intended previously but was not happening.

Doc changes:

* Clarified in PUBLISHING.rst how to handle dependencies and how to publish to bioconda.
* requirements.txt was unneeded, so delete it and remove the reference to it from PUBLISHING.rst

4.1.0
-----

Functionality Changes:

* MergeSTR: Flanking basepairs are now removed from HipSTR records before merging.
  In particular, records with different flank lengths but the same repeat section will now merge instead of being
  marked as incompatible.

* CompareSTR: the tool now only compares records that start and end at the same position. If a partial overlap in records
  is detected, the program will output a warning to the user. This warning contains IDs of the records and their positions.

Misc:

* mergeutils: function GetMinHarmonizedRecords was transformed into GetIncrementAndComparability, which allows the caller
  to define custom predicate that decides whether records are comparable.

4.0.2
-----

Bug fixes:

* https://github.com/gymrek-lab/TRTools/issues/146 fixed record positions being compared twice
* CompareSTR: Decision on which records are comparable is now based on data from harmonized TRRecords,
  and not from the records directly from VCF readers. Thanks to this, HipSTR records which have different starting positions,
  but position of their repeat is at the same position are compared correctly (harmonization step removes this difference).
* MergeSTR failed on mixed ploidy samples (i.e. chrX). Fix one such bug. Note: none of the tools are 
  fully tested for chrX even with this fix.


4.0.1
-----

Bug fixes:

* https://github.com/gymrek-lab/TRTools/issues/143 Fix HipstrMinSuppReads filter when
  there are called samples but none have ALLREADS

4.0.0
-----

Features:

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF parsing.
  This makes both the underlying VCF reading code and the TRTools code
  significantly faster and more memory efficient. For instance, the loading of
  VCFs into memory is now > 15x faster for VCFs with many samples.
  Some tools will still need further updates to be usable for large datasets,
  but those updates should now be possible and much easier.
  (e.g. emitting progress reports to stdout as needed, flags to disable
  computations that cannot be done at such scale)
* DumpSTR has a new flag --zip to produce a bgzipped and tabix-indexed output VCF
* StatSTR now can calculate the entropy of the allele distribution at each locus with the
  --entropy flag
* The `TRTools documentation website <https://trtools.readthedocs.io/en/latest/>`_ now
  displays the release notes.

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
* CompareSTR's docs now explicitly tell the user to order phased calls to
  prevent spurious mismatching. If phasing is not desired, use --ignore-phasing
* CompareSTR will now error if at a single locus both files do not have either all
  unphased calls, or all phased calls. If phasing is not desired, use --ignore-phasing

Output changes:

* DumpSTR call level filters now have the value of the filter and the value 
  which triggered the filter appended to the filter name in the FILTER format field.
  (e.g. GangSTRCallMinDepth20_12 because the field had a depth of 12 and that's lower
  than the required min depth of 20)
* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in the 
  samplog output file
* When running DumpSTR, loci where all the calls were either already nocalls
  or were filtered by call-level filters before the locus-level filters were run are now
  marked as 'NO_CALLS_REMAINING' instead of 'PASS'.
* When DumpSTR filters a call and replaces each of its format fields with the no call
  '.', fields with more than one value are now represented correctly. For example,
  for 2 values '.,.' is used rather than just a single '.'
* MergeSTR header lines are now copied over from the input VCFs instead of
  only copying over a few recognized fields (e.g. ID and Length
  were the only contig fields that were previously retained, but URL wouldn't be)
* MergeSTR output alt alleles for eh and popstr are now ordered by length.
  MergeSTR output alt alleles for advntr, gangstr and hipstr, when there are multiple
  alt alleles of the same length, are now ordered alphabetically instead
  of arbitrarily.
* CompareSTR no longer outputs the file <prefix>-callcompare.tab - the existence
  of that file was never documented, and besides, all its information could
  be seen more easily simply by looking at the input VCFs
* In CompareSTR's overall.tab file, the ranges in the format columns are now written
  [a,b) or [a,b] instead of a-b
* CompareSTR's locuscompare.tab file now outputs loci in the order they were
  encountered in the input VCfs as opposed to an arbitrary order
* The 'sample' column in CompareSTR's locuscompare.tab file has been renamed to
  'numcalls' to match the other two tab files.

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
* MergeSTR now outputs the same phase as the input files instead of always outputting
  unphased data
* MergeSTR now correctly outputs Number=A, G or R (number of entries in this field equal
  to number of alternate alleles at this locus, the number of alleles including the ref,
  or the number of unique polyploid genotypes) correctly in INFO and FORMAT fields instead
  of outputing Number=-1, -2 or -3
* CompareSTR claimed it was outputting the square (Pearson) correlation coefficient
  but was actually outputting the raw (unsquared) correlation coefficient. It is now
  outputting the squared coefficient as documented.
* CompareSTR now correctly compares unphased calls without regard to order in the VCF
  (e.g. 'AAAA/AAA' now matches against 'AAA/AAAA')
* CompareSTR's docs claimed the bubble plots axes were measured in basepair difference
  from the reference, but they were actually measured in number of repeats different
  from the reference. The behavior has not been changed and the claim has been updated
  to match the behavior.
* When using binned format fields in CompareSTR where the range of values did not
  evenly divide into the requested binsize, the highest valued bin used to always
  be the same size as all the other bins and include values over the
  limit specified by the user. Now it caps at that maximum.
  E.g. binsizes 0:210:50 used to create the bins
  [0,50), [50,100), [100,150), [150, 200), [200, 250)
  and now create the bins
  [0,50), [50,100), [100,150), [150, 200), [200, 210]
* When using binned format fields in CompareSTR where the range of values
  evenly divided into the requested binsize, loci which obtained the requested
  maximum would be excluded. They are now included.
  E.g. binsizes 0:200:50 used to create the bins
  [0,50), [50,100), [100,150), [150, 200) and samples with value 200 would
  not fall into any bin. This now creates the bins
  [0,50), [50,100), [100,150), [150, 200] and samples with value 200 fall into
  the last bin

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

3.0.3
-----

Bug fixes:

* Fixed a spot where qcSTR would crash because we passed Pandas a set instead of a list
* MergeSTR now writes out the header for the GT FORMAT field
