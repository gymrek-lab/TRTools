# Changelog

## [6.2.0](https://github.com/gymrek-lab/TRTools/compare/v6.1.0...v6.2.0) (2026-02-19)


### Features

* add `--only-passing` option to statSTR ([#229](https://github.com/gymrek-lab/TRTools/issues/229)) ([51c0481](https://github.com/gymrek-lab/TRTools/commit/51c048191684342737edcc08e58f9d29806d563f))


### Documentation

* add note about bcftools version for merging to annotaTR ([#245](https://github.com/gymrek-lab/TRTools/issues/245)) ([6a2d4fb](https://github.com/gymrek-lab/TRTools/commit/6a2d4fb1291a87add573dd43f4e738bc883563ec))
* **annotaTR:** explain that SNPs aren't supported ([#248](https://github.com/gymrek-lab/TRTools/issues/248)) ([5132f34](https://github.com/gymrek-lab/TRTools/commit/5132f3499de752a880ca8c339b29856a97bc7436))
* **annotaTR:** fix beagpleap typo ([#260](https://github.com/gymrek-lab/TRTools/issues/260)) ([f8ef1e9](https://github.com/gymrek-lab/TRTools/commit/f8ef1e96d7d4cac187217543aca5562eba126623))
* drop support for py 3.8 and update to poetry 2 in our CI ([#249](https://github.com/gymrek-lab/TRTools/issues/249)) ([6bdc596](https://github.com/gymrek-lab/TRTools/commit/6bdc596db02dba1e38f7ab11fba8cb20377e0eaf))

## [6.1.0](https://github.com/gymrek-lab/TRTools/compare/v6.0.2...v6.1.0) (2024-11-06)


### Features

* Add longtr support ([#232](https://github.com/gymrek-lab/TRTools/issues/232)) ([93ff231](https://github.com/gymrek-lab/TRTools/commit/93ff231a74a8b0cd03ee957b8374abfc9451119d))
* Add options to facilitate debugging annotaTR runs on large files ([#234](https://github.com/gymrek-lab/TRTools/issues/234)) ([0d5d96b](https://github.com/gymrek-lab/TRTools/commit/0d5d96b8699d043aa88dc62866ded0465e0ba713))
* adding the annotaTR tool ([#226](https://github.com/gymrek-lab/TRTools/issues/226)) ([59dd532](https://github.com/gymrek-lab/TRTools/commit/59dd5323b6d58b38abe36836fbfe13f789ff4048))
* support for Apple M1 silicon ([#227](https://github.com/gymrek-lab/TRTools/issues/227)) ([dfec9c2](https://github.com/gymrek-lab/TRTools/commit/dfec9c229955f2862fb796d3d3360d1044f9cdc8))


### Bug Fixes

* associatr test files overwritten during pytest ([#230](https://github.com/gymrek-lab/TRTools/issues/230)) ([2b955be](https://github.com/gymrek-lab/TRTools/commit/2b955bed4d7d81ad5ba25a392fdd81df90a534d9))
* bug in beagle dosages with no alt alleles ([#231](https://github.com/gymrek-lab/TRTools/issues/231)) ([ae3ca15](https://github.com/gymrek-lab/TRTools/commit/ae3ca15e64617f8502cb56a71c839ef8a35f524e))
* enable filter_hrun with beagle file in dumpSTR ([be63e79](https://github.com/gymrek-lab/TRTools/commit/be63e79bb6402276e82cd7aea85bb76942217f75))


### Documentation

* generate internal links dynamically ([#228](https://github.com/gymrek-lab/TRTools/issues/228)) ([9402042](https://github.com/gymrek-lab/TRTools/commit/94020423b8667caeb20ea2615008d7441dfe5ac4))
* mention github codespaces and clarify dependency management ([#223](https://github.com/gymrek-lab/TRTools/issues/223)) ([45ccf28](https://github.com/gymrek-lab/TRTools/commit/45ccf28c306bcd98c6bc2b285d659571bf9afb6d))
* update biorxiv citation for prancSTR ([#235](https://github.com/gymrek-lab/TRTools/issues/235)) ([1e45d23](https://github.com/gymrek-lab/TRTools/commit/1e45d233afbf23e8d738786546f7c89ab6fe6805))
* warn about scipy deprecation in dumpSTR ([#221](https://github.com/gymrek-lab/TRTools/issues/221)) ([a1a3a39](https://github.com/gymrek-lab/TRTools/commit/a1a3a39b7b7fd2d2f694d3a41bfe7bcf72ff72c3))

## [6.0.2](https://github.com/gymrek-lab/TRTools/compare/v6.0.1...v6.0.2) (2024-06-24)


### Bug Fixes

* support for numpy 2.0 ([#220](https://github.com/gymrek-lab/TRTools/issues/220)) ([8661226](https://github.com/gymrek-lab/TRTools/commit/86612265fc0302ad86ebed48787530d0b4b95550))


### Documentation

* update hipstr link ([#217](https://github.com/gymrek-lab/TRTools/issues/217)) ([e1f8511](https://github.com/gymrek-lab/TRTools/commit/e1f8511ff070bc74b40318e95f7fe7a7bda57bb8))

## [6.0.1](https://github.com/gymrek-lab/TRTools/compare/v6.0.0...v6.0.1) (2024-02-21)


### Bug Fixes

* clone command in `test_trtools.sh` testing script ([#214](https://github.com/gymrek-lab/TRTools/issues/214)) ([756a516](https://github.com/gymrek-lab/TRTools/commit/756a5169a352c23d0681d70e46d78599a79017a2))

## [6.0.0](https://github.com/gymrek-lab/TRTools/compare/v5.1.1...v6.0.0) (2024-02-19)


### ⚠ BREAKING CHANGES

* remove pybedtools as a dependency. This slightly changes the interpretation of BED files. From now on, the end position of each region will not be included. ([#184](https://github.com/gymrek-lab/TRTools/issues/184))

### Features

* add version constraints, poetry, and `poetry.lock` file, and use `nox-poetry` and py3.12 in CI ([#202](https://github.com/gymrek-lab/TRTools/issues/202)) ([ed9c961](https://github.com/gymrek-lab/TRTools/commit/ed9c9610d0a6d0503efb07ad38a995a1957826ba))


### Bug Fixes

* deprecation of binom_test in scipy 1.12.0 and cryptic error occurring with incorrect VCFs in cyvcf2 0.30.26 ([#208](https://github.com/gymrek-lab/TRTools/issues/208)) ([78ec86a](https://github.com/gymrek-lab/TRTools/commit/78ec86a92dd7bdd8fffc1a88e416681af1b54234))


### Documentation

* reviewer and maintainer responsibilities with automated releases and PRs ([7762d87](https://github.com/gymrek-lab/TRTools/commit/7762d8796144597c4bdb5f68f51436bf14420ef2))
* specify conda install via biocontainers ([2807421](https://github.com/gymrek-lab/TRTools/commit/2807421568ad0a477a42346017820de808d2e5dd))


### Code Refactoring

* remove pybedtools as a dependency. This slightly changes the interpretation of BED files. From now on, the end position of each region will not be included. ([#184](https://github.com/gymrek-lab/TRTools/issues/184)) ([89c1dd3](https://github.com/gymrek-lab/TRTools/commit/89c1dd315a1cb96edb2d3479f7eb36990044a01c))

## [5.1.1](https://github.com/gymrek-lab/TRTools/compare/v5.1.0...v5.1.1) (2023-11-27)

### Bug fixes

* Remove stray files from source distribution

## [5.1.0](https://github.com/gymrek-lab/TRTools/compare/v5.0.2...v5.1.0) (2023-11-22)

### Features

* Added prancSTR for mosaicism detection
* Added simTR for simulating NGS reads with stutter errors at TRs

## 5.0.2

### Bug fixes

* MergeSTR now will no longer sometimes emit an alternate allele
    identical to the ref allele when dealing with flanking base pairs.

## 5.0.1

### Bug fixes

* Remove outdated call in qcSTR to `np.float()`

## 5.0.0

### Features

* associaTR has been released!

Current limitations:

* Does not support binary phenotypes yet
* Does not support producing data for plotting individual loci yet
* Values in the output file aside from the p-values, coefficients and
    standard errors have not been fully tested

## 4.2.1

### Bug fixes

* Fix bioconda build

## 4.2.0

### Features

* TRTools can now read VCFs produced by Beagle imputation.

### Bug fixes

* MergeSTR now successfully merges files containing multiple
    chromosomes instead of emitting a \'stuck in infinite loop\' message
    and crashing
* MergeSTR no longer crashes if run with the \--verbose flag and the
    last position in each of the VCFs being merged isn\'t identical.
* StatSTR now errors out if any of the files listed in \--samples
    contain no samples that are present in the input VCF.
* DumpSTR now reads call depth from the LC format field when the DP
    format field is not present. This was intended previously but was
    not happening.

### Documentation

* Clarified in PUBLISHING.rst how to handle dependencies and how to
    publish to bioconda.
* requirements.txt was unneeded, so delete it and remove the reference
    to it from PUBLISHING.rst

## 4.1.0

### Functionality Changes

* MergeSTR: Flanking basepairs are now removed from HipSTR records
    before merging. In particular, records with different flank lengths
    but the same repeat section will now merge instead of being marked
    as incompatible.
* CompareSTR: the tool now only compares records that start and end at
    the same position. If a partial overlap in records is detected, the
    program will output a warning to the user. This warning contains IDs
    of the records and their positions.

### Misc

* mergeutils: function GetMinHarmonizedRecords was transformed into
    GetIncrementAndComparability, which allows the caller to define
    custom predicate that decides whether records are comparable.

## 4.0.2

### Bug fixes

* <https://github.com/gymrek-lab/TRTools/issues/146> fixed record
    positions being compared twice
* CompareSTR: Decision on which records are comparable is now based on
    data from harmonized TRRecords, and not from the records directly
    from VCF readers. Thanks to this, HipSTR records which have
    different starting positions, but position of their repeat is at the
    same position are compared correctly (harmonization step removes
    this difference).
* MergeSTR failed on mixed ploidy samples (i.e. chrX). Fix one such
    bug. Note: none of the tools are fully tested for chrX even with
    this fix.

## 4.0.1

### Bug fixes

* <https://github.com/gymrek-lab/TRTools/issues/143> Fix
    HipstrMinSuppReads filter when there are called samples but none
    have ALLREADS

## 4.0.0

### Features

* Underlying libraries now use cyvcf2 instead of PyVCF for VCF
    parsing. This makes both the underlying VCF reading code and the
    TRTools code significantly faster and more memory efficient. For
    instance, the loading of VCFs into memory is now \> 15x faster for
    VCFs with many samples. Some tools will still need further updates
    to be usable for large datasets, but those updates should now be
    possible and much easier. (e.g. emitting progress reports to stdout
    as needed, flags to disable computations that cannot be done at such
    scale)
* DumpSTR has a new flag \--zip to produce a bgzipped and
    tabix-indexed output VCF
* StatSTR now can calculate the entropy of the allele distribution at
    each locus with the \--entropy flag
* The [TRTools documentation
    website](https://trtools.readthedocs.io/en/stable/) now displays the
    release notes.

### Command line interface changes

* StatSTR\'s \--region option now requires the input VCF to be
    bgzipped and tabix indexed.
* If DumpSTR is used on an input VCF with unexpectedly typed INFO
    fields \'AC\', \'REFAC\', \'HET\', \'HWEP\', \'HRUN\' or FORMAT
    field \'FILTER\', it now errors out and asks you to rename those
    fields before rerunning DumpSTR. (If they already exist but have the
    correct number and type DumpSTR will overwrite them and issue a
    warning in case that was not intended)
* CompareSTR\'s docs used to claim that when comparing alleles from
    different callers those callers must use the same allele notation
    (e.g. implying that ExpansionHunter\'s \'\<STR12\>\' and GangSTR\'s
    \'ACACACACACAC\' notation would always mismatch). That statement was
    never true for length based comparisons - CompareSTR has always been
    able to do length based comparisons regardless of notation. The
    incorrect claim has been removed from CompareSTR\'s docs.
* CompareSTR\'s docs now explicitly tell the user to order phased
    calls to prevent spurious mismatching. If phasing is not desired,
    use \--ignore-phasing
* CompareSTR will now error if at a single locus both files do not
    have either all unphased calls, or all phased calls. If phasing is
    not desired, use \--ignore-phasing

### Output changes

* DumpSTR call level filters now have the value of the filter and the
    value which triggered the filter appended to the filter name in the
    FILTER format field. (e.g. GangSTRCallMinDepth20_12 because the
    field had a depth of 12 and that\'s lower than the required min
    depth of 20)
* DumpSTR locus filter HRUN is now written as HRUN and not HRUN0 in
    the samplog output file
* When running DumpSTR, loci where all the calls were either already
    nocalls or were filtered by call-level filters before the
    locus-level filters were run are now marked as
    \'NO_CALLS_REMAINING\' instead of \'PASS\'.
* When DumpSTR filters a call and replaces each of its format fields
    with the no call \'.\', fields with more than one value are now
    represented correctly. For example, for 2 values \'.,.\' is used
    rather than just a single \'.\'
* MergeSTR header lines are now copied over from the input VCFs
    instead of only copying over a few recognized fields (e.g. ID and
    Length were the only contig fields that were previously retained,
    but URL wouldn\'t be)
* MergeSTR output alt alleles for eh and popstr are now ordered by
    length. MergeSTR output alt alleles for advntr, gangstr and hipstr,
    when there are multiple alt alleles of the same length, are now
    ordered alphabetically instead of arbitrarily.
* CompareSTR no longer outputs the file \<prefix\>-callcompare.tab -
    the existence of that file was never documented, and besides, all
    its information could be seen more easily simply by looking at the
    input VCFs
* In CompareSTR\'s overall.tab file, the ranges in the format columns
    are now written \[a,b) or \[a,b\] instead of a-b
* CompareSTR\'s locuscompare.tab file now outputs loci in the order
    they were encountered in the input VCfs as opposed to an arbitrary
    order
* The \'sample\' column in CompareSTR\'s locuscompare.tab file has
    been renamed to \'numcalls\' to match the other two tab files.

### Python interface changes

* The trtools.utils.tr_harmonizer module has been reworked to use
    cyvcf2, and in doing so a large portion of its interface has changed
    in small ways.
* The big conceptual change is that instead of repeatedly calling a
    method on a TRRecord object like GetStringGenotype for each sample
    in the VCF, instead you call the new corresponding method
    GetStringGenotypes once, and it returns a numpy array of values
    where the first axis of the array ranges over the samples.
* The way missing calls and samples with lower than maximal ploidy are
    handled is now tested and documented. These representations of these
    genotypes have been aligned with cyvcf2\'s standards. For more info,
    see the docs of the index, length and string genotype getter
    methods.

### Bug fixes

* The AC, REFAC fields that DumpSTR output used to be incorrect, are
    now correct
* If you specify \--drop-filtered DumpSTR will no longer set all
    values in the output .loclog.tab file to zero and instead set them
    to their proper values (which are the same as if you had not
    specified \--drop-filtered)
* DumpSTR now correctly adds ##FILTER=\<ID=PASS,Description=\"All
    filters passed\"\> to the header line
* DumpSTR now no longer says HipSTRCallFlankIndels is applied to
    nocalls
* MergeSTR now outputs the same phase as the input files instead of
    always outputting unphased data
* MergeSTR now correctly outputs Number=A, G or R (number of entries
    in this field equal to number of alternate alleles at this locus,
    the number of alleles including the ref, or the number of unique
    polyploid genotypes) correctly in INFO and FORMAT fields instead of
    outputing Number=-1, -2 or -3
* CompareSTR claimed it was outputting the square (Pearson)
    correlation coefficient but was actually outputting the raw
    (unsquared) correlation coefficient. It is now outputting the
    squared coefficient as documented.
* CompareSTR now correctly compares unphased calls without regard to
    order in the VCF (e.g. \'AAAA/AAA\' now matches against
    \'AAA/AAAA\')
* CompareSTR\'s docs claimed the bubble plots axes were measured in
    basepair difference from the reference, but they were actually
    measured in number of repeats different from the reference. The
    behavior has not been changed and the claim has been updated to
    match the behavior.
* When using binned format fields in CompareSTR where the range of
    values did not evenly divide into the requested binsize, the highest
    valued bin used to always be the same size as all the other bins and
    include values over the limit specified by the user. Now it caps at
    that maximum. E.g. binsizes 0:210:50 used to create the bins
    \[0,50), \[50,100), \[100,150), \[150, 200), \[200, 250) and now
    create the bins \[0,50), \[50,100), \[100,150), \[150, 200), \[200,
    210\]
* When using binned format fields in CompareSTR where the range of
    values evenly divided into the requested binsize, loci which
    obtained the requested maximum would be excluded. They are now
    included. E.g. binsizes 0:200:50 used to create the bins \[0,50),
    \[50,100), \[100,150), \[150, 200) and samples with value 200 would
    not fall into any bin. This now creates the bins \[0,50), \[50,100),
    \[100,150), \[150, 200\] and samples with value 200 fall into the
    last bin

### Quality of life improvements

* StatSTR, when printing output to a file, now prints timing
    diagnostics to stdout.
* DumpSTR will fail faster if output directory does not exist
* When encountering issues with identifying the caller type for each
    input VCF, MergeSTR now prints an error and gracefully returns
    instead of dying to an uncaught exception
* MergeSTR incompatible INFO field warnings now specify which locus
    has an incompatible field

### Regressions

* The \--gangstr-require-support filter has been disabled.

### Outstanding bugs

* The dumpSTR ExpansionHunter ADFL ADIR ADSP filters have never worked
* DumpSTR remains untested on ExpansionHunter filters and files
* DumpSTR remains untested on loci with variable ploidy and/or
    partially genotyped samples (e.g. .\|2)
* When running CompareSTR with the \--stratify options where
    \--stratify-file is either not specified or is explicitly set to
    zero, for each format field all calls where the value of that field
    in vcf1 does not fall into the same bin as the value of that field
    in vcf2 are silently not compared for that format field. The correct
    behavior here is probably to create paired bins based on a range of
    values from vcf1 and a range from vcf2. Regardless, the behavior
    here should be documented.

## 3.0.3

### Bug fixes

* Fixed a spot where qcSTR would crash because we passed Pandas a set
    instead of a list
* MergeSTR now writes out the header for the GT FORMAT field
