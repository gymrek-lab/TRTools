.. overview_directive
.. |associaTR overview| replace:: associaTR performs association testing of the lengths of TRs against phenotypes.
.. overview_directive_done

AssociaTR
==========

|associaTR overview|

This tool is designed for running a GWAS of TRs against a single outcome. While slower, it can also be invoked once for each phenotype
to perform a PHWAS of a group of TRs against many phenotypes. Many tools perform association testing and are more fully featured than
associaTR - the main value of associaTR is its ability to conveniently perform length based tests for multiallelic TR loci instead of
per-allele based tests.

associaTR is in **beta**. In particular, while its main outputs (p-values, coefficients and standard errors) are all fully tested,
it has not been fully code-reviewed, and other outputs still need testing. It also does not currently support binary phenotypes.

Importantly, associaTR works with Beagle dosages, as imputation dosages are more reliable than imputed best-guess genotypes for association
testing.

Usage
-----

Required positional parameters: 

* :code:`outfile` - the location to write the association results tsv to.
* :code:`tr_vcf` - the location of vcf containing the TR genotypes to test
* :code:`phenotype_name` - the name of the phenotype being tested against, used to write appropriately named column headers
* :code:`traits` - At least one `.npy <https://numpy.org/doc/stable/reference/generated/numpy.save.html>`_
  2d float array file containing trait values for samples. The first trait is the phenotype to be regressed
  against, all other traits from that file are used as covariates. If more than one file is provided, then
  all traits in the additional files are added as additional covariates. 

  If :code:`--same-samples` is not specified, the first column of each file must be the numeric sample ID,
  designating which row corresponds to which sample. (Thus the phenotype will correspond to the second column
  from the first file.) If multiple files are specified, they will be joined on sample ID. 

  If :code:`--same-samples` is specified, there must be the same number of rows in each array as the number 
  of samples in the vcf, those will be the traits for those samples.
  (Thus the phenotype will correspond to the first column of the first array). If 
  multiple files are specified, then they will be concatenated horizontally. Since IDs do not need to be stored 
  in the npy arrays, :code:`--same-samples` allows for non-numeric sample IDs. 
 
Optional input parameters:   

* :code:`--vcftype` One of :code:`eh, hipstr, gangstr, popstr, advntr`
  Specify which caller produced the TR VCF, useful when the VCF is ambiguous 
  and the caller cannot be automatically inferred.
* :code:`--same-samples` - see the description of traits above
* :code:`--beagle-dosages` - regress against Beagle dosages from the AP{1,2} fields instead of from the GT field. 
  Note: The GP field that Beagle outputs is not supported

Optional filtering parameters:

* :code:`--sample-list` File containing list of samples to use, one sample ID per line. 
  If not specified, all samples are used.
* :code:`--region` Restrict to chr:start-end
* :code:`--non-major-cutoff` - If not :code:`--beagle-dosages`, then this is the non-major-allele-count cutoff. 
  I.e. filter all loci with :code:`non-major-allele-count < cutoff.`
  If working with dosages, this cutoff is applied to the dosage sums. 
  As with the regression itself, for this purpose alleles are coallesced by length. 
  Note that this is a raw value cutoff, not a percentage cutoff.
  This defaults to 20 per `plink's best practices <https://www.cog-genomics.org/plink/2.0/assoc#glm>`_
  Set to 0 to disable this filter. 

Regression details
------------------

AssociaTR performs ordinary least squares regression, i.e. it fits the model

  y = g*beta + Xb + error

where `y` is the outcome, `g` are the TR length genotypes, `beta` is the coefficient associating 
`y` and `g`, `X` are the other covariates, `b` are the coefficients of those covariates, and `error`
is a normally distributed error term. associaTR reports the p-value corresponding to the hypothesis
:code:`beta != 0`, `beta` and the standard error of `beta`. It also reports the `R^2` value of this linear
model.

If not using dosages, then `g` is simply :code:`sum(len(allele_1) + len(allele_2))` for each sample.
If using dosages, then `g` is :code:`Σ_{lengths L} L*[Pr(len(allele_1) == L) + Pr(len(allele_2) == L)]`  for each sample.

An intercept term is always automatically included in `X` - do not include one yourself.
`y`, `g` and `X` are always standardized (transformed to mean zero and variance 1). So it is not necessary for you to
pre-standardize `y` and `X`. Coefficients and standard errors are reported in the original units of the provided traits,
not on the standardized units. However, if you pre-standardize your traits, then they will by necessity be reported
in those standardized units. The coefficient and standard error are reported in change per number of repeat copies.
To get them measured in change per basepair, divide by the length of the repeat unit.

Outputs
-------

* :code:`chrom` - chromosome of the TR
* :code:`pos` - start position in bp of the TR
* :code:`alleles` - the lengths of the TR in the input VCF
* :code:`n_samples_tested` - number of samples tested. This is the number of samples in the VCF, restricted by:
  
  * If :code:`--same-samples` is not specified, then there can be fewer samples in the trait arrays. If so,
    samples not in those arrays will be omitted.
  * If :code:`--sample-list` is specified, samples not in that list will be omitted.
  * Samples with missing genotypes will be omitted on a per-variant basis.

* :code:`locus_filtered` - False if the locus was not filtered, or a reason for filtering otherwise. Possible reasons:

  * No called samples
  * Only one called length (if dosages, this implies all other lengths have exactly zero dosage)
  * This locus did not pass the optional cutoff specified by :code:`--non-major-cutoff`
  * The number of samples remaining is less than or equal to the number of covariates in the model (including
    the genotypes and the intercept term).

* :code:`p_{phenotype_name}` - the p-value of the association or nan if the locus was filtered
* :code:`coeff_{phenotype_name}` - the coefficient of the association, measured in phenotype units per repeat copy
  or nan if the locus was filtered.
* :code:`se_{phenotype_name}` - the standard error of the coefficient or nan if the locus was filtered
* :code:`regression_R^2` - the R^2 of the fitted model or nan if the locus was filtered
* :code:`motif` - the TR's motif
* :code:`period` - the TR's period
* :code:`ref_len` - the length of the reference allele for the TR
* :code:`allele_frequency` - the frequency of each length allele in the tested samples
* :code:`dosage_esimtated_r2_per_length_allele` - Only applicable to :code:`--beagle-dosages`. This is a dictionary of r^2 values, one for
  each length, correlating each haplotype's imputed probability of obtaining that length vs whether or not the imputed best-guess
  genotype is that length
* :code:`r2_length_dosages_vs_best_guess_lengths` - Only applicable to :code:`--beagle-dosages`. This is the r^2 value for the correlation
  between each haplotype's average imputed length (:code:`Σ_{lengths L} L*[Pr(len(allele) == L)`) with the best-guess
  length :code:`L` of that haplotype.

Example Commands
----------------

Below is an :code:`associaTR` example. For this example no TRs causally impact the simulated phenotype.
Data files for this example can be found at https://github.com/gymrek-lab/TRTools/tree/master/example-files::

  associaTR \
    association_results.tsv \
    ceu_ex.vcf.gz \
    simulated_phenotype \
    simulated_traits_0.npy \
    simulated_traits_1.npy \
    --same-samples

