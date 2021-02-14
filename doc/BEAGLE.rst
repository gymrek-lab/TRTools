.. _working_with_beagle:

Working with Calls Imputed by Beagle
====================================

TRTools has preliminary support for tandem repeats imputed by Beagle_.
Using a reference panel that contains both SNP and TR genotypes, Beagle can impute
TR genotypes into samples on which only SNP genotypes were directly called.
Note that TRTools
only support Beagle imputed calls which were originally called by one of the
:ref:`supported TR callers <callers>`.

To run TRTools with a properly formatted Beagle-generated VCF, simply
pass it to the utility and add the flag ``--vcftype beagle-CALLER``
where ``CALLER`` is the name of the caller that originally called the
TRs. For example, ``--vcftype beagle-hipstr``.
To properly format the Beagle generate VCF, keep reading.

..
    For example, the SNPSTR reference panel located
    `here <https://gymreklab.com/2018/03/05/snpstr_imputation.html>`_
    and published in
    `this paper <https://www.nature.com/articles/s41467-018-06694-0>`_ by Saini et al.
    contains HipSTR-called TR genotypes genome-wide. Beagle can impute those HipSTR
    calls into other samples with only SNP genotypes.

Invoking Beagle
---------------
By default, Beagle only outputs hardcalls and allele dosages for each imputed call
which are insufficient for some types of call level filtering and quality control.
Add either of the flags ``ap=true`` or ``gp=true`` to your call to Beagle to
have it emit allele-specific probability scores that can be filtered on. (
Using ``ap`` will cause the resulting (uncompressed) VCF to increase in size
by about 3 times. Using ``gp`` can cause it to increase in size quadratically in
the number of alleles at each locus. As such, we recommend ``ap``.)

The VCFs Beagle outputs must be preprocessed before being used by TRTools.
These VCFs will have all the loci that were in the reference panel -
both SNP and TR loci. Furthermore, the TR loci will not have any of the ``INFO`` fields
necessary to describe the locus. To fix these issues:

* Use ``bcftools view`` with the flags ``-i ID=@$ID_FILENAME`` or ``-r $LOCUS_FILENAME``
  to output a VCF with only the TR loci (subsetting by either ID or locus,
  respectively).
* Use ``bcftools annotate`` with the flags ``--collapse none``, ``--single-overlaps``
  and ``-c INFO/FIELD1,INFO/FIELD2,...`` to copy the info fields over from the reference
  panel to the VCF file with TR loci.

The ``INFO`` fields that need to be copied over for each TR genotyper are:

* HipSTR: ``START, END, PERIOD``

Citation
--------
If you use Beagle in an analysis you publish, please
`cite it <https://faculty.washington.edu/browning/beagle/beagle.html#citation>`_
per Beagle's authors' request.

Additional work needed
----------------------
* Support Beagle imputed TRs from callers other than HipSTR
* Support Beagle's ``gp`` field
* Add support for QC based on dosage in addition to call probability?
* Publish the INFO fields associated with the SNPSTR reference panel
  so we can refer people to that. Until then, calls imputed by it cannot
  be used in TRTools.

.. _Beagle: https://faculty.washington.edu/browning/beagle/beagle.html
