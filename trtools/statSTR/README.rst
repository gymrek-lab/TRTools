.. overview_directive
.. |statSTR overview| replace:: StatSTR takes in a TR genotyping VCF file and outputs per-locus statistics.
.. overview_directive_done


StatSTR
=======

|statSTR overview|

Usage
-----
To run statSTR use the following command::

   statSTR
      --vcf <vcf file> \
      --out <prefix> \
      --thresh \
      [additional options]

Required parameters:

* :code:`--vcf <string>`: input the STR VCF file
* :code:`--vcftype <string>`: Type of VCF file being processed. Default = :code:`auto` Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.
* :code:`--out <string>`: prefix to name output files. Set to stdout to write to standard output.

Optional filtering group parameters:

* :code:`--samples <string>`: A file containing a list of samples to include in computing statistics. If not given, all samples are used. To compute statistics for multiple groups of samples, you can give a comma-separated list of samples files.
* :code:`--sample-prefixes <string>`: Prefixes to name output for each samples group. By default uses 1,2,3 etc. Must be sample length as :code:`--samples`.
* :code:`--region <string>`: restrict to specific regions (chrom:start-end).

For specific statistics available, see below.

Output file
-----------

StatSTR outputs a tab-delimited file called :code:`<prefix>.tab` with per locus statistics.
It has columns: chrom, start, end, and an additional column for each statistic specified.

Statistics options
------------------

The following options can be specified to compute various per-locus statistics. The column
produced in the output file has the same name as the specified option:

* :code:`--thresh`: Output the maximum observed allele length.
* :code:`--afreq`: Output allele frequences. Comma-separated list of allele:freq
* :code:`--acount`: Output allele counts. Comma-separated list of allele:count
* :code:`--hwep`: Output Hardy Weinberg p-values per locus. 
* :code:`--het`: Output heterozygosity of each locus.
* :code:`--entropy`: Output the entropy of the distribution of genotypes at each locus in bits.
  (Entropy is a measurement of how many genotypes there are at a locus and how 
  evenly distributed they are. Higher entropy means more genotypes more evenly distributed.
  It is a more accurate measurement of the complexity of a locus than heterozygosity.
  See `wikipedia <https://en.wikipedia.org/wiki/Information_content>`_ for more info.)
* :code:`--mean`: Output mean allele length. 
* :code:`--mode`: Output mode allele length.
* :code:`--var`: Output variance of allele length.
* :code:`--numcalled`: Output number of called samples.

For genotypers which output allele sequences, :code:`--use-length` will collapse alleles by length.
(This is implicit for genotypers which only output allele lengths.)

Example command
---------------

If you have installed the TRTools package, you can run an example statSTR command using files in this repository::

  statSTR \
    --vcf $REPODIR/test/common/sample_vcfs/test_hipstr.vcf \
    --out teststats \
    --afreq

where :code:`$REPODIR` points to the root path of this repository.

If you are testing from source, you can run::

  python statSTR.py \
    --vcf $REPODIR/test/common/sample_vcfs/test_hipstr.vcf \
    --out teststats \
    --afreq

This will output :code:`teststats.tab` with columns chrom, start, end, afreq.

