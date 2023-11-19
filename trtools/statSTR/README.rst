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

* :code:`--vcf <string>`: The input TR VCF file
* :code:`--out <string>`: The prefix to name output files. Set to stdout to write to standard output.

Optional general parameters:

* :code:`--vcftype <string>`: The type of VCF file being processed. Default = :code:`auto` Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.
* :code:`--samples <string>`: A file containing a list of samples to include in computing statistics. If not given, all samples are used. To compute statistics for multiple groups of samples, you can give a comma-separated list of samples files. Sample files should list one sample per line, no header line. Samples not found in the VCF are silently ignored.
* :code:`--sample-prefixes <string>`: The prefixes to name output for each samples group. By default uses 1, 2, 3 etc. Must be sample length as :code:`--samples`.
* :code:`--region <string>`: Restrict to specific regions (chrom:start-end). Requires the input VCF to be bgzipped and tabix indexed.
* :code:`--precision <int>`: How much precision to use when writing stats (default = 3)

For specific statistics available, see below.

StatSTR outputs a tab-delimited file :code:`<outprefix>.tab` with per locus statistics. See a description of the output file below.

See `Example Commands`_ below for example statSTR commands for different supported TR genotypers based on example data files in this repository.

Statistics options
------------------

The following options can be specified to compute various per-locus statistics. The column
produced in the output file has the same name as the specified option:

* :code:`--thresh`: Output the maximum observed allele length.
* :code:`--afreq`: Output allele frequences. Comma-separated list of allele:freq.
  Only called alleles are included in the list. If there are no called alleles, '.' is emitted.
* :code:`--acount`: Output allele counts. Comma-separated list of allele:count
  Only called alleles are included in the list. If there are no called alleles, '.' is emitted.
* :code:`--nalleles`: Output the number of alleles per locus with at least a given frequency.
* :code:`--nalleles-thresh`: Threshold for the frequency cutoff for :code:`--nalleles`. Defaults to :code:`1%` if not specified.
* :code:`--hwep`: Output Hardy Weinberg p-values per locus.
* :code:`--het`: Output heterozygosity of each locus.
* :code:`--entropy`: Output the bit-entropy of the distribution of alleles at each locus.
  Entropy is a measurement of how hard it is to predict genotypes at a locus, where higher
  entropy values indicate more complex loci. See
  `wikipedia <https://en.wikipedia.org/wiki/Information_content>`_ for the mathematical definition
  of entropy.
* :code:`--mean`: Output mean allele length.
* :code:`--mode`: Output mode allele length.
* :code:`--var`: Output variance of allele length.
* :code:`--numcalled`: Output number of called samples.

The statistics :code:`thresh, hewp, het, entropy, mean, mode, var` will output :code:`nan` if there are no called alleles.
The statistics :code:`afreq, acount` will output :code:`.` if there are no called alleles.
The statistic :code:`nalleles` will output :code:`0` if there are no called alleles.

For genotypers which output allele sequences, :code:`--use-length` will collapse alleles by length.
(This is implicit for genotypers which only output allele lengths.)

Output file
-----------

StatSTR outputs a tab-delimited file with columns ``chrom``, ``start`` and ``end`` plus an additional column for each statistic specified.
If multiple sample groups are specified, instead there is one additional column for each group-by-statistic pair

Example Commands
----------------

Below are :code:`statSTR` examples using VCFs from supported TR genotypers. Data files can be found at https://github.com/gymrek-lab/TRTools/tree/master/example-files::

  # AdVNTR
  statSTR --vcf NA12878_chr21_advntr.sorted.vcf.gz \
        --out stdout \
        --afreq

  # ExpansionHunter
  statSTR --vcf NA12891_chr21_eh.sorted.vcf.gz \
        --out stats_eh \
        --numcalled

  # GangSTR
  statSTR --vcf trio_chr21_gangstr.sorted.vcf.gz \
         --out stats_gangstr \
         --numcalled \
         --mean

  # HipSTR
  statSTR --vcf trio_chr21_hipstr.sorted.vcf.gz \
        --vcftype hipstr \
        --out stats_hipstr \
        --acount \
        --afreq \
        --mean

  # PopSTR
  statSTR --vcf trio_chr21_popstr.sorted.vcf.gz \
        --out stats_popstr \
        --mean \
        --samples ex-samples.txt
