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
      --out <string> \
      --thresh \
      [additional options]

Required Parameters: 

* :code:`--vcf <string>`: the input TR VCF file
* :code:`--vcftype <string>`: Type of VCF file being processed. Default = :code:`auto` Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.
* :code:`--out <string>`: prefix to name output files. Set to stdout to write to standard output.

Filtering Group parameters: 

* :code:`--samples <string>`: A file containing a list of samples to include in computing statistics. If not given, all samples are used. To compute statistics for multiple groups of samples, you can give a comma-separated list of samples files.
* :code:`--sample-prefixes <string>`: Prefixes to name output for each samples group. By default uses 1,2,3 etc. Must be sample length as :code:`--samples`.
* :code:`--region <string>`: restrict to specific regions (chrom:start-end). 

For specific statistics available, see below.

StatSTR will output a tab delimited file :code:`$out.tab` with per-locus statistics. See a description of the output file below.

See `Example Commands`_ below for example statSTR commands for different supported TR genotypers based on example data files in this repository.

Statistics options
------------------

The following options can be specified to compute various per-locus statistics:

* :code:`--thres`: Output the maximum observed allele length. (Output column="thresh") 
* :code:`--afreq`: Output allele frequences. (Output column="afreq"). Comma-separated list of allele:freq  
* :code:`--acount`: Output allele counts. (Output column="acount"). Comma-separated list of allele:count  
* :code:`--hwep`: Output Hardy Weinberg p-values oer locus. (Output column="hwep") 
* :code:`--het`: Output heterozygosity of each locus. (Output column="het") 
* :code:`--mean`: Output mean allele length. (Output column="mean") 
* :code:`--mode`: Output mode allele length. (Output column="mode") 
* :code:`--var`: Output variance of allele length. (Output column="var") 
* :code:`--numcalled`: Output number of called samples. (Output column="numcalled") 

Specifying :code:`--use-length` will collapse alleles by length (only relevant for HipSTR, which outputs sequence differences in TR alleles).

Output file
-----------

StatSTR outputs a tab-delimited file with columns: chrom, start, end, plus an additional column for each statistic specified using the options described above.

Example Commands
----------------

Below are :code:`statSTR` examples using VCFs from supported TR genotypers. Data files can be found at https://github.com/gymreklab/TRTools/tree/master/example-files::

  # AdVNTR
  statSTR --vcf NA12878_chr21_advntr.sorted.vcf.gz --out stdout --afreq

  # ExpansionHunter
  statSTR --vcf NA12891_chr21_eh.sorted.vcf.gz --out stats_eh --numcalled

  # GangSTR
  statSTR --vcf trio_chr21_gangstr.sorted.vcf.gz --out stats_gangstr --numcalled --mean

  # HipSTR
  statSTR --vcf trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr --out stats_gangstr --acount --afreq --mean

  # PopSTR
  statSTR --vcf trio_chr21_popstr.sorted.vcf.gz --out stats_popstr --mean --samples ex-samples.txt
