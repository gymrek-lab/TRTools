.. overview_directive
.. |mergeSTR overview| replace:: MergeSTR merges multiple VCF files produced by the same TR genotyper into a single VCF file.
.. overview_directive_done

MergeSTR
========

|mergeSTR overview|

If TR genotyping was performed separately on different samples or batches of samples, mergeSTR can be used to combine the resulting VCFs into one file. This is often necessary for downstream steps such as: computing per-locus statistics, performing per-locus filtering, and association testing.

While other VCF libraries have capabilities to merge VCF files, they do not always handle multi-allelic TRs properly, especially if the allele definitions are different across files. MergeSTR is TR-aware and currently handles VCF files obtained by: GangSTR, HipSTR, ExpansionHunter, popSTR, or adVNTR. See below for specific VCF fields supported for each genotyper.

Note: mergeSTR does not support merging VCFs produced by different TR genotypers as the desired outcome of such an operation is highly dependant on the use-case at hand.
If this is your use case and you with it to be supported, consider `filing an issue <https://github.com/gymreklab/TRTools/issues>`_ with us.

Usage
-----
MergeSTR takes as input two or more VCF files with TR genotypes and outputs a combined VCF file. Note, input VCF files must be bgzipped, sorted, indexed, and have the appropriate :code:`##contig` header lines. See `Instructions on Compressing and Indexing VCF files`_ below for commands for preparing tabix-indexed VCF files.

To run mergeSTR use the following command::

	mergeSTR \
  	  --vcfs <file1.vcf,file2.vcf,...> \
  	  --out test \
  	  [additional options]

Required Parameters:

* :code:`--vcf <VCF>`: Comma-separated list of VCF files to merge. All must have been created by the same TR genotyper. Must be bgzipped, sorted, and indexed.
* :code:`--vcftype <string>`: Type of VCF files being merged. Default = :code:`auto`. Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.
* :code:`--out <string>`: prefix to name output files

Special Merge Options:

* :code:`--update-sample-from-file`: Append file names to sample names. Useful if sample names are repeated across VCF files.

Optional Additional Parameters:

* :code:`--verbose`: Prints out extra information
* :code:`--quiet`: Doesn't print out anything

Example MergeSTR command
------------------------

If you have installed the TRTools package, you can run an example mergeSTR command using files in this repository::

	FILE1=${REPODIR}/test/common/sample_vcfs/mergeSTR_vcfs/test_file_gangstr1.vcf.gz
	FILE2=${REPODIR}/test/common/sample_vcfs/mergeSTR_vcfs/test_file_gangstr2.vcf.gz
	mergeSTR \
   	--vcfs ${FILE1},${FILE2} \
   	--out test_run

where :code:`$REPODIR` points to the root path of this repository.

If you are testing from source, you can run::

     python mergeSTR.py \
   	--vcfs ${FILE1},${FILE2} \
   	--out test_run

This command should create a file :code:`test_run.vcf` with the merged genotypes.

Supported VCF fields
--------------------

In addition to proper merging of alleles at multi-allelic sites, MergeSTR supports the following VCF fields for each tool. Fields not listed are currently ignored when merging. INFO fields below are expected to be constant across loci being merged.

**AdVNTR**

* Supported INFO fields: END, RU, RC
* Supported FORMAT fields: DP, SR, FL, ML

**ExpansionHunter**

* Supported INFO fields: END, REF, REPID, RL, RU, SVTYPE
* Supported FORMAT fields: ADFL,ADIR,ADSP,LC,REPCI,REPCN,SO

**GangSTR**

* Supported INFO fields: END, RU, PERIOD, REF, EXPTHRESH
* Supported FORMAT fields: DP,Q,REPCN,REPCI,RC,ENCLREADS,FLNKREADS,ML,INS,STDERR,QEXP

**HipSTR**

* Supported INFO fields: START, END, PERIOD
* Supported FORMAT fields: GB,Q,PQ,DP,DSNP,PSNP,PDP,GLDIFF,DSTUTTER,DFLANKINDEL,AB,FS,DAB,ALLREADS,MALLREADS

**PopSTR**

* Supported INFO fields: Motif
* Supported FORMAT fields: AD, DP, PL

Instructions on Compressing and Indexing VCF files
--------------------------------------------------
MergeSTR requires the input file to be compressed and indexed. Use the following commands to create compressed and indexed vcf file::

  bgzip file.vcf
  tabix -p vcf file.vcf.gz


