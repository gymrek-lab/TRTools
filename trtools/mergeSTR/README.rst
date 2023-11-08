.. overview_directive
.. |mergeSTR overview| replace:: MergeSTR merges multiple VCF files produced by the same TR genotyper into a single VCF file.
.. overview_directive_done

MergeSTR
========

|mergeSTR overview|

If TR genotyping was performed separately on different samples or batches of samples, mergeSTR can be used to combine the resulting VCFs into one file. This is often necessary for downstream steps such as: computing per-locus statistics, performing per-locus filtering, and association testing.

While other VCF libraries have capabilities to merge VCF files, they do not always handle multi-allelic TRs properly, especially if the allele definitions are different across files. MergeSTR is TR-aware. See below for specific VCF fields supported for each genotyper.

mergeSTR does not support merging VCFs produced by different TR genotypers - that is a more complex usecase, and we are designing a separate tool to do that.

Usage
-----
MergeSTR takes as input two or more VCF files with TR genotypes and outputs a combined VCF file. Note, input VCF files must be bgzipped, sorted, indexed, and have the appropriate :code:`##contig` header lines. See `Instructions on Compressing and Indexing VCF files`_ below for commands for preparing tabix-indexed VCF files.

To run mergeSTR use the following command::

	mergeSTR \
  	  --vcfs <file1.vcf.gz,file2.vcf.gz,...> \
  	  --out <string> \
  	  [additional options]

Required Parameters:

* :code:`--vcf <VCF>`: Comma-separated list of VCF files to merge. All must have been created by the same TR genotyper. Must be bgzipped, sorted, and indexed. (See `Instructions on Compressing and Indexing VCF files`_ below)
* :code:`--vcftype <string>`: Type of VCF files being merged. Default = :code:`auto`. Must be one of: :code:`gangstr`, :code:`advntr`, :code:`hipstr`, :code:`eh`, :code:`popstr`.
* :code:`--out <string>`: prefix to name output files

Special Merge Options:

* :code:`--update-sample-from-file`: Append file names to sample names. Useful if sample names are repeated across VCF files.

Optional Additional Parameters:

* :code:`--verbose`: Prints out extra information
* :code:`--quiet`: Doesn't print out anything

MergeSTR outputs a merged VCF file :code:`$out.vcf` with the merged genotypes. See `Example Commands`_ below for example mergeSTR commands for different supported TR genotypers based on example data files in this repository.

Note: when running MergeSTR on HipSTR data, MergeSTR first strips the flanking basepairs from all alleles before merging, allowing records with different flanking basepairs but the same repeat to be comparable.
Be aware that in the infrequent occurences where there are indels in the flanking regions called by HipSTR this will cause incorrect alleles to be outputted by the merge.

Supported VCF fields
--------------------

In addition to proper merging of alleles at multi-allelic sites, MergeSTR supports the following VCF fields for each tool. Fields not listed are currently ignored when merging. INFO fields below are expected to be constant across loci being merged.

**AdVNTR**

* Supported INFO fields: END, RU, RC, VID
* Supported FORMAT fields: DP, SR, FL, ML

**ExpansionHunter**

* Supported INFO fields: END, REF, REPID, RL, RU, SVTYPE, VARID
* Supported FORMAT fields: ADFL, ADIR, ADSP, LC, REPCI, REPCN, SO

**GangSTR**

* Supported INFO fields: END, RU, PERIOD, REF, EXPTHRESH
* Supported FORMAT fields: DP, Q, REPCN, REPCI, RC, ENCLREADS, FLNKREADS, ML, INS, STDERR, QEXP

**HipSTR**

* Supported INFO fields: START, END, PERIOD
* Supported FORMAT fields: GB, Q, PQ, DP, DSNP, PSNP, PDP, GLDIFF, DSTUTTER, DFLANKINDEL, AB, FS, DAB, ALLREADS, MALLREADS

**PopSTR**

* Supported INFO fields: Motif
* Supported FORMAT fields: AD, DP, PL

Instructions on Compressing and Indexing VCF files
--------------------------------------------------
MergeSTR requires the input file to be compressed and indexed. Use the following commands to create compressed and indexed vcf file::

  bgzip file.vcf
  tabix -p vcf file.vcf.gz

Example Commands
----------------

Below are :code:`mergeSTR` examples using VCFs from supported TR genotypers. Data files can be found at https://github.com/gymrek-lab/TRTools/tree/master/example-files::

  # AdVNTR
  # Note, you first need to reheader files to add required contig lines to VCF headers
  for sample in NA12878 NA12891 NA12892; do
      bcftools reheader -f hg19.fa.fai -o ${sample}_advntr_reheader.vcf.gz ${sample}_chr21_advntr.sorted.vcf.gz
      tabix -p vcf ${sample}_advntr_reheader.vcf.gz
  done
  FILE1=NA12878_advntr_reheader.vcf.gz
  FILE2=NA12891_advntr_reheader.vcf.gz
  FILE3=NA12892_advntr_reheader.vcf.gz
  mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out test_merge_advntr --vcftype advntr --update-sample-from-file # outputs test_merge_advntr.vcf

  # ExpansionHunter
  # Note, you first need to reheader files to add required contig lines to VCF headers
  for sample in NA12878 NA12891 NA12892; do 
      bcftools reheader -f hg19.fa.fai -o ${sample}_eh_reheader.vcf.gz ${sample}_chr21_eh.sorted.vcf.gz
      tabix -p vcf ${sample}_eh_reheader.vcf.gz
  done
  FILE1=NA12878_eh_reheader.vcf.gz
  FILE2=NA12891_eh_reheader.vcf.gz
  FILE3=NA12892_eh_reheader.vcf.gz
  mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out test_merge_eh --vcftype eh # outputs test_merge_eh.vcf

  # GangSTR
  FILE1=NA12878_chr21_gangstr.sorted.vcf.gz
  FILE2=NA12891_chr21_gangstr.sorted.vcf.gz
  FILE3=NA12892_chr21_gangstr.sorted.vcf.gz
  mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out test_merge_gangstr --vcftype gangstr # outputs test_merge_gangstr.vcf

  # HipSTR
  FILE1=NA12878_chr21_hipstr.sorted.vcf.gz
  FILE2=NA12891_chr21_hipstr.sorted.vcf.gz
  FILE3=NA12892_chr21_hipstr.sorted.vcf.gz
  mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out test_merge_hipstr --vcftype hipstr # outputs test_merge_hipstr.vcf

  # PopSTR
  FILE1=NA12878_chr21_popstr.sorted.vcf.gz
  FILE2=NA12891_chr21_popstr.sorted.vcf.gz
  FILE3=NA12892_chr21_popstr.sorted.vcf.gz
  mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out test_merge_popstr --vcftype popstr # outputs test_merge_popstr.vcf
