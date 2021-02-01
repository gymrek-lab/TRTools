There is a copy of each of these images in doc/source/images. If changed here
please update them there (as well as the READMEs in both folders)

diffref-bias, diffref-histogram and sample-callnum were produced by running qcSTR on
example-files/trio_chr21_popstr.sorted.vcf.gz

chrom-callnum was produced by running qcSTR on
trtools/testsupport/sample_vcfs/many_samples_multiple_chroms.vcf.gz

quality-per-locus and quality-per-sample were produced by running qcSTR with 
the corresponding options on trtools/testsupport/sample_vcfs/many_samples.vcf.gz

quality-locus-stratified is very similar to running qcSTR with 
the corresponding options on trtools/testsupport/sample_vcfs/qc_vcfs/5by50.vcf

quality-per-call and quality-sample-stratified were generated on a popstr vcf
which we no longer interpret as having quality scores, so they are not currently
reproducible and should be replaced eventually
