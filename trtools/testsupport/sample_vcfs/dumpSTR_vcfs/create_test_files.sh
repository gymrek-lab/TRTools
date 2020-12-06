# Note: this command should be run from the root of the repo

# Note: the thresholds in this file are arbitrary and just used for testing
# do not assume they are reasonable for data analysis

out=trtools/testsupport/sample_vcfs/dumpSTR_vcfs

# locus_filters
echo "locus filters"
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/trio_chr21_hipstr.sorted.vcf.gz --out $out/locus_filters --min-locus-callrate 0.5 --min-locus-hwep 0.5 --min-locus-het 0.05 --max-locus-het 0.45 --filter-regions-names foo_region --filter-regions ../repo/trtools/testsupport/sample_vcfs/dumpSTR_vcfs/sample_region.bed.gz --vcftype hipstr

echo "drop filtered"
# same as above test, only difference should be the vcf file, so delete the other two
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/trio_chr21_hipstr.sorted.vcf.gz --out $out/drop_filtered --min-locus-callrate 0.5 --min-locus-hwep 0.5 --min-locus-het 0.05 --max-locus-het 0.45 --filter-regions-names foo_region --filter-regions ../repo/trtools/testsupport/sample_vcfs/dumpSTR_vcfs/sample_region.bed.gz --vcftype hipstr --drop-filtered
rm $out/drop_filtered.samplog.tab
rm $out/drop_filtered.loclog.tab

# advntr_filters
echo "advntr"
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/NA12878_chr21_advntr.sorted.vcf.gz --out $out/advntr_filters --advntr-min-call-DP 50 --advntr-max-call-DP 2000  --advntr-min-spanning 1 --advntr-min-flanking 20 --advntr-min-ML 0.95

# eh_filters
# TODO some of the EH filters never worked in the first place
# time python -m trtools.dumpSTR.dumpSTR --vcf example-files/NA12878_chr21_eh.sorted.vcf.gz --out $out/eh_filters --eh-min-ADFL 3 --eh-min-ADIR 3 --eh-min-ADSP 1 --eh-min-call-LC 50 --eh-max-call-LC 1000

# gangstr_filters
echo "gangSTR"
# second one has require support removed, reenable after it is recreated
time python -m trtools.dumpSTR.dumpSTR --vcf trtools/testsupport/sample_vcfs/test_gangstr.vcf --out $out/gangstr_filters_expansion --gangstr-expansion-prob-het 0.001 --gangstr-expansion-prob-hom 0.0005 --gangstr-expansion-prob-total 0.001
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/trio_chr21_gangstr.sorted.vcf.gz --out $out/gangstr_filters_most --gangstr-min-call-DP 10 --gangstr-max-call-DP 100 --gangstr-min-call-Q 0.9  --gangstr-filter-span-only --gangstr-filter-spanbound-only --gangstr-filter-badCI # --gangstr-require-support 10 --gangstr-readlen 150

# hipstr_filters
echo "hipSTR"
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/trio_chr21_hipstr.sorted.vcf.gz --out $out/hipstr_filters --filter-hrun --use-length --max-locus-het 0.45 --min-locus-het 0.05 --min-locus-hwep 0.5 --hipstr-max-call-flank-indel 0.05 --hipstr-max-call-stutter 0.3 --hipstr-min-supp-reads 10 --hipstr-min-call-DP 30 --hipstr-max-call-DP 200 --hipstr-min-call-Q 0.9 --vcftype hipstr

# popstr_filters
echo "popSTR"
time python -m trtools.dumpSTR.dumpSTR --vcf example-files/NA12878_chr21_popstr.sorted.vcf.gz --out $out/popstr_filters --popstr-min-call-DP 30 --popstr-max-call-DP 200 --popstr-require-support 15
