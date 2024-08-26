#!/usr/bin/env bash

# This command should be run from the root of the repo

python -m trtools.annotaTR.annotaTR \
	--vcf trtools/testsupport/sample_vcfs/dumpSTR_vcfs/trio_chr21_gangstr.sorted.vcf.gz \
	--dosages bestguess \
	--out trtools/testsupport/sample_vcfs/annotaTR_vcfs/gangstr_bestguess

python -m trtools.annotaTR.annotaTR \
	--vcf trtools/testsupport/sample_vcfs/dumpSTR_vcfs/trio_chr21_gangstr.sorted.vcf.gz \
	--dosages bestguess_norm \
	--out trtools/testsupport/sample_vcfs/annotaTR_vcfs/gangstr_bestguess_norm

python -m trtools.annotaTR.annotaTR \
	--vcf trtools/testsupport/sample_vcfs/dumpSTR_vcfs/trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr \
	--dosages bestguess_norm \
	--out trtools/testsupport/sample_vcfs/annotaTR_vcfs/hipstr_bestguess_norm

python -m trtools.annotaTR.annotaTR \
	--vcf trtools/testsupport/sample_vcfs/beagle/1kg_snpstr_21_first_100k_second_50_STRs_imputed.vcf.gz \
	--vcftype hipstr \
	--ref-panel trtools/testsupport/sample_vcfs/beagle/1kg_snpstr_21_first_100k_first_50_annotated.vcf.gz \
	--dosages bestguess_norm \
	--out trtools/testsupport/sample_vcfs/annotaTR_vcfs/hipstr_beagle

python -m trtools.annotaTR.annotaTR \
	--vcf trtools/testsupport/sample_vcfs/beagle/beagle_imputed_withap.vcf.gz \
	--vcftype hipstr \
	--ref-panel trtools/testsupport/sample_vcfs/beagle/beagle_refpanel.vcf.gz \
	--match-refpanel-on trimmedalleles \
	--dosages beagleap \
	--out trtools/testsupport/sample_vcfs/annotaTR_vcfs/beagleap_trimmed

# Restrict each output to 200 lines to keep files small
for outfile in gangstr_bestguess gangstr_bestguess_norm hipstr_bestguess_norm hipstr_beagle beagleap_trimmed
do
	cat trtools/testsupport/sample_vcfs/annotaTR_vcfs/${outfile}.vcf | head -n 200 > tmp
	mv tmp trtools/testsupport/sample_vcfs/annotaTR_vcfs/${outfile}.vcf
done