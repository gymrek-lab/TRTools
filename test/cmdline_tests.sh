#!/bin/bash

# This script contains command line tests for TRTools utilities

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

EXDATADIR="example-files"
TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

echo "** Checking version **"
# Check version
for tool in mergeSTR dumpSTR qcSTR statSTR compareSTR
do
    sh -c "${tool} --version" >/dev/null 2>&1 || die "Failed to get ${tool} version"
done
python -c "import trtools; print(trtools.__version__)" >/dev/null 2>&1 || die "Failed to get version from TRTools library"

echo "** Checking outprefix options**"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test --mean  || die "Should write statstr to test.tab"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/ --mean >/dev/null 2>&1 && die "Trying to set outprefix to nonexistent dir"
# TODO uncomment these
#statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR} --mean >/dev/null 2>&1 && die "Trying to set outprefix to dirname"
#statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/ --mean >/dev/null 2>&1 && die "Trying to set outprefix to dirname"
qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test >/dev/null 2>&1 || die "Should write qc files to test prefix"
qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx >/dev/null 2>&1 && die "Trying to set outputrefix to nonexistent directory"
# TODO uncomment these
#qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR} >/dev/null 2>&1 && die "Trying to set outputrefix to dirname"
#qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/ >/dev/null 2>&1 && die "Trying to set outputrefix to dirname"
dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test >/dev/null 2>&1 || die "Should write qc files to test prefix"
dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx >/dev/null 2>&1 && die "Trying to set outputrefix to nonexistent directory"
# TODO uncomment these
#dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR} >/dev/null 2>&1 && die "Trying to set outputrefix to dirname"
#dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/ >/dev/null 2>&1 && die "Trying to set outputrefix to dirname"

# TODO add outprefix tests for mergeSTR, compareSTR

# TODO check with bcftools index

echo "** Checking setting --vcftype incorrectly **"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype hipstr >/dev/null 2>&1 && die "Should be gangstr VCF, hipstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype eh >/dev/null 2>&1 && die "Should be gangstr VCF, eh specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype advntr >/dev/null 2>&1 && die "Should be gangstr VCF, advntr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype popstr >/dev/null 2>&1 && die "Should be gangstr VCF, popstr specified"

statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype gangstr >/dev/null 2>&1 && die "Should be hipstr VCF, gangstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype eh >/dev/null 2>&1  && die "Should be hipstr VCF, eh specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype advntr >/dev/null 2>&1 && die "Should be hipstr VCF, advntr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype popstr >/dev/null 2>&1 && die "Should be hipstr VCF, popstr specified"

statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype gangstr >/dev/null 2>&1 && die "Should be EH VCF, gangstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype hipstr >/dev/null 2>&1 && die "Should be EH VCF, hipstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype advntr >/dev/null 2>&1 && die "Should be EH VCF, advntr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype popstr >/dev/null 2>&1 && die "Should be EH VCF, popstr specified"

statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype gangstr >/dev/null 2>&1 && die "Should be popstr VCF, gangstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype hipstr >/dev/null 2>&1 && die "Should be popstr VCF, hipstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype advntr >/dev/null 2>&1 && die "Should be popstr VCF, advntr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype eh >/dev/null 2>&1 && die "Should be popstr VCF, EH specified"

statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype gangstr >/dev/null 2>&1 && die "Should be advntr VCF, gangstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype hipstr >/dev/null 2>&1 && die "Should be advntr VCF, hipstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype popstr >/dev/null 2>&1 && die "Should be advntr VCF, popstr specified"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype eh >/dev/null 2>&1 && die "Should be advntr VCF, EH specified"

echo "** Checking merge **"
# Test mergeSTR on all supported tools
# AdVNTR
# Note, you first need to reheader files to add required contig lines to VCF headers
for sample in NA12878 NA12891 NA12892; do
    bcftools reheader -f ${EXDATADIR}/hg19.fa.fai -o ${TMPDIR}/${sample}_advntr_reheader.vcf.gz ${EXDATADIR}/${sample}_chr21_advntr.sorted.vcf.gz >/dev/null 2>&1 || die "bcftools failed"
    tabix -p vcf ${TMPDIR}/${sample}_advntr_reheader.vcf.gz >/dev/null 2>&1 || die "tabix failed"
done
FILE1=${TMPDIR}/NA12878_advntr_reheader.vcf.gz
FILE2=${TMPDIR}/NA12891_advntr_reheader.vcf.gz
FILE3=${TMPDIR}/NA12892_advntr_reheader.vcf.gz
mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_advntr --vcftype advntr --update-sample-from-file >/dev/null 2>&1 || die "Merge adVNTR failed"

# ExpansionHunter
# Note, you first need to reheader files to add required contig lines to VCF headers
for sample in NA12878 NA12891 NA12892; do
    bcftools reheader -f ${EXDATADIR}/hg19.fa.fai -o ${TMPDIR}/${sample}_eh_reheader.vcf.gz ${EXDATADIR}/${sample}_chr21_eh.sorted.vcf.gz >/dev/null 2>&1 || die "bcftools failed"
    tabix -p vcf ${TMPDIR}/${sample}_eh_reheader.vcf.gz >/dev/null 2>&1 || die "tabix failed"
done
FILE1=${TMPDIR}/NA12878_eh_reheader.vcf.gz
FILE2=${TMPDIR}/NA12891_eh_reheader.vcf.gz
FILE3=${TMPDIR}/NA12892_eh_reheader.vcf.gz
mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_eh --vcftype eh >/dev/null 2>&1 || die "merge EH failed"

# GangSTR
FILE1=${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_gangstr.sorted.vcf.gz
mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_gangstr --vcftype gangstr >/dev/null 2>&1 || die "merge gangstr failed"

# HipSTR
FILE1=${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_hipstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_hipstr.sorted.vcf.gz
mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_hipstr --vcftype hipstr >/dev/null 2>&1 || die "merge hipstr failed"

# PopSTR
FILE1=${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_popstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_popstr.sorted.vcf.gz
mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_popstr --vcftype popstr >/dev/null 2>&1 || die "merge popstr failed"

echo "** Checking statstr **"
statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --afreq >/dev/null 2>&1 || die "statstr advntr failed"
statSTR --vcf ${EXDATADIR}/NA12891_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/stats_eh --numcalled >/dev/null 2>&1 || die "statstr eh failed"
statSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/stats_gangstr --numcalled --mean >/dev/null 2>&1 || die "statsr gangstr failed"
statSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr --out ${TMPDIR}/stats_gangstr --acount --afreq --mean >/dev/null 2>&1 || die "statstr hipstr failed"
statSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/stats_popstr --mean --samples ${EXDATADIR}/ex-samples.txt >/dev/null 2>&1 || die "statstr popstr failed"

echo "** Checking dumpSTR **"
dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --advntr-min-call-DP 5 --out ${TMPDIR}/test_dumpstr_advntr >/dev/null 2>&1 || die "dumpstr advntr failed"
dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_eh --eh-min-call-LC 50 --num-records 10 --drop-filtered >/dev/null 2>&1 || die "dumpstr EH failed"
dumpSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_gangstr --min-locus-callrate 0.9 --num-records 10 >/dev/null 2>&1 || die "dumpstr gangstr failed"
dumpSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr --out ${TMPDIR}/test_dumpstr_hipstr --filter-hrun --num-records 10 >/dev/null 2>&1 || die "dumpstr hipstr failed"
dumpSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_popstr --min-locus-callrate 0.9 --popstr-min-call-DP 10 --num-records 100 >/dev/null 2>&1 || die "dumpstr popstr failed"

echo "** Checking compareSTR **"
FILE1=${TMPDIR}/NA12878_advntr_reheader.vcf.gz
compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out ${TMPDIR}/advntr_vs_advntr --noplot >/dev/null 2>&1 || die "compareSTR advntr failed"
compareSTR \
    --vcf1 ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz \
    --vcf2 ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz \
    --vcftype1 hipstr --vcftype2 eh --out ${TMPDIR}/hipstr_vs_eh >/dev/null 2>&1 || die "compareSTR hipstr eh failed"

FILE1=${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz
compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out ${TMPDIR}/popstr_vs_popstr >/dev/null 2>&1 || die "comparestr popstr failed"

echo "** Checking qcstr **"
qcSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_gangstr --period 4 --quality per-locus >/dev/null 2>&1 || die "QC gangstr failed"
qcSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_hipstr --vcftype hipstr --samples ${EXDATADIR}/ex-samples.txt >/dev/null 2>&1 || die "QC hipstr failed"
qcSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/test_qc_eh >/dev/null 2>&1 || die "QC EH failed"
qcSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out ${TMPDIR}/test_qc_advntr >/dev/null 2>&1 || die "QC advntr failed"
qcSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_popstr >/dev/null 2>&1 || die "QC popstr failed"

echo "** Checking qcstr - merged VCFs**"
qcSTR --vcf ${TMPDIR}/test_merge_gangstr.vcf --out ${TMPDIR}/test_qc_gangstr --period 4 --quality per-locus >/dev/null 2>&1 || die "QC gangstr merged failed"
qcSTR --vcf ${TMPDIR}/test_merge_hipstr.vcf --out ${TMPDIR}/test_qc_hipstr --vcftype hipstr --samples ${EXDATADIR}/ex-samples.txt >/dev/null 2>&1 || die "QC hipstr merged failed"
qcSTR --vcf ${TMPDIR}/test_merge_eh.vcf --out ${TMPDIR}/test_qc_eh >/dev/null 2>&1 || die "QC EH merged failed"
#qcSTR --vcf ${TMPDIR}/test_merge_advntr.vcf --out ${TMPDIR}/test_qc_advntr >/dev/null 2>&1 || die "QC advntr merged failed" # TODO uncomment
#qcSTR --vcf ${TMPDIR}/test_merge_popstr.vcf --out ${TMPDIR}/test_qc_popstr >/dev/null 2>&1 || die "QC popstr merged failed" # TODO uncomment

echo "tests completed successfully!"
exit 0


