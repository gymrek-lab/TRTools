#!/bin/bash

# This script contains command line tests for TRTools utilities

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

runcmd_pass()
{
    echo "[runcmd_pass]: $1"
    sh -c "$1" >/dev/null 2>&1 || die "Error running: $1"
}

runcmd_fail()
{
    echo "[runcmd_fail]: $1"
    sh -c "$1" >/dev/null 2>&1 && die "Command should have failed: $1"
}

EXDATADIR="example-files"
TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

# Check version
for tool in mergeSTR dumpSTR qcSTR statSTR compareSTR
do
    runcmd_pass "${tool} --version"
done

runcmd_pass "python -c 'import trtools; print(trtools.__version__)'"

runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test --mean"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/ --mean"
runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR} --mean"
#runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/ --mean"

runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/"
#runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
#runcmd_fail "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
#runcmd_fail "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
#runcmd_fail "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
#runcmd_fail "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

# TODO check with bcftools index for mergestr, comparestr

runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype hipstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype eh"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype advntr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out stdout --mean --vcftype popstr"

runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype gangstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype eh"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype advntr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz --out stdout --mean --vcftype popstr"

runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype gangstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype hipstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype advntr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out stdout --mean --vcftype popstr"

runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype gangstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype hipstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype advntr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz --out stdout --mean --vcftype eh"

runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype gangstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype hipstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype popstr"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --mean --vcftype eh"

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
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_advntr --vcftype advntr --update-sample-from-file"
runcmd_fail "mergeSTR --vcfs ${FILE1},${FILE1} --out ${TMPDIR}/test_merge_advntr --vcftype advntr" # duplicate samples

# ExpansionHunter
# Note, you first need to reheader files to add required contig lines to VCF headers
for sample in NA12878 NA12891 NA12892; do
    bcftools reheader -f ${EXDATADIR}/hg19.fa.fai -o ${TMPDIR}/${sample}_eh_reheader.vcf.gz ${EXDATADIR}/${sample}_chr21_eh.sorted.vcf.gz >/dev/null 2>&1 || die "bcftools failed"
    tabix -p vcf ${TMPDIR}/${sample}_eh_reheader.vcf.gz >/dev/null 2>&1 || die "tabix failed"
done
FILE1=${TMPDIR}/NA12878_eh_reheader.vcf.gz
FILE2=${TMPDIR}/NA12891_eh_reheader.vcf.gz
FILE3=${TMPDIR}/NA12892_eh_reheader.vcf.gz
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_eh --vcftype eh"

# GangSTR
FILE1=${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_gangstr.sorted.vcf.gz
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_gangstr --vcftype gangstr"

# HipSTR
FILE1=${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_hipstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_hipstr.sorted.vcf.gz
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_hipstr --vcftype hipstr"

# PopSTR
FILE1=${EXDATADIR}/NA12878_chr21_popstr.sorted.vcf.gz
FILE2=${EXDATADIR}/NA12891_chr21_popstr.sorted.vcf.gz
FILE3=${EXDATADIR}/NA12892_chr21_popstr.sorted.vcf.gz
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_popstr --vcftype popstr"

runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out stdout --afreq"
runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12891_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/stats_eh --numcalled"
runcmd_pass "statSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/stats_gangstr --numcalled --mean"
runcmd_pass "statSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr --out ${TMPDIR}/stats_gangstr --acount --afreq --mean"
runcmd_pass "statSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/stats_popstr --mean --samples ${EXDATADIR}/ex-samples.txt"

runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --advntr-min-call-DP 5 --out ${TMPDIR}/test_dumpstr_advntr"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_eh --eh-min-call-LC 50 --num-records 10 --drop-filtered"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_gangstr --min-locus-callrate 0.9 --num-records 10"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --vcftype hipstr --out ${TMPDIR}/test_dumpstr_hipstr --filter-hrun --num-records 10"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/test_dumpstr_popstr --min-locus-callrate 0.9 --popstr-min-call-DP 10 --num-records 100"

FILE1=${TMPDIR}/NA12878_advntr_reheader.vcf.gz
runcmd_pass "compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out ${TMPDIR}/advntr_vs_advntr --noplot"
runcmd_pass "compareSTR \
    --vcf1 ${EXDATADIR}/NA12878_chr21_hipstr.sorted.vcf.gz \
    --vcf2 ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz \
    --vcftype1 hipstr --vcftype2 eh --out ${TMPDIR}/hipstr_vs_eh"

FILE1=${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz
runcmd_pass "compareSTR --vcf1 ${FILE1} --vcf2 ${FILE1} --out ${TMPDIR}/popstr_vs_popstr"

runcmd_pass "qcSTR --vcf ${EXDATADIR}/trio_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_gangstr --period 4 --quality per-locus"
runcmd_pass "qcSTR --vcf ${EXDATADIR}/trio_chr21_hipstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_hipstr --vcftype hipstr --samples ${EXDATADIR}/ex-samples.txt"
runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_eh.sorted.vcf.gz --out ${TMPDIR}/test_qc_eh"
runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --out ${TMPDIR}/test_qc_advntr"
runcmd_pass "qcSTR --vcf ${EXDATADIR}/trio_chr21_popstr.sorted.vcf.gz --out ${TMPDIR}/test_qc_popstr"

runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_gangstr.vcf --out ${TMPDIR}/test_qc_gangstr --period 4 --quality per-locus"
runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_hipstr.vcf --out ${TMPDIR}/test_qc_hipstr --vcftype hipstr --samples ${EXDATADIR}/ex-samples.txt"
runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_eh.vcf --out ${TMPDIR}/test_qc_eh"
#runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_advntr.vcf --out ${TMPDIR}/test_qc_advntr"
#runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_popstr.vcf --out ${TMPDIR}/test_qc_popstr"

echo "tests completed successfully!"
exit 0

