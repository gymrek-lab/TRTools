#!/usr/bin/env bash

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
    bash -c "$1" >/dev/null 2>&1 || die "Error running: $1"
}

runcmd_fail()
{
    echo "[runcmd_fail]: $1"
    bash -c "$1" >/dev/null 2>&1 && die "Command should have failed: $1"
}

if [ $# -eq 0 ]; then
    # use default example location
    EXDATADIR="example-files"
    BEAGLEDIR="trtools/testsupport/sample_vcfs/beagle"
elif (( $# != 2 )) ; then
    echo "usage: cmdline_tests.sh {example_dir} {beagle_dir}" 2>&1
    echo "Expected 2 arguments but recieved $#" 2>&1
    exit 1
else
    EXDATADIR=$1
    BEAGLEDIR=$2
fi

TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

# Check version
for tool in mergeSTR dumpSTR qcSTR statSTR compareSTR associaTR prancSTR simTR
do
    runcmd_pass "${tool} --version"
done

runcmd_pass "python -c 'import trtools; print(trtools.__version__)'"

# Check for valid/invalid output locations

# Example command running prancSTR for only one chromosome with hipstr output file
# --only-passing skips VCF records with non-passing filters
runcmd_pass "prancSTR --vcf ${EXDATADIR}/CEU_subset.vcf.gz --out ${TMPDIR}/CEU_chr1 --vcftype hipstr --only-passing --region chr1"
# Example command running prancSTR for only one sample
runcmd_pass "prancSTR --vcf ${EXDATADIR}/CEU_subset.vcf.gz --only-passing --out ${TMPDIR}/NA12878_chr1 --samples NA12878"

if ! command -v art_illumina &> /dev/null; then
    echo "Skipping simTR tests. art_illumina not found"
else
    # Example command running simTR for a dummy dataset with dummy allele bed file and other input parameters
    mkdir ${TMPDIR}/test-simtr
    runcmd_pass "simTR --coords chr11_CBL:5001-5033 --ref ${EXDATADIR}/CBL.fa --outprefix ${TMPDIR}/test-simtr --tmpdir ${TMPDIR}/test-simtr --repeat-unit CGG --art art_illumina --coverage 1000 --read-length 150 --seed 12345 --u 0.02 --d 0.02 --rho 0.9"
fi

runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test --mean"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx --mean"
runcmd_pass "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR} --mean"
runcmd_fail "statSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/ --mean"

runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
runcmd_fail "qcSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
runcmd_fail "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
runcmd_fail "mergeSTR --vcfs ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz,${EXDATADIR}/NA12891_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/test"
runcmd_fail "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/kittens/xxx"
runcmd_pass "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}"
runcmd_fail "compareSTR --vcf1 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --vcf2 ${EXDATADIR}/NA12878_chr21_gangstr.sorted.vcf.gz --out ${TMPDIR}/"

runcmd_pass "associaTR association_results.tsv ${EXDATADIR}/ceu_ex.vcf.gz simulated_phenotype ${EXDATADIR}/simulated_traits_0.npy --same-samples"
runcmd_pass "associaTR association_results.tsv ${EXDATADIR}/ceu_ex.vcf.gz simulated_phenotype ${EXDATADIR}/simulated_traits_0.npy ${EXDATADIR}/simulated_traits_1.npy --same-samples"
runcmd_fail "associaTR association_results.tsv nonexistant simulated_phenotype ${EXDATADIR}/simulated_traits_0.npy ${EXDATADIR}/simulated_traits_1.npy --same-samples"
runcmd_fail "associaTR association_results.tsv ${EXDATADIR}/ceu_ex.vcf.gz simulated_phenotype nonexistant --same-samples"
runcmd_fail "associaTR association_results.tsv ${EXDATADIR}/ceu_ex.vcf.gz simulated_phenotype ${EXDATADIR}/simulated_traits_0.npy nonexistant --same-samples"

# check for invalid vcftypes

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
runcmd_fail "mergeSTR --vcfs ${FILE1},${FILE1} --out ${TMPDIR}/test_merge_advntr_dup --vcftype advntr" # duplicate samples

# ExpansionHunter
# Note, you first need to reheader files to add required contig lines to VCF headers
for sample in NA12878 NA12891 NA12892; do
    bcftools reheader -f ${EXDATADIR}/hg19.fa.fai -o ${TMPDIR}/${sample}_eh_reheader.vcf.gz ${EXDATADIR}/${sample}_chr21_eh.sorted.vcf.gz >/dev/null 2>&1 || die "bcftools failed"
    tabix -p vcf ${TMPDIR}/${sample}_eh_reheader.vcf.gz >/dev/null 2>&1 || die "tabix failed"

    # Make bcftools index'ed versions
    cp ${TMPDIR}/${sample}_eh_reheader.vcf.gz ${TMPDIR}/${sample}_eh_bcf_reheader.vcf.gz
    bcftools index ${TMPDIR}/${sample}_eh_bcf_reheader.vcf.gz
done
FILE1=${TMPDIR}/NA12878_eh_reheader.vcf.gz
FILE2=${TMPDIR}/NA12891_eh_reheader.vcf.gz
FILE3=${TMPDIR}/NA12892_eh_reheader.vcf.gz
runcmd_pass "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_eh --vcftype eh"

FILE1=${TMPDIR}/NA12878_eh_bcf_reheader.vcf.gz
FILE2=${TMPDIR}/NA12891_eh_bcf_reheader.vcf.gz
FILE3=${TMPDIR}/NA12892_eh_bcf_reheader.vcf.gz
runcmd_fail "mergeSTR --vcfs ${FILE1},${FILE2},${FILE3} --out ${TMPDIR}/test_merge_eh --vcftype eh" # should fail with BCF index

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

runcmd_pass "dumpSTR --vcf ${EXDATADIR}/NA12878_chr21_advntr.sorted.vcf.gz --advntr-min-call-DP 100 --out ${TMPDIR}/test_dumpstr_advntr"
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
runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_advntr.vcf --out ${TMPDIR}/test_qc_advntr"
runcmd_pass "qcSTR --vcf ${TMPDIR}/test_merge_popstr.vcf --out ${TMPDIR}/test_qc_popstr"

echo "--- Running trtools_prep_beagle_vcf.sh tests --- "
prep_beagle_out="$TMPDIR"/test_prep_beagle_vcf.vcf.gz
ref_panel="$BEAGLEDIR"/1kg_snpstr_21_first_100k_first_50_annotated.vcf.gz
imputed_vcf="$BEAGLEDIR"/1kg_snpstr_21_first_100k_second_50_STRs_imputed.vcf.gz

runcmd_fail "trtools_prep_beagle_vcf.sh hipstr nonexistent.vcf.gz $imputed_vcf $prep_beagle_out"
runcmd_fail "trtools_prep_beagle_vcf.sh hipstr $ref_panel nonexistent.vcf.gz $prep_beagle_out"

trtools_prep_beagle_vcf.sh hipstr "$ref_panel" "$imputed_vcf" "$prep_beagle_out"

if ! [[ -f "$prep_beagle_out" ]] ; then
    echo "prep_beagle_vcf test didn't produce output file" >&2
    exit 1
fi

if ! [[ -f "$prep_beagle_out".tbi ]] ; then
    echo "prep_beagle_vcf test didn't produce index file" >&2
    exit 1
fi

if (( 1172 != $(zcat < "$prep_beagle_out" | grep -vc '#') )) ; then
    echo "prep_beagle_vcf outputted a file that didn't have the expected number of lines (1172)"
    exit 1
fi

if (( 1172 != $(zcat < "$prep_beagle_out" | grep -v '#' | grep -c 'START') )) ||
    (( 1172 != $(zcat < "$prep_beagle_out" | grep -v '#' | grep -c 'END') )) ||
    (( 1172 != $(zcat < "$prep_beagle_out" | grep -v '#' | grep -c 'PERIOD') ))
then
    echo "prep_beagle_vcf outputted a file that didn't have the expected number of INFO annotations"
    exit 1
fi
echo '------'

echo "tests completed successfully!"
exit 0

