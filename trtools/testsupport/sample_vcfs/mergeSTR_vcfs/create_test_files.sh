#!/usr/bin/env bash

# This command should be run from the root of the repo

for type in advntr eh gangstr hipstr popstr ; do
	python -m trtools.mergeSTR.mergeSTR \
		--vcfs trtools/testsupport/sample_vcfs/mergeSTR_vcfs/test_file_${type}1.vcf.gz,trtools/testsupport/sample_vcfs/mergeSTR_vcfs/test_file_${type}2.vcf.gz \
		--vcftype ${type} \
		--out trtools/testsupport/sample_vcfs/mergeSTR_vcfs/${type}_merged
done
