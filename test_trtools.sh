#!/usr/bin/env bash

command -v git >/dev/null 2>&1 || { echo >&2 "git is not available, but is required for downloading the test data. Aborting."; exit 1; }
command -v pytest >/dev/null 2>&1 || { echo >&2 "pytest is not available, but is required for running the test suite. Aborting."; exit 1; }

containing_dir=$(dirname $0)
if [ ! -d "$containing_dir/tests/common/sample_regions" ]  ; then
	echo "Downloading test data ..."
	TMP=/tmp/trtools_data_download
	mkdir $TMP
	pushd $TMP
    git init .
	git remote add origin -f https://github.com/gymreklab/TRTools.git
	# Inform git to only download the files we're interested in
	# (this should work for newer versions of git, though not tested)
	# won't work for older ones
	printf "tests/common/sample_regions\ntests/common/sample_vcfs\n" > \
		.git/info/sparse-checkout
	git pull origin master
	popd
	echo "Download done"
fi

pytest . --datadir $TMP/tests/common

