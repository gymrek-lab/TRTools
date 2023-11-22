#!/usr/bin/env bash

# Run the tests for an installed copy of trtools
# If not already present, megabytes of test data will be downloaded before test running

echo "Running test_trtools.sh"

command -v git >/dev/null 2>&1 || { echo >&2 "git is not available, but is required for downloading the test data. Aborting."; exit 1; }
command -v pytest >/dev/null 2>&1 || { echo >&2 "pytest is not available, but is required for running the unit test suite. Aborting."; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo >&2 "tabix is not available, but is required for running the command line test suite. Aborting."; exit 1; }
command -v bcftools >/dev/null 2>&1 || { echo >&2 "bcftools is not available, but is required for running the command line test suite. Aborting."; exit 1; }

TMP=/tmp/trtools_data_download
if [ ! -d "$TMP" ] ; then
	echo "Downloading test data ..."
	mkdir $TMP
	pushd $TMP
	git init .
	git remote add origin -f https://github.com/gymrek-lab/TRTools.git
	git pull origin master
	popd
	echo "Download done"
else
	echo "Test data already downloaded"
fi

echo "Repo with test data located at $TMP"

# Figure out where trtools is installed
# Note: cd to an arbitrary directory before doing this so that in case
# the current local directory is a copy of the trtools repo,
# python reports the location of the trtools installation and not the
# current location
# (python would find the current one because the local directory is always added 
# to the python path on startup, and is on the path before other installations)
mkdir -p /tmp/trtools_tmp
cd /tmp/trtools_tmp || exit 1
loc=$(dirname "$(python -c 'import trtools;print(trtools.__file__)')")
echo "Location of trtools installation: $loc"

echo "Running pytest ..."
cd "$loc" || exit 1
# run unit tests
python -m pytest . -p trtools.testsupport.dataloader --datadir "$TMP"/trtools/testsupport
# run command line tests
$TMP/test/cmdline_tests.sh $TMP/example-files $TMP/trtools/testsupport/sample_vcfs/beagle