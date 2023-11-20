#!/bin/bash

# This script contains command line tests for simTR
# It is separate from cmdline_tests.sh since we
# only want to run if art_illumina is installed,
# but don't want to require that since it is only
# required to run simTR

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

TMPDIR=$(mktemp -d -t tmp-XXXXXXXXXX)

echo "Saving tmp files in ${TMPDIR}"

# Check version
for tool in simTR
do
    runcmd_pass "${tool} --version"
done

runcmd_pass "python -c 'import trtools; print(trtools.__version__)'"

