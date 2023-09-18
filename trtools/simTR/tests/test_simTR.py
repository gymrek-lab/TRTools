import argparse
import os

import pytest

from ..simTR import *

# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.ref = None
    args.repeat_unit = None
    args.outprefix = None
    args.tmpdir = None
    args.u = None
    args.d = None
    args.rho = None
    args.p_thresh = None
    args.coverage = None
    args.read_length = None
    args.insert = None
    args.sd = None
    args.window = None
    args.art = None
    args.single = False
    args.debug = False
    return args

#Test the Coordinates
def test_ParseCoordinates1():
    chrom, start, end = ParseCoordinates("chr1:200-500")
    assert chrom == "chr1"
    assert start == 200
    assert end == 500

def test_ParseCoordinates2():
    chrom, start, end = ParseCoordinates("chrX:1000-1500")
    assert chrom == "chrX"
    assert start == 1000
    assert end == 1500

def test_ParseCoordinates3():
    chrom, start, end = ParseCoordinates("chrY:300-600")
    assert chrom == "chrY"
    assert start == 300
    assert end == 600

#Test the Max Delta
def test_GetMaxDelta1():
    sprob = 0.05
    rho = 0.9
    pthresh = 0.001
    delta = GetMaxDelta(sprob, rho, pthresh)
    assert delta == 3

def test_GetMaxDelta2():
    sprob = 0.1
    rho = 0.8
    pthresh = 0.0001
    delta = GetMaxDelta(sprob, rho, pthresh)
    assert delta == 6

def test_GetMaxDelta3():
    sprob = 0.02
    rho = 0.95
    pthresh = 0.00001
    delta = GetMaxDelta(sprob, rho, pthresh)
    assert delta == 4

#
def test_GetAlleleSeq1():
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = 0
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATATCGTA"
    assert newseq == expected_seq

def test_GetAlleleSeq2():
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = 2
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATATATATCGTA"
    assert newseq == expected_seq

def test_GetAlleleSeq3():
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = -1
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATCGTA"
    assert newseq == expected_seq

#
def test_CreateAlleleFasta1(tmpdir):
    newseq = "AGCTATATATCGTA"
    delta = 3
    fasta_path = CreateAlleleFasta(newseq, delta, tmpdir)
    assert os.path.exists(fasta_path)
    assert os.path.isfile(fasta_path)
    with open(fasta_path, "r") as f:
        lines = f.readlines()
        assert len(lines) == 2
        assert lines[0] == ">seq_3\n"
        assert lines[1] == "AGCTATATATCGTA\n"

def test_CreateAlleleFasta2(tmpdir):
    newseq = "ACGT"
    delta = 0
    fasta_path = CreateAlleleFasta(newseq, delta, tmpdir)
    assert os.path.exists(fasta_path)
    assert os.path.isfile(fasta_path)
    with open(fasta_path, "r") as f:
        lines = f.readlines()
        assert len(lines) == 2
        assert lines[0] == ">seq_0\n"
        assert lines[1] == "ACGT\n"

def test_CreateAlleleFasta3(tmpdir):
    newseq = "TGCATG"
    delta = -2
    fasta_path = CreateAlleleFasta(newseq, delta, tmpdir)
    assert os.path.exists(fasta_path)
    assert os.path.isfile(fasta_path)
    with open(fasta_path, "r") as f:
        lines = f.readlines()
        assert len(lines) == 2
        assert lines[0] == ">seq_-2\n"
        assert lines[1] == "TGCATG\n"