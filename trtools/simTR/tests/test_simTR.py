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
    args.outprefix = str(tmpdir /  "test")
    args.tmpdir = None
    args.coords = None
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 1
    args.coverage = 10
    args.read_length = 100
    args.insert = 300
    args.sd = 100
    args.window = 1000
    args.art = "art_illumina"
    args.single = False
    args.debug = False
    args.seed = 12345
    args.coords = "chr11_CBL:1-2"
    args.repeat_unit = "CGG"
    return args

# Test art_illumina
def check_art():
    if shutil.which("art_illumina") is None:
        common.WARNING("Skipping simTR test. art_illumina not installed")
        return False
    return True

# Test no such file or directory
def test_WrongRefFile(args, simtrdir):
    fname = os.path.join(simtrdir, "test_non_existent.bed")
    if os.path.exists(fname):
        os.remove(fname)
    args.ref = fname
    retcode = main(args)
    assert retcode==1

# Test bad output directory
def test_WrongOutdir(args, simtrdir):
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.outprefix = "bad//x/y/z"
    retcode = main(args)
    assert retcode==1

# Test wrong ART path
def test_WrongARTPath(args, simtrdir):
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.art = "nonexistent_art"
    retcode = main(args)
    assert retcode==1

# Test wrong ART path
def test_WrongARTPath2(args, tmpdir, simtrdir):
    args.art = "ls" # this tool will be found, but it is not ART
    os.mkdir(str(tmpdir /  "test-CBLbad-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBLbad")
    args.tmpdir = str(tmpdir /  "test-CBLbad-tmpdir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 0.01
    args.read_length = 150
    args.insert = 300
    args.sd = 50
    args.coverage = 1000
    args.seed = 12345
    retcode = main(args)
    assert retcode==1

# Test wrong ART path
def test_BadParamCombinations(args, simtrdir):
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.single = True
    retcode = main(args)
    assert retcode==1

# Test params
def test_BadParams(args, simtrdir):
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.u = -1
    retcode = main(args)
    assert retcode == 1

    args.u = 100
    retcode = main(args)
    assert retcode == 1
    args.u = 0.01

    args.d = -5
    retcode = main(args)
    assert retcode == 1

    args.d = 5
    retcode = main(args)
    assert retcode == 1
    args.d = 0.01

    args.rho = -5
    retcode = main(args)
    assert retcode == 1

    args.rho = 5
    retcode = main(args)
    assert retcode == 1
    args.rho = 0.9

    args.p_thresh = -5
    retcode = main(args)
    assert retcode == 1

    args.p_thresh = 5
    retcode = main(args)
    assert retcode == 1
    args.p_thresh = 1

    args.coverage = -1
    retcode = main(args)
    assert retcode == 1
    args.coverage = 100

    args.read_length = -1
    retcode = main(args)
    assert retcode == 1
    args.read_length = 100

    args.insert = -1
    retcode = main(args)
    assert retcode == 1
    args.insert = 100

    args.sd = -1
    retcode = main(args)
    assert retcode == 1
    args.sd = 100

    args.window = -1
    retcode = main(args)
    assert retcode == 1
    args.window = 1000

    args.read_length = 100
    args.insert = 100
    retcode = main(args)
    assert retcode == 1
    args.insert = 300

    args.u = 0.9
    args.d = 0.9
    retcode = main(args)
    assert retcode == 1
    args.u = 0.01
    args.d = 0.01

def test_BadTmpDir(args, tmpdir, simtrdir):
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL1")
    args.tmpdir = str(tmpdir /  "bad-tmp-dir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.read_length = 150
    args.coverage = 100
    args.seed = 12345
    retcode = main(args)
    assert retcode == 1

def test_GoodExampleRun1(args, tmpdir, simtrdir):
    if not check_art(): return
    os.mkdir(str(tmpdir /  "test-CBL1-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL1")
    args.tmpdir = str(tmpdir /  "test-CBL1-tmpdir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 0.01
    args.read_length = 150
    args.insert = 300
    args.sd = 50
    args.coverage = 1000
    args.seed = 12345
    retcode = main(args)
    assert retcode == 0

def test_GoodExampleRun2(args, tmpdir, simtrdir):
    if not check_art(): return
    os.mkdir(str(tmpdir /  "test-CBL2-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL2")
    args.tmpdir = str(tmpdir /  "test-CBL2-tmpdir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 0.01
    args.read_length = 150
    args.coverage = 1000
    args.seed = 12345
    args.single = True
    retcode = main(args)
    assert retcode == 0

def test_GoodExampleRun3(args, tmpdir, simtrdir):
    if not check_art(): return
    args.art = None # make it find ART
    os.mkdir(str(tmpdir /  "test-CBL3-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL3")
    args.tmpdir = str(tmpdir /  "test-CBL3-tmpdir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 0.01
    args.read_length = 150
    args.insert = 300
    args.sd = 50
    args.coverage = 1000
    args.seed = 12345
    retcode = main(args)
    assert retcode == 0

def test_BadCoords(args, tmpdir, simtrdir):
    args.art = None # make it find ART
    os.mkdir(str(tmpdir /  "test-CBL3-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:XXXXXX"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL3")
    args.tmpdir = str(tmpdir /  "test-CBL3-tmpdir")
    args.u = 0.01
    args.d = 0.01
    args.rho = 0.9
    args.p_thresh = 0.01
    args.read_length = 150
    args.insert = 300
    args.sd = 50
    args.coverage = 1000
    args.seed = 12345
    retcode = main(args)
    assert retcode == 1

    args.coords = "chr11_CBL:XXXXXX-YYYY"
    retcode = main(args)
    assert retcode == 1

    args.coords = "chr11_CBL:200-100"
    retcode = main(args)
    assert retcode == 1

    args.coords = "chr11_CBL:5033-5001"
    retcode = main(args)
    assert retcode == 1

    # Bad chrom
    args.coords = "chr11:5000-5033"
    retcode = main(args)
    assert retcode == 1

    # Bad coords
    args.coords = "chr11_CBL:50001-50033"
    retcode = main(args)
    assert retcode == 1

    # Bad repeat unit
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "AT"
    retcode = main(args)
    assert retcode == 1

    # Bad repeat unit
    args.coords = "chr11_CBL:5001-5033"
    args.repeat_unit = "CCG"
    retcode = main(args)
    assert retcode == 1

    # Bad window
    args.window = 10000000
    retcode = main(args)
    assert retcode == 1

    # Bad window
    args.window = 10
    args.insert = 350
    retcode = main(args)
    assert retcode == 1

def test_TooMuchStutter(args, tmpdir, simtrdir):
    args.art = None # make it find ART
    os.mkdir(str(tmpdir /  "test-CBL4-tmpdir"))
    args.ref = os.path.join(simtrdir, "CBL.fa")
    args.coverage = 100
    args.coords = "chr11_CBL:5001-5010"
    args.repeat_unit = "CGG"
    args.outprefix = str(tmpdir /  "test-CBL4")
    args.tmpdir = str(tmpdir /  "test-CBL4-tmpdir")
    args.u = 0.4
    args.d = 0.4
    args.rho = 0.5
    args.p_thresh = 0.01
    args.read_length = 150
    args.insert = 300
    args.sd = 50
    args.coverage = 1000
    args.seed = 12345
    retcode = main(args)
    assert retcode == 1

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

def test_ParseCoordinates4():
    chrom, start, end = ParseCoordinates(0)
    assert chrom is None

    chrom, start, end = ParseCoordinates(":1-100")
    assert chrom is None

    chrom, start, end = ParseCoordinates("xx:-100")
    assert chrom is None

    chrom, start, end = ParseCoordinates("xx:-")
    assert chrom is None    

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

def test_GetMaxDelta4():
    # Rho way too low, delta should be 0
    sprob = 0.02
    rho = 0.01
    pthresh = 0.01
    delta = GetMaxDelta(sprob, rho, pthresh)
    assert delta == 0

def test_GetAlleleSeq1():
    # Keep as is
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = 0
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATATCGTA"
    assert newseq == expected_seq

    # Add a copy
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = 1
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATATATCGTA"
    assert newseq == expected_seq

    # Delete a copy
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = -1
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    expected_seq = "AGCTATATCGTA"
    assert newseq == expected_seq

    # Deletion resulting in a negative sequence length
    seq_preflank = "AGCT"
    seq_postflank = "CGTA"
    seq_repeat = "ATATAT"
    repeat_unit = "AT"
    delta = -5
    newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, repeat_unit, delta)
    assert (newseq is None)

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