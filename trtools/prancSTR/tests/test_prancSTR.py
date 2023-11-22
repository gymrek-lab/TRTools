import argparse
import os

import pytest

from ..prancSTR import *

# Set up base argparser
@pytest.fixture
def args(tmpdir):
    args = argparse.ArgumentParser()
    args.vcf = None
    args.out = str(tmpdir / "test")
    args.region = None
    args.only_passing = False
    args.debug = False
    args.vcftype = "hipstr"
    args.samples = None
    args.quiet = True
    args.output_all = False
    args.readfield = "MALLREADS"
    return args

# Test no such file or directory
def test_WrongFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_non_existent.vcf")
    if os.path.exists(fname):
        os.remove(fname)
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

    # real path but not VCF
    fname = os.path.join(vcfdir, "CEU_test.vcf.gz.tbi")
    args.vcf = fname
    retcode = main(args)
    assert retcode==1

# Test bad output directory
def test_BadOutdir(args, vcfdir, tmpdir):
    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcf = fname
    args.out = str(tmpdir / "bad/test")
    retcode = main(args)
    assert retcode==1

    args.out = str(tmpdir)+os.sep
    print(args.out)
    retcode = main(args)
    assert retcode==1

# Test a good VCF file
def test_RightFile(args, vcfdir):
    fname = os.path.join(vcfdir, "test_hipstr.vcf")
    args.vcftype = "auto"
    args.vcf = fname
    retcode = main(args)
    assert retcode==0

    args.quiet = False
    args.debug = True
    retcode = main(args)
    assert retcode==0    

    # Wrong VCF type
    args.vcftype = "advntr"
    retcode = main(args)
    assert retcode==1

# Test a multi-sample VCF file
def test_MosaicCase(args, vcfdir):
    fname = os.path.join(vcfdir, "CEU_test.vcf.gz")
    args.vcf = fname
    args.quiet = True
    retcode = main(args)
    assert retcode==0

    args.quiet = False
    retcode = main(args)
    assert retcode==0

    # Specific region
    # Note: cyvcf2 handles region parsing and will
    # output a warning if no intervals found
    args.quiet = True
    args.region = "chr1:987287-987288"
    retcode = main(args)
    assert retcode==0

    # Samples list
    args.samples = "NA12878"
    retcode = main(args)
    assert retcode==0

    # Bad samples list should just give no output
    # Note, we do output a warning if we can't 
    # find a sample
    args.samples = "XYZ"
    retcode = main(args)
    assert retcode==0

    # With only passing
    args.samples = "NA12878"
    args.only_passing = True
    args.region = None
    retcode = main(args)
    assert retcode==0

    # With only passing - debug mode
    args.only_passing = True
    args.region = None
    args.debug = True
    retcode = main(args)
    assert retcode==0

    # Write to stdout
    args.samples = "NA12878"
    args.out = "stdout"
    retcode == main(args)
    assert retcode==0

    # With bad readfield
    args.readfield = "badreadfield"
    retcode = main(args)
    assert retcode==1



# Test the probability of observing a certain repeat length
def test_StutterProb1():
    delta = 0
    stutter_u = 0.1
    stutter_d = 0.05
    stutter_rho = 0.2
    expected_prob = 1 - stutter_u - stutter_d
    assert StutterProb(delta, stutter_u, stutter_d, stutter_rho) == expected_prob

def test_StutterProb2():
    delta = 3
    stutter_u = 0.1
    stutter_d = 0.05
    stutter_rho = 0.2
    expected_prob = stutter_u * stutter_rho * (pow((1 - stutter_rho), (delta - 1)))
    assert StutterProb(delta, stutter_u, stutter_d, stutter_rho) == expected_prob

def test_StutterProb3():
    delta = -2
    stutter_u = 0.1
    stutter_d = 0.05
    stutter_rho = 0.2
    expected_prob = stutter_d * stutter_rho * (pow((1 - stutter_rho), (abs(delta) - 1)))
    assert StutterProb(delta, stutter_u, stutter_d, stutter_rho) == expected_prob

def test_StutterProb4():
    delta = 10
    stutter_u = 0.1
    stutter_d = 0.05
    stutter_rho = 0.2
    expected_prob = stutter_u * stutter_rho * (pow((1 - stutter_rho), (delta - 1)))
    assert StutterProb(delta, stutter_u, stutter_d, stutter_rho) == expected_prob

def test_StutterProb5():
    delta = -5
    stutter_u = 0.1
    stutter_d = 0.05
    stutter_rho = 0.2
    expected_prob = stutter_d * stutter_rho * (pow((1 - stutter_rho), (abs(delta) - 1)))
    assert StutterProb(delta, stutter_u, stutter_d, stutter_rho) == expected_prob

#Test the values of C and f
def test_MaximizeMosaicLikelihoodBoth1():
    reads = [10, 11, 10, 11, 10]
    A = 9
    B = 12
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    maxiter = 100
    locname = "None"
    quiet = True
    C, f = MaximizeMosaicLikelihoodBoth(reads, A, B, stutter_probs, maxiter, locname, quiet)
    assert C == 9
    assert f == 0.01

def test_MaximizeMosaicLikelihoodBoth2():
    reads = [-3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2]
    A = -2
    B = -2
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    maxiter = 100
    locname = "None"
    quiet = True
    C, f = MaximizeMosaicLikelihoodBoth(reads, A, B, stutter_probs, maxiter, locname, quiet)
    assert C == -2
    assert f== 0.01

def test_MaximizeMosaicLikelihoodBoth3():
    reads = [-5, -5, -4, -4, -3, -3, -2, -2, -1, -1]
    A = -5
    B = -1
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    maxiter = 100
    locname = "None"
    quiet = True
    C, f = MaximizeMosaicLikelihoodBoth(reads, A, B, stutter_probs, maxiter, locname, quiet)
    assert C == -5
    assert f == pytest.approx(0.0167, abs=1e-2)

# Test values of the alleles into the difference in repetitions with respect to the reference
def test_ExtractReadVector1():
    mallreads=None
    period=3
    reads=ExtractReadVector(mallreads, period)
    assert reads==[]

def test_ExtractReadVector2():
    mallreads="-6|4;-4|28"
    period=1
    reads=ExtractReadVector(mallreads, period)
    assert reads==[-6, -6, -6, -6, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4]

def test_ExtractReadVector3():
    mallreads="9|3;10|5;11|2"
    period=1
    reads=ExtractReadVector(mallreads, period)
    assert reads ==[9, 9, 9, 10, 10, 10, 10, 10, 11, 11]

def test_ExtractReadVector4():
    mallreads="-12|9;-4|16;0|29;4|11"
    period=2
    reads=ExtractReadVector(mallreads, period)
    assert reads ==[-6, -6, -6, -6, -6, -6, -6, -6, -6, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]

# Test minimum and maximum values of a number
def test_ConfineRange1():
    x_cons=ConfineRange(30, 40, 50)
    assert x_cons==40

def test_ConfineRange2():
    x_cons=ConfineRange(60, 40, 50)
    assert x_cons==50

def test_ConfineRange3():
    x_cons=ConfineRange(45, 40, 50)
    assert x_cons==45

#Test the likelihood of observing the reads
def test_Likelihood_mosaic1():
    reads = [10, 11, 10, 11, 10]
    A = 9
    B = 12
    C = 9
    f = 0.01
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    result_likelihood = Likelihood_mosaic(A, B, C, f, reads, stutter_probs)
    assert -2300 <= result_likelihood <= -2290

def test_Likelihood_mosaic2():
    reads = [-3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2]
    A = -2
    B = -2
    C = -2
    f = 0.01
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    result_likelihood = Likelihood_mosaic(A, B, C, f, reads, stutter_probs)
    assert -15000 <= result_likelihood <= -14000 

def test_Likelihood_mosaic3():
    reads = [-5, -5, -4, -4, -3, -3, -2, -2, -1, -1]
    A = -5
    B = -1
    C = -5
    f = pytest.approx(0.0167, 1e-2)
    stutter_probs = [x * 0.001 for x in range(-100, 100)]
    result_likelihood = Likelihood_mosaic(A, B, C, f, reads, stutter_probs)
    assert -4600 <= result_likelihood <= -4550     

#Test the survival function of a point mass at 0
def test_SF1():
    x=10
    result_sf=SF(x)
    assert result_sf==0

def test_SF2():
    x=0
    result_sf=SF(x)
    assert result_sf==1

def test_SF3():
    x=-1
    result_sf=SF(x)
    assert result_sf==1

#Test the C prediction
def test_Just_C_Pred1():
    reads = [10, 11, 10, 11, 10]
    A = 9
    B = 12
    f = 0.01
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    C = Just_C_Pred(reads, A, B, f, stutter_probs)
    assert C == 9

def test_Just_C_Pred2():    
    reads = [-6, -6, -6, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,]
    A = -2
    B = -2
    f = 0.0362320
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    C = Just_C_Pred(reads, A, B, f, stutter_probs)
    assert C == -2

def test_Just_C_Pred3():    
    reads = [-5, -5, -4, -4, -3, -3, -2, -2, -1, -1]
    A = -5
    B = -1
    f = 0.0167
    stutter_probs = [x * 0.001 for x in range(-200, 201)]
    C = Just_C_Pred(reads, A, B, f, stutter_probs)
    assert C == -5

#Test the f prediction
def test_Just_F_Pred1():
    reads = [10, 11, 10, 11, 10]
    A = 9
    B = 12
    C = 9
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    f = Just_F_Pred(reads, A, B, C, stutter_probs)
    assert f == 0.01

def test_Just_F_Pred2():
    reads = [-6, -6, -6, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,]
    A = -2
    B = -2
    C = -2 
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    f = Just_F_Pred(reads, A, B, C, stutter_probs)
    assert f == pytest.approx(0.036, abs=1e-1)

def test_Just_F_Pred3():
    reads = [-5, -5, -4, -4, -3, -3, -2, -2, -1, -1]
    A = -5
    B = -1
    C = -5 
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    f = Just_F_Pred(reads, A, B, C, stutter_probs)
    assert f == pytest.approx(0.0167, abs=1e-2)

#
def test_ComputePvalue1():
    reads = [10, 11, 10, 11, 10]
    A = 9
    B = 12
    best_C = 9
    best_f = 0.01
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    pval = ComputePvalue(reads, A, B, best_C, best_f, stutter_probs)
    assert pval == 1

def test_ComputePvalue2():
    reads = [-6, -6, -6, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,]
    A = -2
    B = -2
    best_C = -2
    best_f = 0.0362320
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    pval = ComputePvalue(reads, A, B, best_C, best_f, stutter_probs)
    assert pval == 1

def test_ComputePvalue3():
    reads = [-3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2]
    A = -5
    B = -1
    best_C = -5
    best_f = 0.0167
    stutter_probs = [x * 0.001 for x in range(-100, 101)]
    pval = ComputePvalue(reads, A, B, best_C, best_f, stutter_probs)
    assert pval == 1