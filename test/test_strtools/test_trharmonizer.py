import os, sys
import pytest
import vcf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','strtools'))

import tr_harmonizer as trh

COMMDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "common")
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")

#### Test TRRecord using dummy info ####
# Set up dummy class with gt_alleles
class DummyVCFSample:
    def __init__(self, gt_alleles):
        self.gt_alleles = gt_alleles
        self.called = True
# Set up dummy VCF records which are just lists of genotypes
dummy_record1 = [] # Example record with real data
dummy_record1.append(DummyVCFSample(['0','1']))
dummy_record1.append(DummyVCFSample(['1','1']))
dummy_record1.append(DummyVCFSample(['1','1']))
dummy_record1.append(DummyVCFSample(['1','2']))
dummy_record1.append(DummyVCFSample(['2','2']))
dummy_record1.append(DummyVCFSample(['0'])) # add a haploid sample
dummy_record2 = [] # Empty record
dummy_record3 = [] # All reference
for i in range(3): dummy_record3.append(DummyVCFSample(['0','0']))
dummy_record3.append(DummyVCFSample(['0','0','0'])) # add a triploid sample

def test_TRRecord_GetGenotypes():
    # Test good example
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec) # To test str function
    true_gts = [[ref_allele, alt_alleles[0]], \
                [alt_alleles[0], alt_alleles[0]], \
                [alt_alleles[0], alt_alleles[0]], \
                [alt_alleles[0], alt_alleles[1]], \
                [alt_alleles[1], alt_alleles[1]], [ref_allele]]
    true_len_gts = [[3, 4], [4, 4], [4, 4], [4, 6], [6, 6], [3]]
    ind = 0
    for sample in rec.vcfrecord:
        stringt = rec.GetStringGenotype(sample)
        lengt = rec.GetLengthGenotype(sample)
        assert(all([(stringt[i]==true_gts[ind][i]) for i in range(len(stringt))]))
        assert(all([(lengt[i]==true_len_gts[ind][i]) for i in range(len(lengt))]))
        ind += 1
    # Test example where alt=[None]
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    print(rec) # To test str function
    for sample in rec.vcfrecord:
        stringt = rec.GetStringGenotype(sample)
        lengt = rec.GetLengthGenotype(sample)
        assert(all([item==ref_allele for item in stringt]))
        assert(all([item==3 for item in lengt]))
    # Test example with discrepancy between alt_alleles and genotypes given
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, [None], "CAG", "")
    
def test_GetGenotypeCounts():
    # TODO Test using dummy records above.
    # Test samplelist, uselength options
    pass 

def test_GetAlleleCounts():
    # TODO Test using dummy records above.
    # Test samplelist, uselength options
    pass

def test_GetAlleleFreqs():
    # TODO Test using dummy records above.
    # Test samplelist, uselength options
    pass

def test_GetMaxAllele():
    # TODO Test using dummy records above.
    # Test samplelist option
    pass

#### Test TRRecordHarmonizer on different files ####
# TODO: test that we can correctly infer the vcf type for gangstr, advntr, hipstr, eh, popstr
# TODO: test that we don't explode if a user gives the wrong input type, but rather fail gracefully (e.g. eh for gangstr)
# TODO: test input an invalid VCF type, e.g. vcftype="notarealformat"
# TODO: test input a VCF that came from something else e.g. SNP calls from samtools
# TODO: test the actual harmonizer output on one or more records from each type
