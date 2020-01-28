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
    def __init__(self, gt_alleles, called, sample=''):
        self.gt_alleles = gt_alleles
        self.called = called
        self.sample = sample
# Set up dummy VCF records which are just lists of genotypes
dummy_record1 = [] # Example record with real data
dummy_record1.append(DummyVCFSample(['0','1'], True, 'S1'))
dummy_record1.append(DummyVCFSample(['1','1'], True, 'S2'))
dummy_record1.append(DummyVCFSample(['1','1'], True, 'S3'))
dummy_record1.append(DummyVCFSample(['1','2'], True, 'S4'))
dummy_record1.append(DummyVCFSample(['2','2'], True, 'S5'))
dummy_record1.append(DummyVCFSample(['0'], True, 'S6')) # add a haploid sample
dummy_record2 = [] # Empty record
dummy_record3 = [] # All reference
for i in range(3): dummy_record3.append(DummyVCFSample(['0','0'], True, 'S7'))
dummy_record3.append(DummyVCFSample(['0','0','0'], True, 'S8')) # add a triploid sample
dummy_record4 = [] # Example record not called (not sure what gt field should look like)
dummy_record4.append(DummyVCFSample(['0'], False, 'S9'))

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
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec) # To test str function
    true_gt_counts = {(ref_allele, alt_alleles[0]): 1, \
                      (alt_alleles[0], alt_alleles[0]): 2, \
                      (alt_alleles[0], alt_alleles[1]): 1, \
                      (alt_alleles[1], alt_alleles[1]): 1, (ref_allele,): 1}
    true_len_gt_counts = {(3, 4): 1, (4, 4): 2, (4, 6): 1, (6, 6): 1, (3,): 1}

    gt_counts_uselength = rec.GetGenotypeCounts()
    gt_counts_nolength = rec.GetGenotypeCounts(uselength = False)
    assert (all(v == true_len_gt_counts[k] for k,v in gt_counts_uselength.items()) and len(gt_counts_uselength) == len(true_len_gt_counts))
    assert (all(v == true_gt_counts[k] for k,v in gt_counts_nolength.items()) and len(gt_counts_nolength) == len(true_gt_counts))
    
    # Test good example with samplelist
    true_gt_counts_slist = {(ref_allele, alt_alleles[0]): 1, \
                            (alt_alleles[0], alt_alleles[0]): 1, \
                            (ref_allele,): 1}
    true_len_gt_counts_slist = {(3, 4): 1, (4, 4): 1, (3,): 1}
    slist = ['S1', 'S3', 'S6']
    gt_counts_uselength_slist = rec.GetGenotypeCounts(samplelist = slist)
    gt_counts_nolength_slist = rec.GetGenotypeCounts(samplelist = slist, uselength = False)
    assert (all(v == true_len_gt_counts_slist[k] for k,v in gt_counts_uselength_slist.items()) and len(gt_counts_uselength_slist) == len(true_len_gt_counts_slist))
    assert (all(v == true_gt_counts_slist[k] for k,v in gt_counts_nolength_slist.items()) and len(gt_counts_nolength_slist) == len(true_gt_counts_slist))
    
    # Test example where alt=[None]
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_gt_counts = {(3, 3, 3): 1, (3, 3): 3}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(v == true_len_gt_counts[k] for k,v in gt_counts_uselength.items()) and len(gt_counts_uselength) == len(true_len_gt_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_gt_counts_slist = {}
    gt_counts_uselength_slist = rec.GetGenotypeCounts(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_gt_counts_slist[k] for k,v in gt_counts_uselength_slist.items()) and len(gt_counts_uselength_slist) == len(true_len_gt_counts_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [None], "CAG", "")
    true_len_gt_counts = {}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(v == true_len_gt_counts[k] for k,v in gt_counts_uselength.items()) and len(gt_counts_uselength) == len(true_len_gt_counts))
    
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
