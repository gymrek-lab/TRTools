import os, sys
import pytest
import vcf
import numpy as np

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
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec) # To test str function
    true_al_counts = {ref_allele: 2, alt_alleles[0]: 6, alt_alleles[1]: 3}
    true_len_al_counts = {3: 2, 4: 6, 6: 3}

    al_counts_uselength = rec.GetAlleleCounts()
    al_counts_nolength = rec.GetAlleleCounts(uselength = False)
    assert (all(v == true_len_al_counts[k] for k,v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))
    assert (all(v == true_al_counts[k] for k,v in al_counts_nolength.items()) and len(al_counts_nolength) == len(true_al_counts))

    
    # Test good example with samplelist
    true_al_counts_slist = {ref_allele: 2, alt_alleles[0]: 3}
    true_len_al_counts_slist = {3: 2, 4: 3}
    slist = ['S1', 'S3', 'S6']
    al_counts_uselength_slist = rec.GetAlleleCounts(samplelist = slist)
    al_counts_nolength_slist = rec.GetAlleleCounts(samplelist = slist, uselength = False)

    assert (all(v == true_len_al_counts_slist[k] for k,v in al_counts_uselength_slist.items()) and len(al_counts_uselength_slist) == len(true_len_al_counts_slist))
    assert (all(v == true_al_counts_slist[k] for k,v in al_counts_nolength_slist.items()) and len(al_counts_nolength_slist) == len(true_al_counts_slist))    

    # Test example where alt=[None]
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_al_counts = {3: 9}
    al_counts_uselength = rec.GetAlleleCounts()
    assert (all(v == true_len_al_counts[k] for k,v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_al_counts_slist = {}
    al_counts_uselength_slist = rec.GetAlleleCounts(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_counts_slist[k] for k,v in al_counts_uselength_slist.items()) and len(al_counts_uselength_slist) == len(true_len_al_counts_slist))    

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [None], "CAG", "")
    true_len_al_counts = {}
    al_counts_uselength = rec.GetAlleleCounts()
    assert (all(v == true_len_al_counts[k] for k,v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))


def test_GetAlleleFreqs():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec) # To test str function
    true_al_freqs = {ref_allele: 0.18181818181818182, alt_alleles[0]: 0.5454545454545454, alt_alleles[1]: 0.2727272727272727}
    true_len_al_freqs = {3: 0.18181818181818182, 4: 0.5454545454545454, 6: 0.2727272727272727}
    al_freqs_uselength = rec.GetAlleleFreqs()
    al_freqs_nolength = rec.GetAlleleFreqs(uselength = False)
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))
    assert (all(v == true_al_freqs[k] for k,v in al_freqs_nolength.items()) and len(al_freqs_nolength) == len(true_al_freqs))
    
    # Test good example with samplelist
    true_al_freqs_slist = {ref_allele: 0.4, alt_alleles[0]: 0.6}
    true_len_al_freqs_slist = {3: 0.4, 4: 0.6}
    slist = ['S1', 'S3', 'S6']
    al_freqs_uselength_slist = rec.GetAlleleFreqs(samplelist = slist)
    al_freqs_nolength_slist = rec.GetAlleleFreqs(samplelist = slist, uselength = False)
    assert (all(v == true_len_al_freqs_slist[k] for k,v in al_freqs_uselength_slist.items()) and len(al_freqs_uselength_slist) == len(true_len_al_freqs_slist))
    assert (all(v == true_al_freqs_slist[k] for k,v in al_freqs_nolength_slist.items()) and len(al_freqs_nolength_slist) == len(true_al_freqs_slist))    

    # Test example where alt=[None]
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_al_freqs = {3: 1.0}
    al_freqs_uselength = rec.GetAlleleFreqs()
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))


    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_len_al_freqs_slist = {}
    al_freqs_uselength_slist = rec.GetAlleleFreqs(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_freqs_slist[k] for k,v in al_freqs_uselength_slist.items()) and len(al_freqs_uselength_slist) == len(true_len_al_freqs_slist))    

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [None], "CAG", "")
    true_len_al_freqs = {}
    al_freqs_uselength = rec.GetAlleleFreqs()
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))

def test_GetMaxAllele():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec) # To test str function
    true_al_max = 6.0
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max
    
    # Test good example with samplelist
    true_al_freqs_slist = {ref_allele: 0.4, alt_alleles[0]: 0.6}
    true_len_al_freqs_slist = {3: 0.4, 4: 0.6}
    slist = ['S1', 'S3', 'S6']
    true_al_max_slist = 4.0
    al_max_slist = rec.GetMaxAllele(samplelist = slist)
    assert al_max_slist == true_al_max_slist

    # Test example where alt=[None]
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_al_max = 3.0
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [None], "CAG", "")
    true_al_max_slist = np.nan
    al_max_slist = rec.GetMaxAllele(samplelist = ['NonExistentSample'])
    assert np.isnan(al_max_slist) == True

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [None], "CAG", "")
    true_al_max = np.nan
    al_max = rec.GetMaxAllele()
    assert np.isnan(al_max) == True



#### Test TRRecordHarmonizer on different files ####
# TODO: test input a VCF that came from something else e.g. SNP calls from samtools

gangstr_path = os.path.join(VCFDIR, "test_gangstr.vcf")
hipstr_path = os.path.join(VCFDIR, "test_hipstr.vcf")
popstr_path = os.path.join(VCFDIR, "test_popstr.vcf")
advntr_path = os.path.join(VCFDIR, "test_advntr.vcf")
ehuntr_path = os.path.join(VCFDIR, "test_ExpansionHunter.vcf")
snps_path = os.path.join(VCFDIR, "snps.vcf")

gangstr_vcf = vcf.Reader(filename=gangstr_path)
hipstr_vcf = vcf.Reader(filename=hipstr_path)
popstr_vcf = vcf.Reader(filename=popstr_path)
advntr_vcf = vcf.Reader(filename=advntr_path)
ehuntr_vcf = vcf.Reader(filename=ehuntr_path)
snps_vcf = vcf.Reader(filename=snps_path)

def test_trh_init():
    # Test example with unknown VCF type given
    with pytest.raises(ValueError):
        trh.TRRecordHarmonizer(gangstr_vcf, vcftype='unknownvcf')
    # Test example with unsupported VCF
    with pytest.raises(ValueError):
        trh.TRRecordHarmonizer(snps_vcf)

    # Test examples with correct preset VCF type
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf, vcftype='gangstr')
    gangstr_infer = gangstr_trh.InferVCFType()
    assert gangstr_infer == trh.VCFTYPES.gangstr

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf, vcftype='hipstr')
    hipstr_infer = hipstr_trh.InferVCFType()
    assert hipstr_infer == trh.VCFTYPES.hipstr

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf, vcftype='popstr')
    popstr_infer = popstr_trh.InferVCFType()
    assert popstr_infer == trh.VCFTYPES.popstr

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf, vcftype='advntr')
    advntr_infer = advntr_trh.InferVCFType()
    assert advntr_infer == trh.VCFTYPES.advntr

    # No implemented yet
    #ehuntr_trh = trh.TRRecordHarmonizer(ehuntr_vcf, vcftype='eh')
    #ehuntr_infer = advntr_trh.InferVCFType()
    #assert ehuntr_infer == trh.VCFTYPES.eh


def test_InferVCFType():
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf)
    gangstr_infer = gangstr_trh.InferVCFType()
    assert gangstr_infer == trh.VCFTYPES.gangstr

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf)
    hipstr_infer = hipstr_trh.InferVCFType()
    assert hipstr_infer == trh.VCFTYPES.hipstr

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf)
    popstr_infer = popstr_trh.InferVCFType()
    assert popstr_infer == trh.VCFTYPES.popstr

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf)
    advntr_infer = advntr_trh.InferVCFType()
    assert advntr_infer == trh.VCFTYPES.advntr

    # No implemented yet
    #ehuntr_trh = trh.TRRecordHarmonizer(ehuntr_vcf)
    #ehuntr_infer = advntr_trh.InferVCFType()
    #assert ehuntr_infer == trh.VCFTYPES.eh

def test_HarmonizeRecord():
    # Gangstr
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf)

    tr_rec1 = next(iter(gangstr_trh)) 
    assert tr_rec1.ref_allele == 'tctgtctgtctg'.upper()
    assert tr_rec1.alt_alleles == [None]
    assert tr_rec1.motif == 'tctg'.upper()
    tr_rec2 = next(iter(gangstr_trh)) 
    tr_rec3 = next(iter(gangstr_trh)) 
    assert tr_rec3.ref_allele == 'tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg'.upper()
    assert tr_rec3.alt_alleles == ['tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg'.upper()]
    assert tr_rec3.motif == 'tg'.upper()



    ## hipstr
    hipstr_vcf = vcf.Reader(filename=hipstr_path)
    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf)
    
    # TODO check if these hipstr examples are correct and if they cover corner cases
    tr_rec1 = next(iter(hipstr_trh)) 
    assert tr_rec1.ref_allele == 'GGTGGTGGTGGGGGCGGTGGGGGTGGTG'
    assert tr_rec1.alt_alleles == ['GGTGGTGGTGGGGGCGGTGGTGGTGCTG']
    assert tr_rec1.motif == 'GGT'
    assert tr_rec1.record_id == 'STR_2'
    tr_rec2 = next(iter(hipstr_trh)) 
    tr_rec3 = next(iter(hipstr_trh)) 
    assert tr_rec3.ref_allele == 'TTTTTTTTTTTTTTT'
    assert tr_rec3.alt_alleles == [None]
    assert tr_rec3.motif == 'T'.upper()
    assert tr_rec3.record_id == 'STR_4'

    ##  popstr
    popstr_vcf = vcf.Reader(filename=popstr_path)
    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf)

    tr_rec1 = next(iter(popstr_trh)) 
    assert tr_rec1.ref_allele == 'GGGGGGGCGGGGGGGGGG'
    assert tr_rec1.alt_alleles == ['G' * 14, 'G' * 17]
    assert tr_rec1.motif == 'G'
    assert tr_rec1.record_id == 'chr21:5020351:M'

    tr_rec2 = next(iter(popstr_trh)) 
    tr_rec3 = next(iter(popstr_trh)) 
    assert tr_rec3.ref_allele == 'TTTTTTTTTTTTTTTTTTTTTT'
    assert tr_rec3.alt_alleles == ['T' * 21]
    assert tr_rec3.motif == 'T'
    assert tr_rec3.record_id == 'chr21:5031126:M'


    ## advntr
    advntr_vcf = vcf.Reader(filename=advntr_path)
    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf)
    
    tr_rec1 = next(iter(advntr_trh)) 
    assert tr_rec1.ref_allele == 'GCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGG'
    assert tr_rec1.alt_alleles == ['GCGCGGGGCGGGGCGCGGGGCGGG', 'GCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGG']
    assert tr_rec1.motif == 'GCGCGGGGCGGG'


    ## TODO Fix this test
    # Test examples with incorrect preset VCF type
    #ic_gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf, vcftype='advntr')
    #with pytest.raises(ValueError):
    #    tr_rec1 = next(iter(ic_gangstr_trh)) # We need to raise a value error if we don't find the correct keys and report that vcf type is probably input incorrectly
    
