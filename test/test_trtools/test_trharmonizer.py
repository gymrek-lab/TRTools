import os
import sys

import numpy as np
import pytest
import tr_harmonizer as trh
import vcf

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'trtools')
)


COMMDIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "common"
)
VCFDIR = os.path.join(COMMDIR, "sample_vcfs")

#### Test TRRecord using dummy info ####


# Set up dummy class with gt_alleles
class DummyVCFSample:
    def __init__(self, gt_alleles, called, sample=''):
        self.gt_alleles = gt_alleles
        self.called = called
        self.sample = sample


class DummyVCFRecord:
    def __init__(self):
        self.samples = []

    def __iadd__(self, other: DummyVCFSample):
        self.samples.append(other)
        return self

    def __iter__(self):
        for sample in self.samples:
            yield sample


# Set up dummy VCF records which are just lists of genotypes
dummy_record1 = DummyVCFRecord()  # Example record with real data
dummy_sample1 = DummyVCFSample(['0', '1'], True, 'S1')
dummy_sample2 = DummyVCFSample(['1', '1'], True, 'S2')
dummy_record1 += dummy_sample1
dummy_record1 += dummy_sample2
dummy_record1 += DummyVCFSample(['1', '1'], True, 'S3')
dummy_record1 += DummyVCFSample(['1', '2'], True, 'S4')
dummy_record1 += DummyVCFSample(['2', '2'], True, 'S5')
dummy_record1 += DummyVCFSample(['0'], True, 'S6')  # add a haploid sample
dummy_record2 = DummyVCFRecord()  # Empty record
dummy_record3 = DummyVCFRecord()  # All reference
for i in range(3):
    dummy_record3 += DummyVCFSample(['0', '0'], True, 'S7')
# add a triploid sample
dummy_record3 += DummyVCFSample(['0', '0', '0'], True, 'S8')
# Example record not called (not sure what gt field should look like)
dummy_record4 = DummyVCFRecord()
dummy_record4 += DummyVCFSample(['0'], False, 'S9')


def test_TRRecord_iter():
    record = trh.TRRecord(dummy_record1, "ACG", ["A", "C", "G", "T"],
                          "FOO", "BAR")
    record_iter = iter(record)
    assert next(record_iter) == dummy_sample1
    assert next(record_iter) == dummy_sample2


def test_TRRecord_allele_lengths():
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    motif = 'FOO'
    ID = 'BAR'

    # alt alleles
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
                     alt_allele_lengths=[4, 6])
    record = trh.TRRecord(dummy_record1, ref_allele, None, motif, ID,
                          alt_allele_lengths=[4, 5.5])
    assert record.alt_alleles == [motif * 4, motif * 5 + "F"]

    # ref allele
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
                     ref_allele_length=5)
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, None, alt_alleles, motif, ID,
                     ref_allele_length=5)
    record = trh.TRRecord(dummy_record1, None, None, motif, ID,
                          ref_allele_length=5.5, alt_allele_lengths=[4, 5.5])
    assert record.ref_allele == motif*5 + 'F'


def test_TRRecord_unique_lengths():
    record = trh.TRRecord(
        dummy_record2,
        "ACGACGACG",
        [
            "ACGAAGACG",
            "ACGACGACGACG",
            "ACGACGACAACG"
        ],
        "ACG",
        "ACG-repeat"
    )

    assert record.UniqueLengthGenotypes() == {0, 2}
    assert record.UniqueLengthGenotypeMapping() == {
        0: 0,
        1: 0,
        2: 2,
        3: 2
    }


def test_TRRecord_full_alleles():
    full_ref = "TCAGCAGCAGA"
    full_alts = [
        "ACAGCAGCAGCAGC",
        "ACAGCAGCAGCAGCAGCAGG",
        "ACAGCAGCAGCAGG",
        "ACAGCAGCAGG",
        "TCAGCAGG",
    ]
    ref_allele = full_ref[1:-1]
    alt_alleles = []
    for allele in full_alts:
        alt_alleles.append(allele[1:-1])
    motif = 'FOO'
    ID = 'BAR'

    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, None, None, motif, ID,
                     full_alleles=(full_ref, full_alts))
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
                     full_alleles=(["CAGCAGCAQQQQQQQQQQQQQQQ"], full_alts))
    with pytest.raises(ValueError):
        bad_alts = [
            "CAGCAGCAQQQQQQQQQQQQQQQ",
            full_alts[1]
        ]
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
                     full_alleles=(ref_allele, bad_alts))

    record = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
                          full_alleles=(ref_allele, alt_alleles))

    assert record.UniqueStringGenotypes() == {0, 1, 2, 5}
    assert record.UniqueStringGenotypeMapping() == {
        0: 0,
        1: 1,
        2: 2,
        3: 1,
        4: 0,
        5: 5
    }


def test_TRRecord_GetGenotypes():
    # Test good example
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec)  # To test str function
    true_gts = [[ref_allele, alt_alleles[0]],
                [alt_alleles[0], alt_alleles[0]],
                [alt_alleles[0], alt_alleles[0]],
                [alt_alleles[0], alt_alleles[1]],
                [alt_alleles[1], alt_alleles[1]], [ref_allele]]
    true_len_gts = [[3, 4], [4, 4], [4, 4], [4, 6], [6, 6], [3]]
    ind = 0
    for sample in rec.vcfrecord:
        stringt = rec.GetStringGenotype(sample)
        lengt = rec.GetLengthGenotype(sample)
        assert(all(
            [(stringt[i] == true_gts[ind][i]) for i in range(len(stringt))]
        ))
        assert(all(
            [(lengt[i] == true_len_gts[ind][i]) for i in range(len(lengt))]
        ))
        ind += 1
    # Test example where alt=[]
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    print(rec)  # To test str function
    for sample in rec.vcfrecord:
        stringt = rec.GetStringGenotype(sample)
        lengt = rec.GetLengthGenotype(sample)
        assert(all([item == ref_allele for item in stringt]))
        assert(all([item == 3 for item in lengt]))
    # Test example with discrepancy between alt_alleles and genotypes given
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, [], "CAG", "")


def test_GetGenotypeCounts():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec)  # To test str function
    true_gt_counts = {(ref_allele, alt_alleles[0]): 1,
                      (alt_alleles[0], alt_alleles[0]): 2,
                      (alt_alleles[0], alt_alleles[1]): 1,
                      (alt_alleles[1], alt_alleles[1]): 1, (ref_allele,): 1}
    true_len_gt_counts = {(3, 4): 1, (4, 4): 2, (4, 6): 1, (6, 6): 1, (3,): 1}

    gt_counts_uselength = rec.GetGenotypeCounts()
    gt_counts_nolength = rec.GetGenotypeCounts(uselength=False)
    assert (all(
        v == true_len_gt_counts[k] for k, v in gt_counts_uselength.items()
    ) and len(gt_counts_uselength) == len(true_len_gt_counts))
    assert (all(
        v == true_gt_counts[k] for k, v in gt_counts_nolength.items()
    ) and len(gt_counts_nolength) == len(true_gt_counts))

    # Test good example with samplelist
    true_gt_counts_slist = {(ref_allele, alt_alleles[0]): 1,
                            (alt_alleles[0], alt_alleles[0]): 1,
                            (ref_allele,): 1}
    true_len_gt_counts_slist = {(3, 4): 1, (4, 4): 1, (3,): 1}
    slist = ['S1', 'S3', 'S6']
    gt_counts_uselength_slist = rec.GetGenotypeCounts(samplelist=slist)
    gt_counts_nolength_slist = rec.GetGenotypeCounts(samplelist=slist,
                                                     uselength=False)
    assert (all(
        v == true_len_gt_counts_slist[k] for k, v in
        gt_counts_uselength_slist.items()
    ) and len(gt_counts_uselength_slist) == len(true_len_gt_counts_slist))
    assert (all(
        v == true_gt_counts_slist[k] for k, v
        in gt_counts_nolength_slist.items()
    ) and len(gt_counts_nolength_slist) == len(true_gt_counts_slist))

    # Test example where alt=[]
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_gt_counts = {(3, 3, 3): 1, (3, 3): 3}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(
        v == true_len_gt_counts[k] for k, v in gt_counts_uselength.items()
    ) and len(gt_counts_uselength) == len(true_len_gt_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_gt_counts_slist = {}
    gt_counts_uselength_slist = \
        rec.GetGenotypeCounts(samplelist=['NonExistentSample'])
    assert (all(
        v == true_len_gt_counts_slist[k] for k, v
        in gt_counts_uselength_slist.items()
    ) and len(gt_counts_uselength_slist) == len(true_len_gt_counts_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "")
    true_len_gt_counts = {}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(
        v == true_len_gt_counts[k] for k, v in gt_counts_uselength.items()
    ) and len(gt_counts_uselength) == len(true_len_gt_counts))


def test_GetAlleleCounts():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "")
    print(rec)  # To test str function
    true_al_counts = {ref_allele: 2, alt_alleles[0]: 6, alt_alleles[1]: 3}
    true_len_al_counts = {3: 2, 4: 6, 6: 3}

    al_counts_uselength = rec.GetAlleleCounts()
    al_counts_nolength = rec.GetAlleleCounts(uselength=False)
    assert (all(
        v == true_len_al_counts[k] for k, v in al_counts_uselength.items()
    ) and len(al_counts_uselength) == len(true_len_al_counts))
    assert (all(
        v == true_al_counts[k] for k, v in al_counts_nolength.items()
    ) and len(al_counts_nolength) == len(true_al_counts))

    # Test good example with samplelist
    true_al_counts_slist = {ref_allele: 2, alt_alleles[0]: 3}
    true_len_al_counts_slist = {3: 2, 4: 3}
    slist = ['S1', 'S3', 'S6']
    al_counts_uselength_slist = rec.GetAlleleCounts(samplelist=slist)
    al_counts_nolength_slist = rec.GetAlleleCounts(samplelist=slist,
                                                   uselength=False)

    assert (all(
        v == true_len_al_counts_slist[k] for k, v
        in al_counts_uselength_slist.items()
    ) and len(al_counts_uselength_slist) == len(true_len_al_counts_slist))
    assert (all(
        v == true_al_counts_slist[k] for k, v
        in al_counts_nolength_slist.items()
    ) and len(al_counts_nolength_slist) == len(true_al_counts_slist))

    # Test example where alt=[]
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_al_counts = {3: 9}
    al_counts_uselength = rec.GetAlleleCounts()
    assert (all(v == true_len_al_counts[k] for k,
                v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_al_counts_slist = {}
    al_counts_uselength_slist = rec.GetAlleleCounts(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_counts_slist[k] for k,v in al_counts_uselength_slist.items()) and len(al_counts_uselength_slist) == len(true_len_al_counts_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "")
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

    # Test example where alt=[]
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_al_freqs = {3: 1.0}
    al_freqs_uselength = rec.GetAlleleFreqs()
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))


    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_len_al_freqs_slist = {}
    al_freqs_uselength_slist = rec.GetAlleleFreqs(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_freqs_slist[k] for k,v in al_freqs_uselength_slist.items()) and len(al_freqs_uselength_slist) == len(true_len_al_freqs_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "")
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

    # Test example where alt=[]
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_al_max = 3.0
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "")
    true_al_max_slist = np.nan
    al_max_slist = rec.GetMaxAllele(samplelist = ['NonExistentSample'])
    assert np.isnan(al_max_slist) == True

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "")
    true_al_max = np.nan
    al_max = rec.GetMaxAllele()
    assert np.isnan(al_max) == True



#### Test TRRecordHarmonizer on different files ####
# TODO: test input a VCF that came from something else e.g. SNP calls from samtools

gangstr_path = os.path.join(VCFDIR, "test_gangstr.vcf")
hipstr_path = os.path.join(VCFDIR, "test_hipstr.vcf")
popstr_path = os.path.join(VCFDIR, "test_popstr.vcf")
advntr_path = os.path.join(VCFDIR, "test_advntr.vcf")
eh_path = os.path.join(VCFDIR, "test_ExpansionHunter.vcf")
snps_path = os.path.join(VCFDIR, "snps.vcf")


def reset_vcfs():
    global gangstr_vcf, hipstr_vcf, popstr_vcf, advntr_vcf, eh_vcf, snps_vcf
    gangstr_vcf = vcf.Reader(filename=gangstr_path)
    hipstr_vcf = vcf.Reader(filename=hipstr_path)
    popstr_vcf = vcf.Reader(filename=popstr_path)
    advntr_vcf = vcf.Reader(filename=advntr_path)
    eh_vcf = vcf.Reader(filename=eh_path)
    snps_vcf = vcf.Reader(filename=snps_path)


def test_trh_init_and_type_infer():
    reset_vcfs()

    # Test example with unknown VCF type given
    with pytest.raises(ValueError):
        trh.TRRecordHarmonizer(gangstr_vcf, vcftype='unknownvcf')
    # Test example with unsupported VCF
    with pytest.raises(ValueError):
        trh.TRRecordHarmonizer(snps_vcf)

    with pytest.raises(ValueError):
        trh.MayHaveImpureRepeats('foo')
    with pytest.raises(ValueError):
        trh.HasLengthRefGenotype('foo')
    with pytest.raises(ValueError):
        trh.HasLengthAltGenotypes('foo')

    # Test examples with correct preset VCF type
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf, vcftype='gangstr')
    assert gangstr_trh.vcftype == trh.VCFTYPES.gangstr
    assert trh.InferVCFType(gangstr_vcf) == trh.VCFTYPES.gangstr
    assert (not gangstr_trh.MayHaveImpureRepeats()
            and not trh.MayHaveImpureRepeats(trh.VCFTYPES.gangstr))
    assert (not gangstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VCFTYPES.gangstr))
    assert (not gangstr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VCFTYPES.gangstr))

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf, vcftype='hipstr')
    assert hipstr_trh.vcftype == trh.VCFTYPES.hipstr
    assert trh.InferVCFType(hipstr_vcf) == trh.VCFTYPES.hipstr
    assert (hipstr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VCFTYPES.hipstr))
    assert (not hipstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VCFTYPES.hipstr))
    assert (not hipstr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VCFTYPES.hipstr))

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf, vcftype='popstr')
    assert popstr_trh.vcftype == trh.VCFTYPES.popstr
    assert trh.InferVCFType(popstr_vcf) == trh.VCFTYPES.popstr
    assert (popstr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VCFTYPES.popstr))
    assert (not popstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VCFTYPES.popstr))
    assert (popstr_trh.HasLengthAltGenotypes()
            and trh.HasLengthAltGenotypes(trh.VCFTYPES.popstr))

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf, vcftype='advntr')
    assert advntr_trh.vcftype == trh.VCFTYPES.advntr
    assert trh.InferVCFType(advntr_vcf) == trh.VCFTYPES.advntr
    assert (advntr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VCFTYPES.advntr))
    assert (not advntr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VCFTYPES.advntr))
    assert (not advntr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VCFTYPES.advntr))

    eh_trh = trh.TRRecordHarmonizer(eh_vcf, vcftype='eh')
    assert eh_trh.vcftype == trh.VCFTYPES.eh
    assert trh.InferVCFType(eh_vcf) == trh.VCFTYPES.eh
    assert (not eh_trh.MayHaveImpureRepeats()
            and not trh.MayHaveImpureRepeats(trh.VCFTYPES.eh))
    assert (eh_trh.HasLengthRefGenotype()
            and trh.HasLengthRefGenotype(trh.VCFTYPES.eh))
    assert (eh_trh.HasLengthAltGenotypes()
            and trh.HasLengthAltGenotypes(trh.VCFTYPES.eh))


def test_HarmonizeRecord():
    reset_vcfs()

    # Unknown type
    with pytest.raises(ValueError):
        trh.HarmonizeRecord("foo", next(snps_vcf))

    # Gangstr
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf)

    tr_rec1 = next(iter(gangstr_trh))
    assert tr_rec1.ref_allele == 'tctgtctgtctg'.upper()
    assert tr_rec1.alt_alleles == []
    assert tr_rec1.motif == 'tctg'.upper()
    assert not tr_rec1.HasFullStringGenotypes()
    assert not tr_rec1.HasFabricatedRefAllele()
    assert not tr_rec1.HasFabricatedAltAlleles()
    tr_rec2 = next(iter(gangstr_trh))
    tr_rec3 = next(iter(gangstr_trh))
    assert (tr_rec3.ref_allele ==
            'tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg'.upper())
    assert (tr_rec3.alt_alleles ==
            ['tgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtgtg'.upper()])
    assert tr_rec3.motif == 'tg'.upper()

    # hipstr
    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf)

    str_iter = iter(hipstr_trh)
    tr_rec1 = next(str_iter)
    assert tr_rec1.ref_allele == 'GGTGGTGGTGGGGGCGGTGGGGGTGGTG'
    assert tr_rec1.alt_alleles == ['GGTGGTGGTGGGGGCGGTGGTGGTGCTG']
    assert tr_rec1.motif == 'GGT'
    assert tr_rec1.record_id == 'STR_2'
    assert not tr_rec1.HasFullStringGenotypes()
    assert not tr_rec1.HasFabricatedRefAllele()
    assert not tr_rec1.HasFabricatedAltAlleles()
    tr_rec2 = next(str_iter)
    tr_rec3 = next(str_iter)
    assert tr_rec3.ref_allele == 'TTTTTTTTTTTTTTT'
    assert tr_rec3.alt_alleles == []
    assert tr_rec3.motif == 'T'.upper()
    assert tr_rec3.record_id == 'STR_4'
    record = next(str_iter)
    while record.record_id != "STR_125":
        record = next(str_iter)
    assert record.HasFullStringGenotypes()
    # TODO this isn't really the correct behavior -
    # we're trimming off an extra repeat from the alt allele
    assert record.full_alleles == (
        "TGCATATATGTATAATATATATTATATATGGA",
        ["TCCATATATGCATAATATATATTATATATATG"],
    )
    assert (record.ref_allele ==
            "ATATATGTATAATATATATTATATAT")
    assert (record.alt_alleles ==
            ["ATATATGCATAATATATATTATATAT"])

    # popstr
    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf)

    tr_rec1 = next(iter(popstr_trh))
    assert tr_rec1.ref_allele == 'GGGGGGGCGGGGGGGGGG'
    assert tr_rec1.alt_alleles == ['G' * 14, 'G' * 17]
    assert tr_rec1.motif == 'G'
    assert tr_rec1.record_id == 'chr21:5020351:M'
    assert not tr_rec1.HasFullStringGenotypes()
    assert not tr_rec1.HasFabricatedRefAllele()
    assert tr_rec1.HasFabricatedAltAlleles()
    tr_rec2 = next(iter(popstr_trh))
    tr_rec3 = next(iter(popstr_trh))
    assert tr_rec3.ref_allele == 'TTTTTTTTTTTTTTTTTTTTTT'
    assert tr_rec3.alt_alleles == ['T' * 21]
    assert tr_rec3.motif == 'T'
    assert tr_rec3.record_id == 'chr21:5031126:M'

    # advntr
    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf)

    tr_rec1 = next(iter(advntr_trh))
    assert tr_rec1.ref_allele == 'GCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGG'
    assert tr_rec1.alt_alleles == ['GCGCGGGGCGGGGCGCGGGGCGGG', 'GCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGGGCGCGGGGCGGG']
    assert tr_rec1.motif == 'GCGCGGGGCGGG'
    assert not tr_rec1.HasFullStringGenotypes()
    assert not tr_rec1.HasFabricatedRefAllele()
    assert not tr_rec1.HasFabricatedAltAlleles()

    # advntr
    eh_trh = trh.TRRecordHarmonizer(eh_vcf)

    tr_rec1 = next(iter(eh_trh))
    motif = 'CAG'
    assert tr_rec1.ref_allele == motif*19
    assert tr_rec1.alt_alleles == [motif*16, motif*18]
    assert tr_rec1.motif == motif
    assert tr_rec1.record_id == 'HTT'
    assert not tr_rec1.HasFullStringGenotypes()
    assert tr_rec1.HasFabricatedRefAllele()
    assert tr_rec1.HasFabricatedAltAlleles()

    ## TODO Fix this test
    # Test examples with incorrect preset VCF type
    #ic_gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf, vcftype='advntr')
    #with pytest.raises(ValueError):
    #    tr_rec1 = next(iter(ic_gangstr_trh)) # We need to raise a value error if we don't find the correct keys and report that vcf type is probably input incorrectly

