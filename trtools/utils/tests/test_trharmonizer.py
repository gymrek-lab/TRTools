import os
import types
from typing import List

import cyvcf2
import numpy as np
import pytest
from pytest import approx

import trtools.utils.tr_harmonizer as trh

# pylint: disable=C0116,C0103

#### Test TRRecord using dummy info ####

# TODO add tests that check this for quality scores?

class DummyCyvcf2Record:
    def __init__(self,
                 gts, # arraylike (n x 2). All entries in a row should contain -1s for no calls
                      # For calls with different ploidy, the array should
                      # contain columns equal to the greatest ploidy and lower
                      # ploidy calls should be padded with -1s
                      # if this is a vcf with no samples in it (just loci
                      # descriptions) then this should be None
                 ref,
                 alt
                 ):
        self.POS = 42
        self.CHROM = '1984'
        self.FORMAT = {}
        self.INFO = {}
        self.ALT = list(alt)
        self.REF = ref

        if gts is not None:
            self.genotype = types.SimpleNamespace()
            self._gts = np.array(gts)
            self._gts = np.concatenate(
                (self._gts, np.zeros((self._gts.shape[0], 1))),
                axis=1
            ) # add the phasing axis, we're not testing that here
            self.genotype.array = lambda: self._gts
        else:
            self.genotype = None


# Set up dummy VCF records which are just lists of genotypes
dummy_record_gts = [
    [0, 1],
    [1, 1],
    [1, 1],
    [1, 2],
    [2, 2],
    [0, -1]]
def get_dummy_record():
    return DummyCyvcf2Record(
        gts=dummy_record_gts, #last sample is haploid
        ref="CAGCAGCAG",
        alt=("CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG")
    )

triploid_gts = np.array([
    [0, 0, -2],
    [0, 0, -2],
    [0, 0, -2],
    [0, 0, 0]
])
def get_triploid_record():
    return DummyCyvcf2Record(
        gts=triploid_gts,
        ref="CAGCAGCAG",
        alt=[]
    )

def get_nocall_record():
    return DummyCyvcf2Record(
        gts=[[-1, -1]],
        ref="CAGCAGCAG",
        alt=[]
    )

def get_record_with_nosamples():
    return DummyCyvcf2Record(
        gts=None,
        ref="CAGCAGCAG",
        alt=["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    )


def test_unexpected_vcf_type():
    with pytest.raises(ValueError):
        trh._UnexpectedTypeError(trh.VcfTypes.gangstr)  # pylint: disable=W0212


def test_TRRecord_print():
    ref = "ABC"
    alt = ["DEF", "GHI"]
    motif = "foo"
    ID = "bar"
    test_rec = DummyCyvcf2Record(None, ref, alt)

    # default format
    record = trh.TRRecord(test_rec, ref, alt, motif, ID, None)
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    # having a quality field doesn't change the format
    record = trh.TRRecord(test_rec, ref, alt, motif, ID, "some_field")
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    # use CHROM:POS if records do not have Ids
    record = trh.TRRecord(test_rec, ref, alt, motif, None, None)
    assert str(record) == "{}:{} {} {} {},{}".format(test_rec.CHROM,
                                                     test_rec.POS,
                                                     motif, ref, alt[0],
                                                     alt[1])


    # full alleles are printed over short alleles when present
    record = trh.TRRecord(test_rec, "B", ["E", "H"], motif, ID, None,
                          full_alleles=(ref, alt))
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    # alt allele lens are printed when seqs are not present
    record = trh.TRRecord(test_rec, ref, None, motif, ID, None,
                          alt_allele_lengths=[3, 5.5])
    assert str(record) == "{} {} {} n_reps:3,n_reps:5.5".format(ID, motif, ref)

    # ref allele len is printed when seq is not present
    record = trh.TRRecord(test_rec, None, None, motif, ID, None,
                          ref_allele_length=7,
                          alt_allele_lengths=[3, 5.5])
    assert str(record) == ("{} {} n_reps:7 n_reps:3,n_reps:5.5"
                           .format(ID, motif))


def test_TRRecord_output_shape():
    rec = get_dummy_record()
    record = trh.TRRecord(rec, rec.REF, rec.ALT,
                          "CAG", "STR1", None)
    assert (record.GetGenotypeIndicies().shape ==
            rec.genotype.array().shape)


def test_TRRecord_allele_seqs_from_lens():
    rec = get_record_with_nosamples()
    ref_allele = rec.REF
    alt_alleles = rec.ALT
    motif = 'CAG'
    ID = 'BAR'

    # alt alleles
    # confirm cannot use alt_allele_lengths and alt_allele sequences at the same time
    with pytest.raises(ValueError):
        trh.TRRecord(rec, ref_allele, alt_alleles, motif, ID,
                     "some_field",
                     alt_allele_lengths=[4, 6])

    # confirm accessing fabricated alt_alleles from lengths
    record = trh.TRRecord(rec, ref_allele, None, motif, ID,
                          "some_field",
                          alt_allele_lengths=[4, 5.5])
    assert record.alt_alleles == [motif * 4, motif * 5 + "C"]

    # ref allele
    # confirm cannot use ref_allele_length and ref_allele sequence at the same time
    with pytest.raises(ValueError):
        trh.TRRecord(rec, ref_allele, alt_alleles, motif, ID, None,
                     ref_allele_length=5)

    # confirm cannot use ref_allele_length and alt_allele sequences at the same time
    with pytest.raises(ValueError):
        trh.TRRecord(rec, None, alt_alleles, motif, ID, None,
                     ref_allele_length=5)

    # confirm accessing fabricated ref_allele from length
    record = trh.TRRecord(rec, None, None, motif, ID, None,
                          ref_allele_length=5.5, alt_allele_lengths=[4, 5.5])
    assert record.ref_allele == motif * 5 + 'C'


def test_TRRecord_unique_lengths():
    alts = [
        "ACGAAGACG",
        "ACGACGACGACG",
        "ACGACGACAACG"
    ]
    ref = "ACGACGACG"
    test_rec = DummyCyvcf2Record(None, ref, alts)
    record = trh.TRRecord(
        test_rec,
        ref,
        alts,
        "ACG",
        "ACG-repeat",
        None
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
    motif = 'CAG'
    ID = 'STR1'
    test_rec = DummyCyvcf2Record(None, full_ref, full_alts)

    # confirm that you need to set ref_allele if you wish to set full alleles
    with pytest.raises(ValueError):
        trh.TRRecord(test_rec, None, alt_alleles, motif, ID, None,
                     full_alleles=(full_ref, full_alts))
    # confirm that you need to set alt_alleles if you wish to set full alleles
    with pytest.raises(ValueError):
        trh.TRRecord(test_rec, ref_allele, None, motif, ID, None,
                     full_alleles=(full_ref, full_alts))
    # confirm that each allele in full_alleles must contain its corresponding
    # non-full allele as a substring
    with pytest.raises(ValueError):
        trh.TRRecord(test_rec, ref_allele, alt_alleles, motif, ID, None,
                     full_alleles=(["CAGCAGCAQQQQQQQQQQQQQQQ"], full_alts))
    with pytest.raises(ValueError):
        bad_alts = [
            "CAGCAGCAQQQQQQQQQQQQQQQ",
            full_alts[1]
        ]
        trh.TRRecord(test_rec, ref_allele, alt_alleles, motif, ID, None,
                     full_alleles=(ref_allele, bad_alts))

    record = trh.TRRecord(test_rec, ref_allele, alt_alleles, motif, ID,
                          None,
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
    dummy_record = get_dummy_record()
    # Test good example
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    true_gts = [[ref_allele, alt_alleles[0]],
                [alt_alleles[0], alt_alleles[0]],
                [alt_alleles[0], alt_alleles[0]],
                [alt_alleles[0], alt_alleles[1]],
                [alt_alleles[1], alt_alleles[1]],
                [ref_allele, '.']]
    true_gts = np.array(true_gts)
    true_len_gts = np.array([[3, 4], [4, 4], [4, 4], [4, 6], [6, 6], [3, -1]])
    assert np.all(rec.GetGenotypeIndicies()[:, :-1] ==
                  np.array(dummy_record_gts))
    assert np.all(rec.GetLengthGenotypes()[:, :-1] == true_len_gts)
    assert np.all(rec.GetStringGenotypes()[:, :-1] == true_gts)

    # Test example where alt=[]
    triploid_record = get_triploid_record()
    rec = trh.TRRecord(triploid_record, ref_allele, [], "CAG", "", None)
    true_len_gts = [[3, 3, -2],
                [3, 3, -2],
                [3, 3, -2],
                [3, 3, 3]]
    true_len_gts = np.array(true_len_gts)
    true_idx_gts = [[0, 0, -2],
                    [0, 0, -2],
                    [0, 0, -2],
                    [0, 0, 0]]
    true_idx_gts = np.array(true_idx_gts)
    true_gts = [[ref_allele, ref_allele, ','],
                [ref_allele, ref_allele, ','],
                [ref_allele, ref_allele, ','],
                [ref_allele, ref_allele, ref_allele]]
    true_gts = np.array(true_gts)

    assert np.all(rec.GetGenotypeIndicies()[:, :-1] == true_idx_gts)
    assert np.all(rec.GetLengthGenotypes()[:, :-1] == true_len_gts)
    assert np.all(rec.GetStringGenotypes()[:, :-1] == true_gts)

    # Test example with fewer alt_alleles than the max genotype index
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record, ref_allele, [], "CAG", "", None)


def _dicts_equal(dict1, dict2):
    return (all(k in dict2 and v == dict2[k] for k,v in dict1.items()) and
            len(dict1) == len(dict2))

def test_GetGenotypeCounts():
    # Test working example, no sample_index,
    # also a call with a missing haplotype
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    true_idx_gt_counts = {(0, 1) : 1,
                          (1, 1) : 2,
                          (1, 2) : 1,
                          (2, 2) : 1}
    true_gt_counts = {(ref_allele, alt_alleles[0]): 1,
                      (alt_alleles[0], alt_alleles[0]): 2,
                      (alt_alleles[0], alt_alleles[1]): 1,
                      (alt_alleles[1], alt_alleles[1]): 1}
    true_len_gt_counts = {
        (3, 4): 1,
        (4, 4): 2,
        (4, 6): 1,
        (6, 6): 1
    }

    gt_counts_uselength = rec.GetGenotypeCounts()
    gt_counts_nolength = rec.GetGenotypeCounts(uselength=False)
    gt_counts_idx = rec.GetGenotypeCounts(index=True)
    assert _dicts_equal(true_idx_gt_counts, gt_counts_idx)
    assert _dicts_equal(true_len_gt_counts, gt_counts_uselength)
    assert _dicts_equal(true_gt_counts, gt_counts_nolength)

    # Test example where alt=[]
    triploid_record = get_triploid_record()
    rec = trh.TRRecord(triploid_record, ref_allele, [], "CAG", "", None)
    true_idx_gt_counts = {(0, 0, 0): 1, (-2, 0, 0): 3}
    true_len_gt_counts = {(3, 3, 3): 1, (-2, 3, 3): 3}
    true_gt_counts = {(ref_allele, ref_allele, ref_allele): 1,
                      (',', ref_allele, ref_allele): 3}
    gt_counts_idx = rec.GetGenotypeCounts(index=True)
    gt_counts_uselength = rec.GetGenotypeCounts()
    gt_counts_nolength = rec.GetGenotypeCounts(uselength=False)
    assert _dicts_equal(true_idx_gt_counts, gt_counts_idx)
    assert _dicts_equal(true_len_gt_counts, gt_counts_uselength)
    assert _dicts_equal(true_gt_counts, gt_counts_nolength)

    # Test example where there are no samples
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    # expect a entry for counting no-calls
    assert len(rec.GetGenotypeCounts(index=True)) == 0
    assert len(rec.GetGenotypeCounts()) == 0
    assert len(rec.GetGenotypeCounts(uselength=True)) == 0
    assert len(rec.GetGenotypeCounts(index=True, include_nocalls=True)) == 1
    assert len(rec.GetGenotypeCounts(include_nocalls=True)) == 1
    assert len(rec.GetGenotypeCounts(uselength=True, include_nocalls=True)) == 1

    # Test working example with sample_index
    true_idx_gt_counts_sindex = {(0, 1): 1, (1, 1): 1}
    true_gt_counts_sindex = {(ref_allele, alt_alleles[0]): 1,
                            (alt_alleles[0], alt_alleles[0]): 1}
    true_len_gt_counts_sindex = {(3, 4): 1, (4, 4): 1}
    sindex = [0, 2, 5]
    dummy_record = get_dummy_record()
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    gt_counts_idx_sindex = rec.GetGenotypeCounts(sample_index=sindex, index=True)
    gt_counts_uselength_sindex = rec.GetGenotypeCounts(sample_index=sindex)
    gt_counts_nolength_sindex = rec.GetGenotypeCounts(sample_index=sindex,
                                                     uselength=False)
    assert _dicts_equal(true_idx_gt_counts_sindex, gt_counts_idx_sindex)
    assert _dicts_equal(true_len_gt_counts_sindex, gt_counts_uselength_sindex)
    assert _dicts_equal(true_gt_counts_sindex, gt_counts_nolength_sindex)


def test_GetAlleleCounts():
    # Test working example, no sample_index
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    true_al_counts = {ref_allele: 2, alt_alleles[0]: 6, alt_alleles[1]: 3}
    true_len_al_counts = {3: 2, 4: 6, 6: 3}
    true_idx_al_counts = {0: 2, 1: 6, 2: 3}

    al_counts_uselength = rec.GetAlleleCounts()
    al_counts_nolength = rec.GetAlleleCounts(uselength=False)
    al_counts_idx = rec.GetAlleleCounts(index=True)
    assert _dicts_equal(true_idx_al_counts, al_counts_idx)
    assert _dicts_equal(true_len_al_counts, al_counts_uselength)
    assert _dicts_equal(true_al_counts, al_counts_nolength)

    # Test example where alt=[]
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    true_len_al_counts = {3: 9}
    true_idx_al_counts = {0: 9}
    true_al_counts = {ref_allele: 9}
    al_counts_idx = rec.GetAlleleCounts(index=True)
    al_counts_uselength = rec.GetAlleleCounts()
    al_counts_nolength = rec.GetAlleleCounts(uselength=False)
    assert _dicts_equal(true_idx_al_counts, al_counts_idx)
    assert _dicts_equal(true_len_al_counts, al_counts_uselength)
    assert _dicts_equal(true_al_counts, al_counts_nolength)

    # Test example where there are no samples
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    assert len(rec.GetAlleleCounts(index=True)) == 0
    assert len(rec.GetAlleleCounts()) == 0
    assert len(rec.GetAlleleCounts(uselength=True)) == 0

    # Test working example with sample_index
    true_idx_al_counts_sindex = {0: 2, 1: 3}
    true_len_al_counts_sindex = {3: 2, 4: 3}
    true_al_counts_sindex = {ref_allele: 2, alt_alleles[0]: 3}
    sindex = [0, 2, 5]
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    al_counts_idx_sindex = rec.GetAlleleCounts(sample_index=sindex, index=True)
    al_counts_uselength_sindex = rec.GetAlleleCounts(sample_index=sindex)
    al_counts_nolength_sindex = rec.GetAlleleCounts(sample_index=sindex,
                                                   uselength=False)
    assert _dicts_equal(true_idx_al_counts_sindex, al_counts_idx_sindex)
    assert _dicts_equal(true_len_al_counts_sindex, al_counts_uselength_sindex)
    assert _dicts_equal(true_al_counts_sindex, al_counts_nolength_sindex)


def test_GetAlleleFreqs():
    # Test working example, no sample_index
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    true_al_freqs = {
        ref_allele: 0.1818181,
        alt_alleles[0]: 0.5454545,
        alt_alleles[1]: 0.2727272
    }
    true_len_al_freqs = {3: 0.1818181, 4: 0.5454545, 6: 0.2727272}
    true_idx_al_freqs = {0: 0.1818181, 1: 0.5454545, 2: 0.2727272}

    al_freqs_uselength = rec.GetAlleleFreqs()
    al_freqs_nolength = rec.GetAlleleFreqs(uselength=False)
    al_freqs_idx = rec.GetAlleleFreqs(index=True)
    assert (all(
        v == approx(true_idx_al_freqs[k]) for k, v in al_freqs_idx.items()
    ) and len(al_freqs_idx) == len(true_idx_al_freqs))
    assert (all(
        v == approx(true_len_al_freqs[k]) for k, v in al_freqs_uselength.items()
    ) and len(al_freqs_uselength) == len(true_len_al_freqs))
    assert (all(
        v == approx(true_al_freqs[k]) for k, v in al_freqs_nolength.items()
    ) and len(al_freqs_nolength) == len(true_al_freqs))

    # Test example where alt=[]
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    true_len_al_freq = {3: 1}
    true_idx_al_freq = {0: 1}
    true_al_freq = {ref_allele: 1}
    al_freq_idx = rec.GetAlleleFreqs(index=True)
    al_freq_uselength = rec.GetAlleleFreqs()
    al_freq_nolength = rec.GetAlleleFreqs(uselength=False)
    assert (all(
        v == true_idx_al_freq[k] for k, v in al_freq_idx.items()
    ) and len(al_freq_idx) == len(true_idx_al_freq))
    assert (all(
        v == true_len_al_freq[k] for k, v in al_freq_uselength.items()
    ) and len(al_freq_uselength) == len(true_len_al_freq))
    assert (all(
        v == true_al_freq[k] for k, v in al_freq_nolength.items()
    ) and len(al_freq_nolength) == len(true_al_freq))

    # Test example where there are no samples
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    assert len(rec.GetAlleleFreqs(index=True)) == 0
    assert len(rec.GetAlleleFreqs()) == 0
    assert len(rec.GetAlleleFreqs(uselength=True)) == 0

    # Test working example with sample_index
    true_idx_al_freqs_sindex = {0: 0.4, 1: 0.6}
    true_len_al_freqs_sindex = {3: 0.4, 4: 0.6}
    true_al_freqs_sindex = {ref_allele: 0.4, alt_alleles[0]: 0.6}
    sindex = [0, 2, 5]
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    al_freqs_idx_sindex = rec.GetAlleleFreqs(sample_index=sindex, index=True)
    al_freqs_uselength_sindex = rec.GetAlleleFreqs(sample_index=sindex)
    al_freqs_nolength_sindex = rec.GetAlleleFreqs(sample_index=sindex,
                                                 uselength=False)
    assert (all(
        v == true_idx_al_freqs_sindex[k] for k, v in
        al_freqs_idx_sindex.items()
    ) and len(al_freqs_idx_sindex) == len(true_idx_al_freqs_sindex))
    assert (all(
        v == true_len_al_freqs_sindex[k] for k, v in
        al_freqs_uselength_sindex.items()
    ) and len(al_freqs_uselength_sindex) == len(true_len_al_freqs_sindex))
    assert (all(
        v == true_al_freqs_sindex[k] for k, v
        in al_freqs_nolength_sindex.items()
    ) and len(al_freqs_nolength_sindex) == len(true_al_freqs_sindex))


def test_GetMaxAllele():
    # Test working example
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    true_al_max = 6
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max

    # Test example where alt=[]
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    true_al_max = 3
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max

    # Test example where there are no called samples
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    true_al_max = np.nan
    al_max = rec.GetMaxAllele()
    assert np.isnan(al_max)

    # Test working example with sample_index
    sindex = [0, 2, 5]
    true_al_max_sindex = 4
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    al_max_sindex = rec.GetMaxAllele(sample_index=sindex)
    assert al_max_sindex == true_al_max_sindex


def test_GetCalledSamples():
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    assert np.all(rec.GetCalledSamples() == [True] * 5 + [False])
    assert np.all(rec.GetCalledSamples(strict=False))

    # Test differences in ploidy
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    assert np.all(rec.GetCalledSamples(strict=True))

    # Test a true no call
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    assert np.all(~rec.GetCalledSamples())


def test_GetSamplePloidies():
    # All samples have the same ploidy
    # even a partial nocall
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    assert np.all(rec.GetSamplePloidies() == 2)

    # Test differences in ploidy
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    assert np.all(rec.GetSamplePloidies() == [2,2,2,3])

    # Test a no call, sample ploidy should not change
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    assert np.all(rec.GetSamplePloidies() == 2)


def test_GetCallRate():
    dummy_record = get_dummy_record()
    ref_allele = dummy_record.REF
    alt_alleles = dummy_record.ALT
    rec = trh.TRRecord(dummy_record, ref_allele, alt_alleles, "CAG", "", None)
    assert rec.GetCallRate() == pytest.approx(5/6)
    assert rec.GetCalledSamples(strict=False) == pytest.approx(1)

    # Test differences in ploidy
    rec = trh.TRRecord(get_triploid_record(), ref_allele, [], "CAG", "", None)
    assert rec.GetCalledSamples(strict=True) == pytest.approx(1)

    # Test a true no call
    rec = trh.TRRecord(get_nocall_record(), ref_allele, [], "CAG", "", None)
    assert rec.GetCalledSamples(strict=False) == pytest.approx(0)


#### Test TRRecordHarmonizer on different files ####


def reset_vcfs(vcfdir):
    global gangstr_vcf, hipstr_vcf, popstr_vcf, advntr_vcf, eh_vcf, snps_vcf
    gangstr_vcf = cyvcf2.VCF(os.path.join(vcfdir, "test_gangstr.vcf"))
    hipstr_vcf = cyvcf2.VCF(os.path.join(vcfdir, "test_hipstr.vcf"))
    popstr_vcf = cyvcf2.VCF(os.path.join(vcfdir, "test_popstr.vcf"))
    advntr_vcf = cyvcf2.VCF(os.path.join(vcfdir, "test_advntr.vcf"))
    eh_vcf = cyvcf2.VCF(os.path.join(vcfdir, "test_ExpansionHunter.vcf"))
    snps_vcf = cyvcf2.VCF(os.path.join(vcfdir, "snps.vcf"))


def test_multitype_vcf(vcfdir):
    with pytest.raises(TypeError):
        reader = cyvcf2.VCF(os.path.join(vcfdir, "test_multitype.vcf"))
        trh.TRRecordHarmonizer(reader)


def test_trh_init_and_type_infer(vcfdir):
    reset_vcfs(vcfdir)

    # Test example with a meaningless VCF type given
    with pytest.raises(ValueError):
        trh.TRRecordHarmonizer(gangstr_vcf, vcftype='unknownvcf')
    # Test example with unsupported VCF
    with pytest.raises(TypeError):
        trh.TRRecordHarmonizer(snps_vcf)

    # Test methods not passing in a vcftype
    with pytest.raises(ValueError):
        trh.MayHaveImpureRepeats('foo')
    with pytest.raises(ValueError):
        trh.HasLengthRefGenotype('foo')
    with pytest.raises(ValueError):
        trh.HasLengthAltGenotypes('foo')
    with pytest.raises(TypeError):
        trh.MayHaveImpureRepeats({}) # give a meaningless argument type

    # Test examples with correct preset VCF type
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf, vcftype='gangstr')
    assert gangstr_trh.vcftype == trh.VcfTypes.gangstr
    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf,
                                         vcftype=trh.VcfTypes.gangstr)
    assert gangstr_trh.vcftype == trh.VcfTypes.gangstr
    assert trh.InferVCFType(gangstr_vcf) == trh.VcfTypes.gangstr
    assert (not gangstr_trh.MayHaveImpureRepeats()
            and not trh.MayHaveImpureRepeats(trh.VcfTypes.gangstr))
    assert (not gangstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VcfTypes.gangstr))
    assert (not gangstr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VcfTypes.gangstr))
    assert not gangstr_trh.IsBeagleVCF()

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf, vcftype='hipstr')
    assert hipstr_trh.vcftype == trh.VcfTypes.hipstr
    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf,
                                        vcftype=trh.VcfTypes.hipstr)
    assert hipstr_trh.vcftype == trh.VcfTypes.hipstr
    assert trh.InferVCFType(hipstr_vcf) == trh.VcfTypes.hipstr
    assert (hipstr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VcfTypes.hipstr))
    assert (not hipstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VcfTypes.hipstr))
    assert (not hipstr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VcfTypes.hipstr))
    assert not hipstr_trh.IsBeagleVCF()

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf, vcftype='popstr')
    assert popstr_trh.vcftype == trh.VcfTypes.popstr
    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf,
                                        vcftype=trh.VcfTypes.popstr)
    assert popstr_trh.vcftype == trh.VcfTypes.popstr
    assert trh.InferVCFType(popstr_vcf) == trh.VcfTypes.popstr
    assert (popstr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VcfTypes.popstr))
    assert (not popstr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VcfTypes.popstr))
    assert (popstr_trh.HasLengthAltGenotypes()
            and trh.HasLengthAltGenotypes(trh.VcfTypes.popstr))

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf, vcftype='advntr')
    assert advntr_trh.vcftype == trh.VcfTypes.advntr
    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf,
                                        vcftype=trh.VcfTypes.advntr)
    assert advntr_trh.vcftype == trh.VcfTypes.advntr
    assert trh.InferVCFType(advntr_vcf) == trh.VcfTypes.advntr
    assert (advntr_trh.MayHaveImpureRepeats()
            and trh.MayHaveImpureRepeats(trh.VcfTypes.advntr))
    assert (not advntr_trh.HasLengthRefGenotype()
            and not trh.HasLengthRefGenotype(trh.VcfTypes.advntr))
    assert (not advntr_trh.HasLengthAltGenotypes()
            and not trh.HasLengthAltGenotypes(trh.VcfTypes.advntr))
    assert not advntr_trh.IsBeagleVCF()

    eh_trh = trh.TRRecordHarmonizer(eh_vcf, vcftype='eh')
    assert eh_trh.vcftype == trh.VcfTypes.eh
    eh_trh = trh.TRRecordHarmonizer(eh_vcf, vcftype=trh.VcfTypes.eh)
    assert eh_trh.vcftype == trh.VcfTypes.eh
    assert trh.InferVCFType(eh_vcf) == trh.VcfTypes.eh
    assert (not eh_trh.MayHaveImpureRepeats()
            and not trh.MayHaveImpureRepeats(trh.VcfTypes.eh))
    assert (eh_trh.HasLengthRefGenotype()
            and trh.HasLengthRefGenotype(trh.VcfTypes.eh))
    assert (eh_trh.HasLengthAltGenotypes()
            and trh.HasLengthAltGenotypes(trh.VcfTypes.eh))
    assert not eh_trh.IsBeagleVCF()

def test_imputed_vcf_types(vcfdir):
    imputed_gangstr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/gangstr_imputed.vcf.gz")), vcftype='gangstr')
    assert imputed_gangstr_trh.vcftype == trh.VcfTypes.gangstr
    assert imputed_gangstr_trh.IsBeagleVCF()
    assert not next(imputed_gangstr_trh).HasQualityScores()

    imputed_advntr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/advntr_imputed.vcf.gz")), vcftype='advntr')
    assert imputed_advntr_trh.vcftype == trh.VcfTypes.advntr
    assert imputed_advntr_trh.IsBeagleVCF()
    assert not next(imputed_advntr_trh).HasQualityScores()

    imputed_hipstr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/hipstr_imputed.vcf.gz")), vcftype='hipstr')
    assert imputed_hipstr_trh.vcftype == trh.VcfTypes.hipstr
    assert imputed_hipstr_trh.IsBeagleVCF()
    assert not next(imputed_hipstr_trh).HasQualityScores()

    imputed_eh_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/eh_imputed.vcf.gz")), vcftype='eh')
    assert imputed_eh_trh.vcftype == trh.VcfTypes.eh
    assert imputed_eh_trh.IsBeagleVCF()
    assert not next(imputed_eh_trh).HasQualityScores()

def test_missing_infos_imputed_vcfs_fail(vcfdir):
    missing_infos_imputed_gangstr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/gangstr_imputed_missing_infos.vcf.gz")), vcftype='gangstr')
    with pytest.raises(TypeError):
        next(missing_infos_imputed_gangstr_trh)

    missing_infos_imputed_advntr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/advntr_imputed_missing_infos.vcf.gz")), vcftype='advntr')
    with pytest.raises(TypeError):
        next(missing_infos_imputed_advntr_trh)

    missing_infos_imputed_hipstr_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/hipstr_imputed_missing_infos.vcf.gz")), vcftype='hipstr')
    with pytest.raises(TypeError):
        next(missing_infos_imputed_hipstr_trh)

    missing_infos_imputed_eh_trh = trh.TRRecordHarmonizer(cyvcf2.VCF(os.path.join(vcfdir, "beagle/eh_imputed_missing_infos.vcf.gz")), vcftype='eh')
    with pytest.raises(TypeError):
        next(missing_infos_imputed_eh_trh)

def test_string_or_vcftype(vcfdir):
    assert (trh.HasLengthAltGenotypes("gangstr")
            == trh.HasLengthAltGenotypes(trh.VcfTypes.gangstr))
    assert (trh.HasLengthRefGenotype("gangstr")
            == trh.HasLengthRefGenotype(trh.VcfTypes.gangstr))
    assert (trh.MayHaveImpureRepeats("gangstr")
            == trh.MayHaveImpureRepeats(trh.VcfTypes.gangstr))
    reset_vcfs(vcfdir)
    assert (trh.HarmonizeRecord("gangstr", next(gangstr_vcf)).GetMaxAllele()
            == len("tctgtctgtctg") / len("tctg"))
    assert (trh.HarmonizeRecord(trh.VcfTypes.gangstr,
                                next(gangstr_vcf)).GetMaxAllele()
            == len("aaaacaaaacaaaacaaaac") / len("aaaac"))


def all_types():
    for vcftype in trh.VcfTypes:
        yield vcftype
    yield "snps"


def get_vcf(vcftype):
    if vcftype == trh.VcfTypes.gangstr:
        return gangstr_vcf
    if vcftype == trh.VcfTypes.hipstr:
        return hipstr_vcf
    if vcftype == trh.VcfTypes.popstr:
        return popstr_vcf
    if vcftype == trh.VcfTypes.advntr:
        return advntr_vcf
    if vcftype == trh.VcfTypes.eh:
        return eh_vcf
    if vcftype == "snps":
        return snps_vcf
    raise ValueError("Unexpected vcftype")
        # TODO add Beagle


def test_wrong_vcftype(vcfdir):
    # an iterator that includes both tr caller types
    # and error file types
    for correct_type in trh.VcfTypes:
        reset_vcfs(vcfdir)
        for incorrect_type in all_types():
            if incorrect_type == correct_type:
                # make sure the incorrect_type is actually incorrect
                continue

            invcf = get_vcf(incorrect_type)
            with pytest.raises(TypeError):
                trh.TRRecordHarmonizer(invcf, vcftype=correct_type)

        reset_vcfs(vcfdir)
        for incorrect_type in all_types():
            if incorrect_type == correct_type:
                # make sure the incorrect_type is actually incorrect
                continue

            invcf = get_vcf(incorrect_type)
            record = next(invcf)
            with pytest.raises(TypeError):
                trh.HarmonizeRecord(correct_type, record)


def test_HarmonizeRecord(vcfdir):
    reset_vcfs(vcfdir)

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
    assert tr_rec3.ref_allele_length == 26
    assert tr_rec3.alt_allele_lengths == [24]

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
    assert tr_rec1.pos == tr_rec1.full_alleles_pos
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
    # record in the test file has flanking bp at the start of the ref. allele
    assert record.pos > record.full_alleles_pos
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


def assertFEquals(f1: float, f2: float):
    epsilon = 1e-6
    assert abs(f1 - f2) < epsilon

#def test_PHREDtoProb():
#    # pylint: disable=W0212
#    assertFEquals(trh._PHREDtoProb(0), 1)
#    assertFEquals(trh._PHREDtoProb(20), .01)
#    assertFEquals(trh._PHREDtoProb(2), 0.63095734448)

#def test_ConvertPLToQualityProb():
#    # pylint: disable=W0212
#    assertFEquals(trh._ConvertPLtoQualityProb([0]), 1)
#    assertFEquals(trh._ConvertPLtoQualityProb([10]), .1)
#    assertFEquals(trh._ConvertPLtoQualityProb([255, 10, 246]), .1)
#    assertFEquals(trh._ConvertPLtoQualityProb([10, 0, 10]), .8)
#    # confirm that PHRED scores of 0 don't drop below phred
#    # score of 1 despite our rebinning approach
#    assertFEquals(trh._ConvertPLtoQualityProb([0, 1, 1, 1]),
#                  trh._PHREDtoProb(1))

def _getVariantFromHarominzer(harmonizer, nvar=1):
    itr = iter(harmonizer)
    while nvar > 0:
        nvar -= 1
        var = next(itr)
    return var

def test_TRRecord_Quality(vcfdir):
    reset_vcfs(vcfdir)

    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf)
    assert gangstr_trh.HasQualityScore()
    var = _getVariantFromHarominzer(gangstr_trh)
    assert var.HasQualityScores()
    assert var.GetQualityScores()[0] == 0.999912

    gangstr_vcf_noqual = cyvcf2.VCF(
        os.path.join(vcfdir, "test_gangstr_noqual.vcf")
    )
    gangstr_trh_noqual = trh.TRRecordHarmonizer(gangstr_vcf_noqual)
    assert not gangstr_trh_noqual.HasQualityScore()
    var = _getVariantFromHarominzer(gangstr_trh_noqual)
    assert not var.HasQualityScores()
    with pytest.raises(TypeError):
        var.GetQualityScores()

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf)
    assert hipstr_trh.HasQualityScore()
    var = _getVariantFromHarominzer(hipstr_trh, nvar=18)
    assert var.HasQualityScores()
    assert var.GetQualityScores()[0] == 0.93

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf)
    assert not popstr_trh.HasQualityScore()
    var = _getVariantFromHarominzer(popstr_trh)
    assert not var.HasQualityScores()

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf)
    assert advntr_trh.HasQualityScore()
    var = _getVariantFromHarominzer(advntr_trh)
    assert var.HasQualityScores()
    assert var.GetQualityScores()[0] == 0.863

    eh_trh = trh.TRRecordHarmonizer(eh_vcf)
    assert not eh_trh.HasQualityScore()
    var = _getVariantFromHarominzer(eh_trh)
    assert not var.HasQualityScores()
    with pytest.raises(TypeError):
        var.GetQualityScores()
