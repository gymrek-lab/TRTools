import os

import numpy as np
import pytest
import vcf

import trtools.utils.tr_harmonizer as trh

# pylint: disable=C0116,C0103

#### Test TRRecord using dummy info ####


# Set up dummy class with gt_alleles
class DummyVCFSample:
    def __init__(self,
                 gt_alleles,
                 called,
                 sample='',
                 quality_field=None,
                 quality_field_val=None):
        self.gt_alleles = gt_alleles
        self.called = called
        self.sample = sample
        self.formats = {}
        if quality_field is not None:
            self.formats[quality_field] = quality_field_val


class DummyVCFRecord:
    def __init__(self,
                 quality_field: str = None
                 ):
        self.samples = []
        self.POS = 42
        self.CHROM = '1984'
        self.formats = {}
        self.quality_field = quality_field

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


def test_unexpected_vcf_type():
    with pytest.raises(ValueError):
        trh._UnexpectedTypeError(trh.VcfTypes.gangstr)  # pylint: disable=W0212


def test_TRRecord_print():
    ref = "ABC"
    alt = ["DEF", "GHI"]
    motif = "foo"
    ID = "bar"
    record = trh.TRRecord(dummy_record1, ref, alt, motif, ID, "some_field")
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    record = trh.TRRecord(dummy_record1, ref, alt, motif, None, None)
    assert str(record) == "{}:{} {} {} {},{}".format(dummy_record1.CHROM,
                                                     dummy_record1.POS,
                                                     motif, ref, alt[0],
                                                     alt[1])

    record = trh.TRRecord(dummy_record1, ref, alt, motif, ID, None)
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    record = trh.TRRecord(dummy_record1, "B", ["E", "H"], motif, ID, None,
                          full_alleles=(ref, alt))
    assert str(record) == "{} {} {} {},{}".format(ID, motif, ref, alt[0],
                                                  alt[1])

    record = trh.TRRecord(dummy_record1, ref, None, motif, ID, None,
                          alt_allele_lengths=[3, 5.5])
    assert str(record) == "{} {} {} n_reps:3,n_reps:5.5".format(ID, motif, ref)

    record = trh.TRRecord(dummy_record1, None, None, motif, ID, None,
                          ref_allele_length=7,
                          alt_allele_lengths=[3, 5.5])
    assert str(record) == ("{} {} n_reps:7 n_reps:3,n_reps:5.5"
                           .format(ID, motif))


def test_TRRecord_iter():
    record = trh.TRRecord(dummy_record1, "ACG", ["A", "C", "G", "T"],
                          "FOO", "BAR", "some_field")
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
                     "some_field",
                     alt_allele_lengths=[4, 6])

    record = trh.TRRecord(dummy_record1, ref_allele, None, motif, ID,
                          "some_field",
                          alt_allele_lengths=[4, 5.5])
    assert record.alt_alleles == [motif * 4, motif * 5 + "F"]

    # ref allele
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID, None,
                     ref_allele_length=5)

    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, None, alt_alleles, motif, ID, None,
                     ref_allele_length=5)

    record = trh.TRRecord(dummy_record1, None, None, motif, ID, None,
                          ref_allele_length=5.5, alt_allele_lengths=[4, 5.5])
    assert record.ref_allele == motif * 5 + 'F'


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
    motif = 'FOO'
    ID = 'BAR'

    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, None, None, motif, ID, None,
                     full_alleles=(full_ref, full_alts))
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID, None,
                     full_alleles=(["CAGCAGCAQQQQQQQQQQQQQQQ"], full_alts))
    with pytest.raises(ValueError):
        bad_alts = [
            "CAGCAGCAQQQQQQQQQQQQQQQ",
            full_alts[1]
        ]
        trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID, None,
                     full_alleles=(ref_allele, bad_alts))

    record = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, motif, ID,
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
    # Test good example
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "", None)
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
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    print(rec)  # To test str function
    for sample in rec.vcfrecord:
        stringt = rec.GetStringGenotype(sample)
        lengt = rec.GetLengthGenotype(sample)
        assert(all([item == ref_allele for item in stringt]))
        assert(all([item == 3 for item in lengt]))
    # Test example with discrepancy between alt_alleles and genotypes given
    with pytest.raises(ValueError):
        trh.TRRecord(dummy_record1, ref_allele, [], "CAG", "", None)


def test_GetGenotypeCounts():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "", None)
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
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_gt_counts = {(3, 3, 3): 1, (3, 3): 3}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(
        v == true_len_gt_counts[k] for k, v in gt_counts_uselength.items()
    ) and len(gt_counts_uselength) == len(true_len_gt_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_gt_counts_slist = {}
    gt_counts_uselength_slist = \
        rec.GetGenotypeCounts(samplelist=['NonExistentSample'])
    assert (all(
        v == true_len_gt_counts_slist[k] for k, v
        in gt_counts_uselength_slist.items()
    ) and len(gt_counts_uselength_slist) == len(true_len_gt_counts_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "", None)
    true_len_gt_counts = {}
    gt_counts_uselength = rec.GetGenotypeCounts()
    assert (all(
        v == true_len_gt_counts[k] for k, v in gt_counts_uselength.items()
    ) and len(gt_counts_uselength) == len(true_len_gt_counts))


def test_GetAlleleCounts():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG", "CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "", None)
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
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_al_counts = {3: 9}
    al_counts_uselength = rec.GetAlleleCounts()
    assert (all(v == true_len_al_counts[k] for k,
                v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_al_counts_slist = {}
    al_counts_uselength_slist = rec.GetAlleleCounts(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_counts_slist[k] for k,v in al_counts_uselength_slist.items()) and len(al_counts_uselength_slist) == len(true_len_al_counts_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "", None)
    true_len_al_counts = {}
    al_counts_uselength = rec.GetAlleleCounts()
    assert (all(v == true_len_al_counts[k] for k,v in al_counts_uselength.items()) and len(al_counts_uselength) == len(true_len_al_counts))


def test_GetAlleleFreqs():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "", None)
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
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_al_freqs = {3: 1.0}
    al_freqs_uselength = rec.GetAlleleFreqs()
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))


    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_len_al_freqs_slist = {}
    al_freqs_uselength_slist = rec.GetAlleleFreqs(samplelist = ['NonExistentSample'])
    assert (all(v == true_len_al_freqs_slist[k] for k,v in al_freqs_uselength_slist.items()) and len(al_freqs_uselength_slist) == len(true_len_al_freqs_slist))

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "", None)
    true_len_al_freqs = {}
    al_freqs_uselength = rec.GetAlleleFreqs()
    assert (all(v == true_len_al_freqs[k] for k,v in al_freqs_uselength.items()) and len(al_freqs_uselength) == len(true_len_al_freqs))

def test_GetMaxAllele():
    # Test good example, no samplelist, uselength=True (default)
    ref_allele = "CAGCAGCAG"
    alt_alleles = ["CAGCAGCAGCAG","CAGCAGCAGCAGCAGCAG"]
    rec = trh.TRRecord(dummy_record1, ref_allele, alt_alleles, "CAG", "", None)
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
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_al_max = 3.0
    al_max = rec.GetMaxAllele()
    assert al_max == true_al_max

    # Test example with non of samples in samplelist in VCF
    rec = trh.TRRecord(dummy_record3, ref_allele, [], "CAG", "", None)
    true_al_max_slist = np.nan
    al_max_slist = rec.GetMaxAllele(samplelist = ['NonExistentSample'])
    assert np.isnan(al_max_slist) == True

    # Test example where that has one uncalled sample only
    rec = trh.TRRecord(dummy_record4, ref_allele, [], "CAG", "", None)
    true_al_max = np.nan
    al_max = rec.GetMaxAllele()
    assert np.isnan(al_max) == True


#### Test TRRecordHarmonizer on different files ####


def reset_vcfs(vcfdir):
    global gangstr_vcf, hipstr_vcf, popstr_vcf, advntr_vcf, eh_vcf, snps_vcf
    gangstr_vcf = vcf.Reader(filename=os.path.join(vcfdir, "test_gangstr.vcf"))
    hipstr_vcf = vcf.Reader(filename=os.path.join(vcfdir, "test_hipstr.vcf"))
    popstr_vcf = vcf.Reader(filename=os.path.join(vcfdir, "test_popstr.vcf"))
    advntr_vcf = vcf.Reader(filename=os.path.join(vcfdir, "test_advntr.vcf"))
    eh_vcf = vcf.Reader(filename=os.path.join(vcfdir, "test_ExpansionHunter.vcf"))
    snps_vcf = vcf.Reader(filename=os.path.join(vcfdir, "snps.vcf"))


def test_multitype_vcf(vcfdir):
    with pytest.raises(TypeError):
        reader = vcf.Reader(filename=os.path.join(vcfdir,
                                                  "test_multitype.vcf"))
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
                print(correct_type, incorrect_type)
                trh.TRRecordHarmonizer(invcf, vcftype=correct_type)

        reset_vcfs(vcfdir)
        for incorrect_type in all_types():
            if incorrect_type == correct_type:
                # make sure the incorrect_type is actually incorrect
                continue

            invcf = get_vcf(incorrect_type)
            record = next(invcf)
            with pytest.raises(TypeError):
                print(correct_type, incorrect_type)
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


def assertFEquals(f1: float, f2: float):
    epsilon = 1e-6
    assert abs(f1 - f2) < epsilon

def test_PHREDtoProb():
    # pylint: disable=W0212
    assertFEquals(trh._PHREDtoProb(0), 1)
    assertFEquals(trh._PHREDtoProb(20), .01)
    assertFEquals(trh._PHREDtoProb(2), 0.63095734448)

def test_ConvertPLToQualityProb():
    # pylint: disable=W0212
    assertFEquals(trh._ConvertPLtoQualityProb([0]), 1)
    assertFEquals(trh._ConvertPLtoQualityProb([10]), .1)
    assertFEquals(trh._ConvertPLtoQualityProb([255, 10, 246]), .1)
    assertFEquals(trh._ConvertPLtoQualityProb([10, 0, 10]), .8)
    # confirm that PHRED scores of 0 don't drop below phred
    # score of 1 despite our rebinning approach
    assertFEquals(trh._ConvertPLtoQualityProb([0, 1, 1, 1]),
                  trh._PHREDtoProb(1))

def _getVariantAndSampleFromHarominzer(harmonizer, nvar=1):
    itr = iter(harmonizer)
    while nvar > 0:
        nvar -= 1
        var = next(itr)
    samp = next(iter(var))
    return var, samp

def test_TRRecord_Quality(vcfdir):
    reset_vcfs(vcfdir)

    gangstr_trh = trh.TRRecordHarmonizer(gangstr_vcf)
    assert gangstr_trh.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(gangstr_trh)
    assert var.HasQualityScores()
    assert var.GetQualityScore(samp) == 0.999912

    gangstr_vcf_noqual = vcf.Reader(
        filename=os.path.join(vcfdir, "test_gangstr_noqual.vcf")
    )
    gangstr_trh_noqual = trh.TRRecordHarmonizer(gangstr_vcf_noqual)
    assert not gangstr_trh_noqual.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(gangstr_trh_noqual)
    assert not var.HasQualityScores()
    with pytest.raises(TypeError):
        var.GetQualityScore(samp)

    hipstr_trh = trh.TRRecordHarmonizer(hipstr_vcf)
    assert hipstr_trh.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(hipstr_trh, nvar=18)
    assert var.HasQualityScores()
    assert var.GetQualityScore(samp) == 0.93

    popstr_trh = trh.TRRecordHarmonizer(popstr_vcf)
    assert not popstr_trh.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(popstr_trh)
    assert not var.HasQualityScores()

    advntr_trh = trh.TRRecordHarmonizer(advntr_vcf)
    assert advntr_trh.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(advntr_trh)
    assert var.HasQualityScores()
    assert var.GetQualityScore(samp) == 0.863

    eh_trh = trh.TRRecordHarmonizer(eh_vcf)
    assert not eh_trh.HasQualityScore()
    var, samp = _getVariantAndSampleFromHarominzer(eh_trh)
    assert not var.HasQualityScores()
    with pytest.raises(TypeError):
        var.GetQualityScore(samp)

