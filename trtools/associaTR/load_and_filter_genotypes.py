#!/bin/env python3

"""
Read VCFs of TRs
apply filters to them for quality control
and then return yield them
"""

import sys
from typing import Dict, Optional, Union

import cyvcf2
import numpy as np

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

allele_len_precision = 2
allele_frequency_precision = 2
dosage_precision = 2
r2_precision = 2

def dict_str(d):
    out = '{'
    first = True
    for key in sorted(d.keys()):
        if not first:
            out +=', '
        first = False
        # make sure keys are quoted so that the resulting string
        # is valid JSON and can be parsed as a dictionary in javascript
        # using JSON.parse
        out += '{}: {}'.format(repr(str(key)), repr(d[key]))
    out += '}'
    return out.replace("'", '"').replace('(', '[').replace(')', ']').replace('nan', '"NaN"')

def clean_len_alleles(d):
    new_d = {}
    for key, val in d.items():
        new_key = round(key, allele_len_precision)
        if new_key not in new_d:
            new_d[new_key] = val
        else:
            new_d[new_key] += val
    return new_d

def clean_len_allele_pairs(d):
    new_d = {}
    for (k1, k2), val in d.items():
        new_key = (round(k1, allele_len_precision), round(k2, allele_len_precision))
        if new_key not in new_d:
            new_d[new_key] = val
        else:
            new_d[new_key] += val
    return new_d

def round_vals(d, precision):
    return {key: round(val, precision) for key, val in d.items()}

def load_trs(vcf_fname: str,
              samples: Union[np.ndarray, slice],
              region: Optional[str] = None,
              non_major_cutoff: float = 20,
              beagle_dosages: bool = False ,
              vcftype: Optional[str] = None,
              _imputed_ukb_strs_paper_period_check: bool = False
             ):
    """
    Iterate over a region returning genotypes at TR loci.

    First yield is a tuple of names of the fields in details.

    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
    vcf_fname:
    samples:
        A boolean array of length nsamples determining which samples are included
        (True) and which are not. Or a slice(None) object as a convenience for specifying
        all samples
    region:
        chr:start-end if not None
    non_major_cutoff : the cutoff value for the non-major allele count/dosage filter
    beagle_dosages: work with dosages from the AP format fields instead of the GT field
    vcftype: override for the vcftype, if not specified then will be inferred

    Yields
    ------
    gts: Union[Dict[float, np.ndarray], np.ndarray]
        If not beagle_dosages, then a float array n_samplesx2 of allele lengths

        If beagle_dosages, then a dictionary from unique length alleles to 2D arrays of size
        (n_samples, 2) which contain the dosages of those alleles for each haplotype
        Length dosage are measured in number of repeats.

        None if locus_filtered is not None
    unique_alleles: np.ndarray
        1D array of unique length alleles (measured in number of repeats),
        If beagle_dosages, this is the same length as the dosages dict
    chrom: str
        e.g. '13'
    pos: int
    called_samples_filter:
        after applying samples above, filters out samples which were not called for this locus
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        Currently the only filter is for non major allele count/dosage,
        filtered loci will look like 'non-major allele dosage/count<VAL'
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.
    """

    vcf = cyvcf2.VCF(vcf_fname)
    inferred_vcftype = trh.InferVCFType(vcf, vcftype if vcftype else 'auto')

    if region is not None:
        region_start = int(region.split(':')[1].split('-')[0])
        vcf = vcf(region)

    deets = [
        'motif',
        'period',
        'ref_len',
        'allele_frequency'
    ]
    if beagle_dosages :
        deets.extend([
            'dosage_estimated_r2_per_length_allele',
            'r2_length_dosages_vs_best_guess_lengths'
        ])
    yield deets

    first = True
    for record in vcf:
        if first and beagle_dosages and "AP1" not in record.FORMAT:
            print("--beagle-dosages specified, missing required field AP1 for the TR")
            if "GP" in record.FORMAT:
                print("We could support the GP field, but currently only support the AP fields")
            print("Erroring out")
            sys.exit(1)

        first = False

        if region is not None and record.POS < region_start:
            # records that overlap this region but started before this region
            # should be considered part of the pervious region and not returned here
            continue

        if _imputed_ukb_strs_paper_period_check and record.INFO.get('PERIOD') is None:
            # there are a few duplicate loci which I didn't handle
            # properly, this identifies and removes them
            continue

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=inferred_vcftype)

        if isinstance(samples, slice):
            assert samples == slice(None)
            called_samples_filter = trrecord.GetCalledSamples()
            curr_samples = called_samples_filter
        else:
            called_samples_filter = trrecord.GetCalledSamples()[samples]
            curr_samples = samples & trrecord.GetCalledSamples()
        
        n_samples = int(np.sum(curr_samples))

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        len_alleles = [round(allele_len, allele_len_precision) for allele_len in len_alleles]

        if not beagle_dosages:
            gts = trrecord.GetLengthGenotypes()[curr_samples, :-1]
            allele_frequency = clean_len_alleles(trrecord.GetAlleleFreqs(curr_samples))
        else:
            gts = {
                _len: np.zeros((n_samples, 2)) for _len in np.unique(len_alleles)
            }

            for p in (1, 2):
                ap = trrecord.format['AP{}'.format(p)]
                gts[len_alleles[0]][:, (p-1)] += \
                        np.maximum(0, 1 - np.sum(ap[curr_samples, :], axis=1))
                for i in range(ap.shape[1]):
                    gts[len_alleles[i+1]][:, (p-1)] += ap[curr_samples, i]

            allele_frequency = {
                _len: np.sum(gts[_len])/(2*n_samples) for _len in gts
            }

            # https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
            # Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.
            # appendix 1
            allele_dosage_r2 = {}

            best_guesses = trrecord.GetLengthGenotypes()[curr_samples, :-1]
            rounded_best_guesses = np.around(best_guesses, allele_len_precision)
            for length in len_alleles:
                # calculate allele dosage r**2 for this length
                if length in allele_dosage_r2:
                    continue

                calls = rounded_best_guesses == length

                allele_dosage_r2[length] = np.corrcoef(
                    calls.reshape(-1), gts[length].reshape(-1)
                )[0,1]**2

            length_r2 = np.corrcoef(
                best_guesses.flatten(),
                np.add.reduce([
                    len_*dosages for len_, dosages in gts.items()
                ]).flatten()
            )[0,1]**2

        locus_details = [
            trrecord.motif,
            str(len(trrecord.motif)),
            str(round(trrecord.ref_allele_length, allele_len_precision)),
            dict_str({key: '{:.2g}'.format(val) for key, val in allele_frequency.items()})
        ]
        if beagle_dosages:
            locus_details.extend([
                dict_str(round_vals(allele_dosage_r2, r2_precision)),
                str(round(length_r2, r2_precision))
            ])

        if len(allele_frequency) == 0:
            filter_reason = 'No called samples'
        elif len(allele_frequency) == 1:
            filter_reason = 'Only one called allele'
        else:
            af = list(allele_frequency.values())
            af.pop(np.argmax(af))
            if np.sum(af)*n_samples*2 < non_major_cutoff:
                filter_reason = 'non-major allele {}<{}'.format("dosage" if beagle_dosages else "count", non_major_cutoff)
            else:
                filter_reason = None

        if filter_reason:
            yield (
                None,
                np.unique(len_alleles),
                trrecord.chrom,
                trrecord.pos,
                called_samples_filter,
                filter_reason,
                locus_details
            )
        else:
            yield (
                gts,
                np.unique(len_alleles),
                trrecord.chrom,
                trrecord.pos,
                called_samples_filter,
                None,
                locus_details
            )

