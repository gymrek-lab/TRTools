#!/bin/env python3

"""
Read VCFs of STRs or SNPs,
apply filters to them for quality control
and then return yield them
"""

from typing import Optional, Set, Union

import cyvcf2
import numpy as np

import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

allele_len_precision = 2
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
        out += f'{repr(str(key))}: {repr(d[key])}'
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

def load_strs(vcf_fname: str,
              region: Optional[str] = None,
              samples: Union[np.ndarray, slice],
              details: bool = True,
              hardcalls = False,
              _imputed_ukb_strs_paper_period_check: bool = False
             ):
    """
    Iterate over a region returning genotypes at STR loci.

    First yield is a tuple of names of the fields in details.
    Every subsequent yield is described in the yields section below.

    Parameters
    ----------
    vcf_fname:
    samples:
        A boolean array of length nsamples determining which samples are included
        (True) and which are not
    region:
        chr:start-end if not None

    Yields
    ------
    dosages: Dict[float, np.ndarray]
        A dictionary from unique length alleles to 2D arrays of size (n_samples, 2)
        which contain the dosages of those alleles for each haplotype
        Length dosage are measured in number of repeats.

        None if locus_filtered is not None

        If hardcalls, then instead of that just an array nx2 of length alleles
    unique_alleles: np.ndarray
        Array of unique length alleles (measured in number of repeats),
        same length as the dosages dict
    chrom: str
        e.g. '13'
    pos: int
    locus_filtered:
        None if the locus is not filtered, otherwise
        a string explaining why.
        'MAC<20' if the minor allele dosage is less than 20
        after sample subsetting, per plink's standard.

        None if hardcalls.
    locus_details:
        tuple of strings with the same length as the first yield
        with the corresponding order.

        None if hardcalls.

    Notes
    -----
    Hardcalls mentioned in the locus details are phased hardcalls, and in some
    corner cases will not correspond to the maximum likelihood unphased allele.
    """

    # TODO make hardcalls work
    # TODO rename hardcalls as length genotypes?

    vcf = cyvcf2.VCF(vcf_fname)

    if region is not None:
        region_start = int(region.split(':')[1].split('-')[0])
        vcf = cyvcf2(region)

    if details:
        yield (
            'motif',
            'period',
            'ref_len',
            #'total_per_allele_dosages',
            #'total_hardcall_alleles',
            #'total_hardcall_genotypes',
            'subset_total_per_allele_dosages', # TODO rename these
            'subset_total_hardcall_alleles', # TODO only one of these
            #'subset_total_hardcall_genotypes',
            #'subset_het',
            #'subset_entropy',
            #'subset_HWEP',
            'subset_allele_dosage_r2'
            # TODO more imputation metrics
        )

    for record in vcf:
        if region is not None and record.POS < region_start:
            # records that overlap this region but started before this region
            # should be considered part of the pervious region and not returned here
            continue

        if _imputed_ukb_strs_paper_period_check and record.INFO.get('PERIOD') is None:
            # there are a few duplicate loci which I didn't handle
            # properly, this identifies and removes them
            continue

        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype='hipstr')

        len_alleles = [trrecord.ref_allele_length] + trrecord.alt_allele_lengths
        len_alleles = [round(allele_len, allele_len_precision) for allele_len in len_alleles]

        #if hardcalls:
        #    yield (trrecord.GetLengthGenotypes()[samples, :-1], np.unique(len_alleles), trrecord.chrom, trrecord.pos, None, None)

        if isinstance(samples, slice):
            assert samples == slice(None)
            n_subset_samples = trrecord.GetNumSamples()
        else:
            n_subset_samples = int(np.sum(samples))

        subset_dosage_gts = {
            _len: np.zeros((n_subset_samples, 2)) for _len in np.unique(len_alleles)
        }

        for p in (1, 2):
            # TODO this isn't right - these best guess calls are only best guess if there's only one
            # allele per length, doesn't take into account prob splitting over imperfections
            # todo genotype dosages
            ap = trrecord.format[f'AP{p}']
            subset_dosage_gts[len_alleles[0]][:, (p-1)] += \
                    np.maximum(0, 1 - np.sum(ap[samples, :], axis=1))
            for i in range(ap.shape[1]):
                subset_dosage_gts[len_alleles[i+1]][:, (p-1)] += ap[samples, i]

        subset_total_dosages = {
            _len: np.sum(subset_dosage_gts[_len]) for _len in subset_dosage_gts
        }

        if details:
            # Is there an issue here? See TODO above
            #subset_total_hardcall_alleles = clean_len_alleles(trrecord.GetAlleleCounts(samples))
            #subset_total_hardcall_genotypes = clean_len_allele_pairs(trrecord.GetGenotypeCounts(samples))
            subset_hardcall_allele_freqs = clean_len_alleles(trrecord.GetAlleleFreqs(samples))

#            subset_het = utils.GetHeterozygosity(subset_hardcall_allele_freqs)
#            subset_entropy = utils.GetEntropy(subset_hardcall_allele_freqs)
#            subset_hwep = utils.GetHardyWeinbergBinomialTest(
#                subset_hardcall_allele_freqs,
#                subset_total_hardcall_genotypes
#            )

            # https://www.cell.com/ajhg/fulltext/S0002-9297(09)00012-3#app1
            # Browning, Brian L., and Sharon R. Browning. "A unified approach to genotype imputation and haplotype-phase inference for large data sets of trios and unrelated individuals." The American Journal of Human Genetics 84.2 (2009): 210-223.
            # appendix 1
            subset_allele_dosage_r2 = {}

            subset_hardcalls = np.around(trrecord.GetLengthGenotypes()[samples, :-1], allele_len_precision)
            for length in len_alleles:
                # calculate allele dosage r**2 for this length
                if length in subset_allele_dosage_r2:
                    continue

                calls = subset_hardcalls == length

                subset_allele_dosage_r2[length] = np.corrcoef(
                    calls.reshape(-1), subset_dosage_gts[length].reshape(-1)
                )[0,1]**2

            locus_details = (
                trrecord.motif,
                str(len(trrecord.motif)),
                str(round(trrecord.ref_allele_length, allele_len_precision)),
                #dict_str(round_vals(total_dosages, dosage_precision)),
                #dict_str(total_hardcall_alleles),
                #dict_str(total_hardcall_genotypes),
                dict_str(round_vals(subset_total_dosages, dosage_precision)),
                dict_str(subset_total_hardcall_freqs), #dict_str(subset_total_hardcall_alleles),
                #dict_str(subset_total_hardcall_genotypes),
                #str(subset_het),
                #str(subset_entropy),
                #str(subset_hwep),
                dict_str(round_vals(subset_allele_dosage_r2, r2_precision))
            )
        else:
            locus_details = None

        mac = list(subset_total_dosages.values())
        mac.pop(np.argmax(mac))

        if np.sum(mac) < 20:
            yield (
                None,
                np.unique(len_alleles),
                trrecord.chrom,
                trrecord.pos,
                'MAC<20',
                locus_details
            )
            continue

        yield (
            subset_dosage_gts,
            np.unique(len_alleles),
            trrecord.chrom,
            trrecord.pos,
            None,
            locus_details
        )

