"""
Util functions for calculating summary STR statistics
and performing basic string operations on STR alleles.
"""
import argparse
import itertools
import math
import os
from typing import Any, Dict, List, Optional

import cyvcf2
import numpy as np
import scipy.stats

import trtools.utils.common as common # pragma: no cover

nucToNumber={"A":0,"C":1,"G":2,"T":3}

def LoadSingleReader(
        vcf_loc: str,
        checkgz: bool = True) -> Optional[cyvcf2.VCF]:
    """
    Return a VCF reader

    Parameters
    ----------
    vcf_loc :
        The location of the VCF file to read
    checkgz:
        Check whether VCF file is gzipped and indexed,
        if not return None

    Returns
    -------
    reader : Optional[cyvcf2.VCF]
        The cyvcf2.VCF instance, or None if the VCF is not present
        or could not be opened
    """
    # check that vcf_loc is a file or file descriptor (ex: '/dev/stdin')
    if not os.path.exists(vcf_loc) or os.path.isdir(vcf_loc):
        common.WARNING("Could not find VCF file %s"%vcf_loc)
        return None
    if checkgz:
        if not vcf_loc.endswith(".vcf.gz") and not vcf_loc.endswith(".vcf.bgz"):
            common.WARNING("Make sure %s is bgzipped and indexed"%vcf_loc)
            return None
        if not os.path.isfile(vcf_loc+".tbi"):
            common.WARNING("Could not find VCF index %s.tbi"%vcf_loc)
            return None
    try:
        return cyvcf2.VCF(vcf_loc)
    except OSError:
        common.WARNING("Could not open VCF file %s. Is it really VCF?"%vcf_loc)
        return None

def LoadReaders(
        vcf_locs: List[str],
        checkgz: bool = True) -> Optional[List[cyvcf2.VCF]]:
    """
    Return a list of VCF readers

    Parameters
    ----------
    vcf_locs :
        A list of vcf locations
    checkgz:
        Check if each VCF file is gzipped and indexed

    Returns
    -------
    readers : Optional[List[cvycf2.VCF]]
        A list of VCF readers, or None if any of them could not be found.
        If checkgz, then also return None if any of the VCF readers
        were not gzipped and tabix indexed.
    """
    readers = []
    for f in vcf_locs:
        rdr = LoadSingleReader(f, checkgz)
        if rdr is None:
            return None
        readers.append(rdr)

    return readers

def GetContigs(vcf: cyvcf2.VCF) -> List[str]:
    '''
    Returns the contig IDs in the VCF.

    Parameters
    ----------
    vcf :
        The vcf to get contigs from

    Returns
    -------
    List[str] :
        A list of contig IDs
    '''
    contigs = []
    for header_line in vcf.header_iter():
        if header_line['HeaderType'].lower() == 'contig':
            contigs.append(header_line['ID'])
    return contigs

def ValidateAlleleFreqs(allele_freqs):
    r"""Check that the allele frequency distribution is valid.

    Allele frequencies must sum to 1.

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles are typically strings (sequences) or floats (num. repeats)

    Returns
    -------
    is_valid : bool
          Return True if the distribution is valid, else False

    Examples
    --------
    >>> ValidateAlleleFreqs({0:0.5, 1:0.5})
    True
    """
    if len(allele_freqs.keys()) == 0: return False
    return abs(1-sum(allele_freqs.values())) <= 0.001

def GetHeterozygosity(allele_freqs):
    r"""Compute heterozygosity of a locus

    Heterozygosity is defined as the probability
    that two randomly drawn alleles are different.

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles are typically strings (sequences) or floats (num. repeats)

    Returns
    -------
    heterozygosity : float
          The heterozygosity of the locus.
          If the allele frequencies dictionary is invalid, return np.nan

    Notes
    -----
    Heterzygosity is computed as:

    .. math:: H = 1-\sum_{i=1..n} p_i^2

    where `p_i` is the frequency of allele `i` and `n` is the number of alleles.

    Examples
    --------
    >>> GetHeterozygosity({0:0.5, 1:0.5})
    0.5
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    return 1-sum([freq**2 for freq in allele_freqs.values()])


def GetEntropy(allele_freqs: Dict[Any, float]) -> float:
    r"""Compute the (bit) entropy of a locus

    Entropy is defined as the
    `entropy <https://en.wikipedia.org/wiki/Information_content>_`
    of the distribution of allele frequencies.

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles are typically strings (sequences) or floats (num. repeats)

    Returns
    -------
    entropy: float
          The entropy of the locus.
          If the allele frequencies dictionary is invalid, return np.nan

    Notes
    -----
    Entropy is computed as:

    .. math:: E = -\sum_{i=1..n} -p_i*log_2(p_i)

    where `p_i` is the frequency of allele `i` and `n` is the number of alleles.

    Examples
    --------
    >>> GetEntropy({0:0.5, 1:0.5})
    1.0
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    return scipy.stats.entropy(list(x for x in allele_freqs.values()), base=2)


def GetMean(allele_freqs):
    r"""Compute the mean allele length

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles must be given in lengths (numbers, not strings)

    Returns
    -------
    mean: float
          Return mean if allele frequencies dictionary is valid

    Examples
    --------
    >>> GetMean({0:0.5, 1:0.5})
    0.5
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    return sum([key*allele_freqs[key] for key in allele_freqs])

def GetMode(allele_freqs):
    """
    Compute the mode allele length.

    If more than one allele has the maximum number of copies out of all
    alleles, choose one at random (but reproducibly)

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles must be given in lengths (numbers, not strings)

    Returns
    -------
    mode: float
          Return mode if allele frequencies dictionary is valid

    Examples
    --------
    >>> GetMode({0:0.1, 1:0.9})
    1
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    mode_freq = -1
    modes = set()
    for allele, freq in allele_freqs.items():
        if freq > mode_freq:
            modes = {allele}
            mode_freq = freq
        if freq == mode_freq:
            modes.add(allele)
    return min(modes)  # use min to make this arbitrary selection reproducible

def GetVariance(allele_freqs):
    r"""Compute the variance of the allele lengths

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles must be given in lengths (numbers, not strings)

    Returns
    -------
    variance: float
          Return variance if allele frequencies dictionary is valid
          np.nan otherwise.

    Examples
    --------
    >>> GetVariance({0:1})
    0
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    mean = GetMean(allele_freqs)
    return sum([allele_freqs[key]*(key-mean)**2 for key in allele_freqs.keys()])

def GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts):
    r"""Compute Hardy Weinberg p-value

    Tests whether the number of observed heterozygous vs.
    homozygous individuals is different than expected
    under Hardy Weinberg Equilibrium given the observed
    allele frequencies, based on a binomial test.

    Parameters
    ----------
    allele_freqs : dict of object: float
          Dictionary of allele frequencies for each allele.
          Alleles are typically strings (sequences) or floats (num. repeats)
    genotype_counts : dict of (object, object): int
          Dictionary of counts of each genotype. Genotypes are defined
          as tuples of alleles. Alleles must be the same as those given in
          allele_freqs

    Returns
    -------
    p-value : float
          The two-sided p-value returned by a binomial test (scipy.stats.binom_test)
          If the allele frequencies dictionary is invalid, return np.nan
          If genotype alleles not included in frequencies dictionary, return np.nan
    """
    if not ValidateAlleleFreqs(allele_freqs):
        return np.nan
    exp_hom_frac = sum([val**2 for val in allele_freqs.values()])
    total_samples = sum(genotype_counts.values())
    num_hom = 0
    for gt in genotype_counts:
        if gt[0] not in allele_freqs.keys():
            return np.nan
        if gt[1] not in allele_freqs.keys():
            return np.nan
        if gt[0] == gt[1]: num_hom += genotype_counts[gt]
    return scipy.stats.binom_test(num_hom, n=total_samples, p=exp_hom_frac)

def GetHomopolymerRun(seq):
    r"""Compute the maximum homopolymer run length in a sequence

    Parameters
    ----------
    seq : str
          String giving a sequence of nucleotides

    Returns
    -------
    runlength : int
          The length of the longest homopolymer run

    Examples
    --------
    >>> GetHomopolymerRun("AATAAAATAAAAAT")
    5
    """
    if len(seq) == 0: return 0
    seq = seq.upper()
    return max(len(list(y)) for (c,y) in itertools.groupby(seq))

def GetCanonicalMotif(repseq):
    r"""Get canonical STR sequence, considering both strands

    The canonical sequence is the first alphabetically
    out of all possible rotations on + and - strands
    of the sequence. e.g. "TG" canonical sequence is "AC".

    Parameters
    ----------
    repseq : str
          String giving a STR motif (repeat unit sequence)

    Returns
    -------
    canon : str
          The canonical sequence of the STR motif

    Examples
    --------
    >>> GetCanonicalMotif("TG")
    'AC'
    """
    repseq = repseq.upper()
    # Get canonical sequence of each strand
    repseq_f = GetCanonicalOneStrand(repseq)
    repseq_r = GetCanonicalOneStrand(ReverseComplement(repseq))
    # choose first seq alphabetically
    for i in range(len(repseq_f)):
        if nucToNumber[repseq_f[i]] < nucToNumber[repseq_r[i]]:
            return repseq_f
        if nucToNumber[repseq_r[i]] < nucToNumber[repseq_f[i]]:
            return repseq_r
    return repseq_f

def GetCanonicalOneStrand(repseq):
    r"""Get canonical STR sequence, considering one strand

    The canonical sequence is the first alphabetically
    out of all possible rotations. e.g. CAG -> AGC.

    Parameters
    ----------
    repseq : str
          String giving a STR motif (repeat unit sequence)

    Returns
    -------
    canon : str
          The canonical sequence of the STR motif

    Examples
    --------
    >>> GetCanonicalOneStrand("CAG")
    'AGC'
    """
    repseq = repseq.upper()
    size = len(repseq)
    canonical = repseq
    for i in range(size):
        newseq = repseq[size-i:]+repseq[0:size-i]
        for j in range(size):
            if nucToNumber[newseq[j]] < nucToNumber[canonical[j]]:
                canonical = newseq
            elif nucToNumber[newseq[j]] > nucToNumber[canonical[j]]:
                break
    return canonical

def ReverseComplement(seq):
    r"""Get reverse complement of a sequence.

    Converts everything to uppsercase.

    Parameters
    ----------
    seq : str
          String of nucleotides.

    Returns
    -------
    revcompseq : str
          Reverse complement of the input sequence.

    Examples
    --------
    >>> ReverseComplement("AGGCT")
    'AGCCT'
    """
    seq = seq.upper()
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A":
            newseq += "T"
        elif char == "G":
            newseq += "C"
        elif char == "C":
            newseq += "G"
        elif char == "T":
            newseq += "A"
        else: newseq += "N"
    return newseq

def InferRepeatSequence(seq, period):
    """
    TODO change to dynamic programming approach
    Infer the repeated sequence in a string

    Parameters
    ----------
    seq : str
        A string of nucleotides
    period : int
        Length of the repeat unit

    Returns
    -------
    repseq : str
        The inferred repeat unit (motif).
        If the input sequence is smaller than the period,
        repseq consists of all N's.

    Examples
    --------
    >>> InferRepeatSequence('ATATATAT', 2)
    'AT'
    """
    if period > len(seq):
        return "N"*period
    best_kmer = None
    best_copies = 0
    for offset in range(0, period):
        kmers = {}
        start_idx = 0
        while start_idx + period <= len(seq):
            kmer = seq[start_idx:(start_idx + period)]
            if kmer not in kmers:
                kmers[kmer] = 1
            else:
                kmers[kmer] += 1
            start_idx += period
            current_best_kmer = max(kmers, key = lambda k: kmers[k])
            current_best_copies = kmers[current_best_kmer]
            if current_best_copies > best_copies:
                best_kmer = current_best_kmer
                best_copies = current_best_copies
    return GetCanonicalOneStrand(best_kmer)

def LongestPerfectRepeat(seq, motif, check_reverse=True):
    r"""
    Determine the length (bp) of the longest 
    perfect repeat stretch

    Credit: function originally written by
    Helyaneh Ziaei-Jam

    Parameters
    ----------
    seq : str
       Repeat region sequence
    motif : str
       Repeat unit sequence
    check_reverse : bool (optional)
       If False, don't check reverse complement

    Returns
    -------
    max_match : int
       Number of bp of longest perfect stretch
    """
    max_matches = []
    seq = seq.upper()
    checkseqs = [seq]
    if check_reverse:
        strand_seq = ReverseComplement(seq)
        checkseqs.append(strand_seq)
    for ref_ in checkseqs:
        for mot in [motif, motif[::-1]]:
                i = 0
                match = 0
                max_match = 0  
                while True:
                    if i >= len(ref_):
                        break
                    for j in range(0,len(motif)):
                        k = i
                        while True:
                            while j < len(mot) and k < len(ref_) and ref_[k] == mot[j]:
                                k += 1
                                j += 1
                                match += 1
                            max_match = max(max_match, match)
                            if j == len(motif):
                                j = 0
                                i = k
                            else:
                                if j == len(motif) - 1:
                                    i += 1
                                match = 0
                                break
                        j = 0
                max_matches.append(max_match)
    return max(max_matches)

def FabricateAllele(motif, length):
    """
    Fabricate an allele with the given motif and length.

    Parameters
    ----------
    motif : str
        the motif to build the allele from
    length : float
        the number of times to copy the motif
        (noninteger implies partial repeats).
        This does NOT specify the desired length of the
        returned string.

    Return
    ------
    str
        the fabricated allele string

    Notes
    -----
    The allele is fabricated with the given motif orientation
    (e.g. motif = 'ACG' could produce 'ACGACGACG' but not 'CGACGACGA').
    Fabricated alleles will contain no impurities nor flanking base pairs.
    In the case of length being a noninteger float (because of partial
    repeats) and where it is unclear if the last nucleotide should be
    included in the fabricated repeat or not due to imprecision in the
    length float, the last nucleotide will always be left off (the
    length of the returned string will always be rounded down).
    """
    fab = math.floor(length) * motif
    idx = 0
    while (len(fab) + 1) / len(motif) < length:
        fab += motif[idx]
        idx += 1

    return fab


class ArgumentDefaultsHelpFormatter(argparse.HelpFormatter):  # pragma: no cover
    """
    Build a custom argument parser that works just like
    argparse.ArgumentDefaultsHelpFormatter except that it
    doesn't display None defaults. Everything below
    is copied from the python library source except for that.

    Help message formatter which adds default values to argument help.

    Only the name of this class is considered a public API. All the methods
    provided by the class are considered an implementation detail.
    """

    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if (action.default is not argparse.SUPPRESS and
                    action.default is not None):
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    help += ' (default: %(default)s)'
        return help


