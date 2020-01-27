"""
Util functions for calculating summary STR statistics
and performing basic string operations on STR alleles.
"""
import itertools
import numpy as np
import scipy.stats
import sys

nucToNumber={"A":0,"C":1,"G":2,"T":3}

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
    return sum(allele_freqs.values()) == 1

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
    "AC"
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
    "AGC"
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
    "AGCCT"
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
        The inferred repeat unit (motif)

    Examples
    --------
    >>> InferRepeatSequence('ATATATAT')
    'AT'
    """
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
