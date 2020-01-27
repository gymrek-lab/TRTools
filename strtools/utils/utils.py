"""
Util functions for calculating summary STR statistics
and performing basic string operations on STR alleles.
"""
import itertools
import scipy.stats
import sys

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

    Examples
    --------
    >>> GetHeterozygosity({0:0.5, 1:0.5})
    0.5
    """
    exp_hom_frac = sum([val**2 for val in allele_freqs.values()])
    total_samples = sum(genotype_counts.values())
    num_hom = 0
    for gt in genotype_counts:
        if gt[0] == gt[1]: num_hom += genotype_counts[gt]
    return scipy.stats.binom_test(num_hom, n=total_samples, p=exp_hom_frac)

def GetHomopolymerRun(seq):
    seq = seq.upper()
    return max(len(list(y)) for (c,y) in itertools.groupby(seq))

nucToNumber={"A":0,"C":1,"G":2,"T":3}

def CheckMotif(motif):
    # If all consist of the same nucleotide, error
    if len(motif) > 1 and len(''.join(set(motif))) == 1: return False
    return True

def GetCanonicalMotif(repseq):
    """ Get canonical STR sequence, considering both strands """
    repseq = repseq.upper()
    repseq_f = getCanonicalMS(repseq)
    repseq_r = getCanonicalMS(reverseComplement(repseq))
    repseq = compareString(repseq_f, repseq_r)
    return repseq

def getCanonicalMS(repseq):
    """ Get canonical STR sequence """
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

def compareString(seq1,seq2):
    """ Compare two strings alphabetically """
    size = len(seq1)
    for i in range(size):
        if nucToNumber[seq1[i]] < nucToNumber[seq2[i]]:
            return seq1
        if nucToNumber[seq1[i]] > nucToNumber[seq2[i]]:
            return seq2
    return seq1

def reverseComplement(seq):
    """ Get the reverse complement of a nucleotide string """
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A": newseq += "T"
        if char == "G": newseq += "C"
        if char == "C": newseq += "G"
        if char == "T": newseq += "A"
    return newseq
