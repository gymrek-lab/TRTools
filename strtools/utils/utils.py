"""
Util functions for calculating STR locus stats from VCF record
"""

import warnings
warnings.filterwarnings("ignore")

import itertools
import scipy.stats
import sys

def GetHeterozygosity(allele_freqs):
    return 1-sum([freq**2 for freq in allele_freqs.values()])

def GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts):
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
