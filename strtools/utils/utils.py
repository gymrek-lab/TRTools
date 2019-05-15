"""
Util functions for calculating STR locus stats from VCF record
"""

import warnings
warnings.filterwarnings("ignore")

import itertools
import scipy.stats
import sys

def GetLengthHet(record):
    len_to_count = {}
    for sample in record:
        if sample["GT"] is None or sample["GT"] == "." or sample["GT"] == "./.": continue
        alleles = sample.gt_bases.split(sample.gt_phase_char())
        l1 = len(alleles[0])
        l2 = len(alleles[1])
        len_to_count[l1] = len_to_count.get(l1, 0) + 1
        len_to_count[l2] = len_to_count.get(l2, 0) + 1
    total = sum(len_to_count.values())
    x = 0
    for al in len_to_count.keys():
        x += (len_to_count[al]*1.0/total)**2
    return 1-x

def GetHomopolymerRun(seq):
    seq = seq.upper()
    return max(len(list(y)) for (c,y) in itertools.groupby(seq))

def GetSTRHWE(record, samples=[], uselength=False, het_output=False):
    hwe_p = 0
    het = 0
    # Get genotypes, allele frequencies
    allele_counts = {}
    obs_het = 0
    obs_hom = 0
    total = 0
    for sample in record:
        if len(samples)>0 and sample.sample not in samples: continue
        if sample["GT"] == "." or sample["GT"] == "./." or sample["GT"] is None: continue
        if uselength:
            try:
                gt = map(int, sample["GB"].split("|")) #HipSTR
            except:
                gt = sample["REPCN"] #GangSTR
        else:
            gt = sample.gt_alleles
        if gt[0] == gt[1]: obs_hom += 1
        else:
            obs_het += 1
        total += 1
        for al in gt:
            allele_counts[al] = allele_counts.get(al, 0) + 1
    # Get Allele frequencies
    allele_freqs = {}
    for key in allele_counts.keys():
        allele_freqs[key] = allele_counts[key]*1.0/sum(allele_counts.values())
    # Get expected num homs/hets
    exp_hom_frac = 0
    for al in allele_freqs.keys():
        exp_hom_frac += allele_freqs[al]**2

    # Binomial test for HWE
    hwe_p = scipy.stats.binom_test(obs_het, n=obs_het+obs_hom, p=1-exp_hom_frac)
    if het_output: # Extended output is HWE P and obs_het
        return hwe_p, obs_het
    else: # Standard output is HWE P
        return hwe_p

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
