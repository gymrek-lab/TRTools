"""
Util functions for calculating STR locus stats from VCF record
"""

import warnings
warnings.filterwarnings("ignore")

import itertools
import scipy.stats
import sys

def ERROR(msg):
    sys.stderr.write(msg.strip()+"\n")
    sys.exit(1)

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

def GetSTRHWE(record, samples=[], uselength=False):
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
            gt = map(int, sample["GB"].split("|"))
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
    return hwe_p
