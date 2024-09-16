"""
This file contains helper functions for
performing ABC (approximate Bayesian computation) in SISTR
"""

import os
import trtools.utils.common as common

def GetSummStats(allele_freqs_list, minfreq, numbins):
    abc_het = 1-sum([item**2 for item in allele_freqs_list])
    abc_common = len([item for item in allele_freqs_list if item >= minfreq])
    abc_bins = GetAlleleBins(allele_freqs_list, numbins)
    return {
        "het": abc_het,
        "common": abc_common,
        "bin": abc_bins
    }

def GetAlleleBins(allele_freqs_list, numbins):
    bins = [0]*numbins
    middle_index = int(len(allele_freqs_list)/2)

    # Everything below boundary_low is combined to a bin
    boundary_low = middle_index - int((numbins-1)/2)
    bins[0] = sum(allele_freqs_list[0:boundary_low+1])

    # Everything above boundary_high is combined to a bin
    boundary_high = middle_index + int((numbins-1)/2)
    bins[numbins-1] = sum(allele_freqs_list[boundary_high:len(allele_freqs_list)])

    # Fill in the rest between lower and upper boundary
    bins_index = 1
    for i in range(boundary_low+1, boundary_high):
        bins[bins_index] = allele_freqs_list[i]
        bins_index += 1
    return bins

def GetABCList(abc_file, minfreq, numbins):
    abc_list = []
    linenum = 0
    with open(abc_file, "r") as f:
        for line in f:
            if linenum == 0:
                header = line.strip().split("\t")
                freqs_col = header.index("freqs")
            else:
                info = line.strip().split("\t")
                s = float(info[0])
                allele_freqs_list = [float(freq) for freq in info[freqs_col].split(",")]
                summstats = GetSummStats(allele_freqs_list, minfreq, numbins)
                abc_list.append([s, summstats])
            linenum += 1
    return abc_list

class SistrABC:
    def __init__(self,
            sistr_index=None,
            minfreq=0.05,
            numbins=5
        ):
        self.minfreq = minfreq
        self.numbins = numbins
        self.LoadABC(sistr_index)

    def LoadABC(self, sistr_index):
        self.ABC_tables = {} # period->opt->table
        abc_files = [item for item in os.listdir(sistr_index) if \
            not item.endswith("_freqs.txt") and not item.endswith("_zero_freqs.txt") and \
            item.endswith(".txt")]
        for abcf in abc_files:
            period = int(abcf.replace(".txt","").split("_")[-2])
            opt_allele = int(abcf.replace(".txt","").split("_")[-1])
            common.MSG(f"Loading ABC index for period={period} opt={opt_allele}", debug=True)
            if period not in self.ABC_tables:
                self.ABC_tables[period] = {}
            self.ABC_tables[period][opt_allele] = GetABCList(os.path.join(sistr_index, abcf), \
                self.minfreq, self.numbins)

    def CheckForModel(self, period, opt_allele):
        if period not in self.ABC_tables.keys():
            return False
        if opt_allele not in self.ABC_tables[period].keys():
            return False
        return True
    