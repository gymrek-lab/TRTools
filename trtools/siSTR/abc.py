"""
This file contains helper functions for
performing ABC (approximate Bayesian computation) in SISTR
"""

import numpy as np
import os
import trtools.utils.common as common
from scipy.stats.distributions import chi2

PRECISION = 7 # for rounding LRT numbers

# Survival function for mixture distribution
def SF(x):
    if x > 0:
        return 0
    if x <= 0:
        return 1

def GetSummStats(allele_freqs_list, minfreq, numbins):
    abc_het = 1-sum([item**2 for item in allele_freqs_list])
    abc_common = len([item for item in allele_freqs_list if item >= minfreq])
    abc_bins = GetAlleleBins(allele_freqs_list, numbins)
    return {
        "het": abc_het,
        "common": abc_common,
        "bin": abc_bins
    }

def GetSummStatsFromLRTLine(line, minfreq, numbins):
    hets = []
    commons = []
    bins = []
    freqs_list = line.split("\t")[1].split(";")
    for fl in freqs_list:
        fl_list = [float(item) for item in fl.split(",")]
        summstats = GetSummStats(fl_list, minfreq, numbins)
        hets.append(summstats["het"])
        commons.append(summstats["common"])
        bins.append(summstats["bin"])
    return {"het": hets, "common": commons, "bin": bins}

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
            abc_index,
            lrt_index,
            minfreq=0.05,
            numbins=5,
            eps_het_numerator=0.005,
            eps_het_denominator=3.0,
            eps_bins=0.3,
            min_abc_acceptance=10
        ):
        self.minfreq = minfreq
        self.numbins = numbins
        self.eps_het_numerator = eps_het_numerator
        self.eps_het_denominator = eps_het_denominator
        self.eps_bins = eps_bins
        self.min_abc_acceptance = min_abc_acceptance

        # Which summary stats to use?
        # For now only use het, bins as in
        # SISTR_v1 https://github.com/BonnieCSE/SISTR/blob/master/sistr/SISTR_v1.py#L32
        self.use_het = True
        self.use_common = False
        self.use_bins = True

        # These get set when loading tables
        self.LRT_tables = {}
        self.LRT_tables_ZERO = {}
        self.ABC_tables = {}

        # Load the tables
        self.LoadABC(abc_index)
        self.LoadLRT(lrt_index)

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

    def GetEpsilons(self, obs_summ_stats):
        eps_het = (obs_summ_stats["het"]-self.eps_het_numerator)/self.eps_het_denominator
        eps_common = (obs_summ_stats["common"])+1
        eps_bins = self.eps_bins
        return eps_het, eps_common, eps_bins

    def ABCCompare(self, obs_summ_stats, sim_summ_stats,
            eps_het, eps_common, eps_bins):
        # If any of the stats we are using fails, continue and don't add
        if self.use_het:
            if abs(obs_summ_stats["het"]-sim_summ_stats["het"]) >= eps_het:
                return False
        if self.use_common:
            if abs(obs_summ_stats["common"]-sim_summ_stats["common"]) >= eps_common:
                return False
        if self.use_bins:
            if sum([abs(obs_summ_stats["bin"][i]-sim_summ_stats["bin"][i]) for i \
                    in range(len(obs_summ_stats["bin"]))]) >= eps_bins:
                return False
        return True

    def RunABC(self, obs_summ_stats, period, opt_allele):
        eps_het, eps_common, eps_bins = self.GetEpsilons(obs_summ_stats)
        s_accepted = []
        num_tested = 0
        for sim in self.ABC_tables[period][opt_allele]:
            sval, sim_summ_stats = sim
            if self.ABCCompare(obs_summ_stats, sim_summ_stats, \
                eps_het, eps_common, eps_bins):
                s_accepted.append(sval)
            num_tested += 1
        abc_results = {"passed": len(s_accepted)>=self.min_abc_acceptance}
       
        # Quit if too few accepted
        if not abc_results["passed"]: return abc_results

        # Compile results and return
        abc_results["median_s"] = round(np.median(s_accepted), PRECISION)
        abc_results["lower_bound"] = round(np.percentile(s_accepted, 2.5), PRECISION)
        abc_results["upper_bound"] = round(np.percentile(s_accepted, 97.5), PRECISION)
        abc_results["num_accepted"] = len(s_accepted)
        abc_results["num_tested"] = num_tested
        return abc_results

    def LoadLRT(self, sistr_index):
        self.LRT_tables_ZERO = {} # period->opt->table for s=0
        self.LRT_tables = {} # period->opt->s->table for s>0
        lrt_zero_files = [item for item in os.listdir(sistr_index) if \
            item.endswith("_zero_freqs.txt")]
        lrt_files = [item for item in os.listdir(sistr_index) if \
            item.endswith("_freqs.txt") and not item.endswith("_zero_freqs.txt")]

        # Load LRT summary statistics for s = 0
        for lzf in lrt_zero_files:
            period = int(lzf.replace("_zero_freqs.txt","").split("_")[-2])
            opt_allele = int(lzf.replace("_zero_freqs.txt","").split("_")[-1])
            if period not in self.LRT_tables_ZERO.keys():
                self.LRT_tables_ZERO[period] = {}
            self.LRT_tables_ZERO[period][opt_allele] = \
                GetSummStatsFromLRTLine(open(os.path.join(sistr_index, lzf), "r").readlines()[1].strip(), \
                    self.minfreq, self.numbins)

        # Load LRT summary statistics for other s values
        for lf in lrt_files:
            period = int(lf.replace("_freqs.txt","").split("_")[-2])
            opt_allele = int(lf.replace("_freqs.txt","").split("_")[-1])
            if period not in self.LRT_tables.keys():
                self.LRT_tables[period] = {}
            self.LRT_tables[period][opt_allele] = {}
            s_list_available = set()
            with open(os.path.join(sistr_index, lf), "r") as f:
                for line in f:
                    if line.startswith("s"): continue # header
                    sval = float(line.strip().split("\t")[0])
                    s_list_available.add(sval)
                    self.LRT_tables[period][opt_allele][sval] = \
                        GetSummStatsFromLRTLine(line.strip(), self.minfreq, self.numbins)
                    self.LRT_tables[period][opt_allele]["s_avail"] = s_list_available

    def GetClosestSValue(self, s_ABC, period, opt_allele):
        min_dist = 100000000
        nearest_s = -2
        for elem in self.LRT_tables[period][opt_allele]["s_avail"]:
            dist = abs(s_ABC - elem)
            if dist < min_dist:
                min_dist = dist
                nearest_s = elem
        return nearest_s

    def GetLikelihood(self, lrt_data, obs_summ_stats):
        eps_het, eps_common, eps_bins = self.GetEpsilons(obs_summ_stats)
        num_accepted = 0
        num_sims = 0
        for i in range(len(lrt_data["het"])):
            sim_summ_stats = {
                "het": lrt_data["het"][i],
                "common": lrt_data["common"][i],
                "bin": lrt_data["bin"][i]
            }
            if self.ABCCompare(obs_summ_stats, sim_summ_stats, \
                eps_het, eps_common, eps_bins):
                num_accepted += 1
            num_sims += 1
        return (num_accepted+1)/num_sims

    def LikelihoodRatioTest(self, s_ABC, obs_summ_stats, period, opt_allele):
        s_ABC_round = self.GetClosestSValue(s_ABC, period, opt_allele)

        # Get likelihood for s=0 and s=ABC_s
        likelihood_0 = self.GetLikelihood(self.LRT_tables_ZERO[period][opt_allele], obs_summ_stats)
        likelihood_s_ABC = self.GetLikelihood(self.LRT_tables[period][opt_allele][s_ABC_round], obs_summ_stats)

        # Get statistics and return
        LR = likelihood_0/likelihood_s_ABC
        LogLR = -2*np.log(LR)
        lrt_results = {
            "likelihood_0": round(likelihood_0, PRECISION),
            "likelihood_s_ABC": round(likelihood_s_ABC, PRECISION),
            "LR": round(LR, PRECISION),
            "LogLR": round(LogLR, PRECISION),
            "pval": round(0.5*SF(LogLR) + 0.5*chi2.sf(LogLR, 1), PRECISION)
        }
        return lrt_results
