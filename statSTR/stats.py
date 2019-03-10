"""
Locus-level and Call-level VCF filters
"""

import strtools.utils.common as common
import strtools.utils.utils as utils

from pybedtools import BedTool
import os
import sys
import vcf.filters

###################################
# Locus level filters
###################################

class Stat_Thresh(vcf.filters.Base):
    'Filter site by call rate'

    name = 'THRESH'

    def __init__(self, thresh):
        self.threshold = thresh

    def __call__(self, record):
        if record.call_rate < self.threshold: return record.call_rate
        else: return None

class Stat_mean(vcf.filters.Base):
    'Filter site by HWE'

    name = 'MEAN'

    def __init__(self, min_locus_hwep, uselength=False):
        self.threshold = min_locus_hwep
        self.uselength = uselength

    def __call__(self, record):
        hwep = utils.GetSTRHWE(record, uselength=self.uselength)
        if hwep < self.threshold: return hwep
        else: return None

class Stat_sd(vcf.filters.Base):
    'Filter sites with low heterozygosity'

    name = 'SD'

    def __init__(self, min_locus_het, uselength=False):
        self.threshold = min_locus_het
        self.uselength = uselength

    def __call__(self, record):
        if self.uselength:
            het = utils.GetLengthHet(record)
        else: het = record.heterozygosity
        if het < self.threshold:
            return het
        return None

def create_region_filter(name, filename):
    class Filter_Regions(vcf.filters.Base):
        'Filter regions from file'
        def __init__(self, name, filename):
            self.threshold = ""
            self.name = name
            self.LoadRegions(filename)
        def LoadRegions(self, filename):
            if not os.path.exists(filename):
                common.ERROR("%s not found"%filename)
            self.regions = BedTool(filename)
            if not self.regions._tabixed():
                sys.stderr.write("Creating tabix index for %s\n"%filename)
                self.regions.tabix(force=True)
        def __call__(self, record):
            interval = "%s:%s-%s"%(record.CHROM, record.POS, record.POS+len(record.REF))
            if "chr" in interval:
                interval2 = interval.replace("chr","")
            else: interval2 = "chr"+interval
            try:
                tb1 = self.regions.tabix_intervals(interval)
                if tb1.count() > 0: return self.name
            except ValueError: pass
            try:
                tb2 = self.regions.tabix_intervals(interval2)
                if tb2.count() > 0: return self.name
            except ValueError: pass
            return None
    f = Filter_Regions(name, filename)
    return Filter_Regions(name, filename)

###################################
# Call level filters
###################################

class Reason:
    name = ""
    def __init__(self):
        pass
    def GetReason(self):
        return self.name

class LowCallDepth(Reason):
    name = "LowCallDepth"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if sample["DP"] < self.threshold: return sample["DP"]
        else: return None

class CallMinSuppReads(Reason):
    name = "MinSuppReads"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if sample["ALLREADS"] is None: return 0
        delim = "|"
        if "/" in sample["GB"]: delim = "/"
        gb = map(int, sample["GB"].split(delim))
        allreads = sample["ALLREADS"].split(";")
        r1 = 0
        r2 = 0
        for item in allreads:
            allele, readcount = map(int, item.split("|"))
            if allele == gb[0]: r1 += readcount
            if allele == gb[1]: r2 += readcount
        min_read_count = min([r1, r2])
        if min_read_count < self.threshold: return min_read_count
        else: return None

###############################
# GangSTR filters
###############################

class ProbHom(Reason):
    name = "ProbHom"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob hom expansion
        if sample["QEXP"][2] < self.threshold: return sample["QEXP"][2]
        else: return None

class ProbHet(Reason):
    name = "ProbHet"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob het expansion
        if sample["QEXP"][1] < self.threshold: return sample["QEXP"][1]
        else: return None


class BadCI(Reason):
    name = "BadCI"
    def __init__(self):
        pass
    def __call__(self, sample):
        #### If ML estimate outside of CI
        ml = [int(item) for item in sample["REPCN"]]
        ci = sample["REPCI"].split(",")
        for i in range(len(ml)):
            ml_est = ml[i]
            ci_low = int(ci[i].split("-")[0])
            ci_high = int(ci[i].split("-")[1])
            if ml_est<ci_low or ml_est>ci_high: return ml_est
        return None
