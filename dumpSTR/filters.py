"""
Locus-level and Call-level VCF filters
"""

from pybedtools import BedTool
import os
import sys
import vcf.filters

if __name__ == "filters":
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import common
    import tr_harmonizer as trh
    import utils
else:
    import strtools.utils.common as common
    import strtools.utils.tr_harmonizer as trh
    import strtools.utils.utils as utils

###################################
# Locus level filters
###################################

class Filter_MinLocusCallrate(vcf.filters.Base):
    'Filter site by call rate'

    name = 'CALLRATE'

    def __init__(self, min_locus_callrate):
        self.threshold = min_locus_callrate

    def __call__(self, record):
        if record.call_rate < self.threshold: return record.call_rate
        else: return None

class Filter_MinLocusHWEP(vcf.filters.Base):
    'Filter site by HWE'

    name = 'HWE'

    def __init__(self, min_locus_hwep, uselength=False):
        self.threshold = min_locus_hwep
        self.uselength = uselength

    def __call__(self, record):
        trrecord = TRRecord(record, record.REF, record.ALT, "", "") # TODO handle different VCF types
        allele_freqs = trrecord.GetAlleleFreqs(uselength=self.uselength)
        genotype_counts = trrecord.GetGenotypeCounts(uselength=self.uselength)
        hwep = utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)
        if hwep < self.threshold: return hwep
        else: return None

class Filter_MinLocusHet(vcf.filters.Base):
    'Filter sites with low heterozygosity'

    name = 'HETLOW'

    def __init__(self, min_locus_het, uselength=False):
        self.threshold = min_locus_het
        self.uselength = uselength

    def __call__(self, record):
        trrecord = TRRecord(record, record.REF, record.ALT, "", "") # TODO handle different VCF types
        het = utils.GetHeterozygosity(trrecord.GetAlleleFreqs(uselength=self.uselength))
        if het < self.threshold:
            return het
        return None

class Filter_MaxLocusHet(vcf.filters.Base):
    'Filter sites with high heterozygosity'

    name = 'HETHIGH'

    def __init__(self, max_locus_het, uselength=False):
        self.threshold = max_locus_het
        self.uselength = uselength

    def __call__(self, record):
        trrecord = TRRecord(record, record.REF, record.ALT, "", "") # TODO handle different VCF types
        het = utils.GetHeterozygosity(trrecord.GetAlleleFreqs(uselength=self.uselength))
        if het > self.threshold:
            return het
        return None

class Filter_LocusHrun(vcf.filters.Base):
    'Filter penta/hexa sites with long homopolymer runs'
    
    name = 'HRUN'

    def __init__(self):
        self.threshold = 0 # e.g. pentamers with hruns of 5+threshold, hexa with 6+threshold

    def __call__(self, record):
        hrun = utils.GetHomopolymerRun(record.REF)
        if record.INFO["PERIOD"] in [5,6] and hrun >= self.threshold+record.INFO["PERIOD"]:
            return hrun
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
                self.regions = None
                common.WARNING("%s not found"%filename) # TODO raise exception
            self.regions = BedTool(filename)
            if not self.regions._tabixed():
                sys.stderr.write("Creating tabix index for %s\n"%filename)
                self.regions.tabix(force=True)
        def __call__(self, record):
            interval = "%s:%s-%s"%(record.CHROM, record.POS, record.POS+len(record.REF))
            if self.regions is None: return None
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
# Call level filters - General
###################################

class Reason:
    name = ""
    def __init__(self):
        pass
    def GetReason(self):
        return self.name

class CallFilterMinValue(Reason):
    def __init__(self, name, field, threshold):
        self.name = name
        self.field = field
        self.threshold = threshold
    def __call__(self, sample):
        if sample[self.field] < self.threshold: return sample[self.field]
        else: return None

class CallFilterMaxValue(Reason):
    def __init__(self, name, field, threshold):
        self.name = name
        self.field = field
        self.threshold = threshold
    def __call__(self, sample):
        if sample[self.field] > self.threshold: return sample[self.field]
        else: return None

###################################
# Call level filters - HipSTR
###################################

class HipSTRCallFlankIndels(Reason):
    name = "HipSTRCallFlankIndels"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if 1.0*sample['DFLANKINDEL']/sample['DP'] > self.threshold:
            return 1.0*sample['DFLANKINDEL']/sample['DP']
        else: return None

class HipSTRCallStutter(Reason):
    name = "HipSTRCallStutter"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if 1.0*sample['DSTUTTER']/sample['DP'] > self.threshold:
            return 1.0*sample['DSTUTTER']/sample['DP']
        else: return None

class HipSTRCallMinSuppReads(Reason):
    name = "HipSTRMinSuppReads"
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

class GangSTRCallExpansionProbHom(Reason):
    name = "GangSTRCallExpansionProbHom"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if sample["QEXP"][2] < self.threshold: return sample["QEXP"][2]
        else: return None

class GangSTRCallExpansionProbHet(Reason):
    name = "GangSTRCallExpansionProbHet"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob het expansion
        if sample["QEXP"][1] < self.threshold: return sample["QEXP"][1]
        else: return None

class GangSTRCallExpansionProbTotal(Reason):
    name = "GangSTRCallExpansionProbTotal"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob het and hom expansion
        if (sample["QEXP"][1]+sample["QEXP"][2]) < self.threshold: return sample["QEXP"][1]+sample["QEXP"][2]
        else: return None

class GangSTRCallSpanOnly(Reason):
    name = "GangSTRCallSpanOnly"
    def __init__(self):
        pass
    def __call__(self, sample):
        #### Only spanning reads
        rcvals = [int(item) for item in sample["RC"].split(",")]
        if sample["DP"] == rcvals[1] : return rcvals[1]
        else: return None

class GangSTRCallSpanBoundOnly(Reason):
    name = "GangSTRCallSpanBoundOnly"
    def __init__(self):
        pass
    def __call__(self, sample):
        #### Only spanning and bounded
        rcvals = [int(item) for item in sample["RC"].split(",")]
        if sample["DP"] == rcvals[1]+rcvals[3] : return rcvals[1]+rcvals[3]
        else: return None

class GangSTRCallBadCI(Reason):
    name = "GangSTRCallBadCI"
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

class GangSTRCallRequireSupport(Reason):
    name = "GangSTRCallRequireSupport"
    def __init__(self, threshold, readlen):
        self.threshold = threshold
        self.readlen = readlen
    def __call__(self, sample):
        ### Require x number of reads supporting
        try:
            encl = dict([(int(item.split(",")[0]), int(item.split(",")[1])) for item in sample["ENCLREADS"].split("|")])
        except: encl = {}
        try:
            flank = dict([(int(item.split(",")[0]), int(item.split(",")[1])) for item in sample["FLNKREADS"].split("|")])
        except:
            flank = {}
        repcn = [int(item) for item in sample["REPCN"]]
        repcn_len = [len(item) for item in sample.gt_bases.split("/")]
        frrcount = int(sample["RC"].split(",")[2])
        for i in range(len(repcn)):
            allele = repcn[i]
            alen = repcn_len[i]
            # First, check enclosing. In all cases if we see good enclosings we're good to go
            if encl.get(allele, 0) >= self.threshold:
                continue
            else: # Fail if we *should* have seen enclosing
                if alen < 0.2*self.readlen: # TODO should this be stricter?
                    return self.threshold
            if alen > 2*self.readlen and frrcount < self.threshold:
                return self.threshold
            # For middle range, need to at least see some flanking
            numflank = sum(flank.values())
            if numflank < self.threshold:
                return self.threshold
        return None

###############################
# AdVNTR filters
###############################

###############################
# ExpansionHunter call-level filters
###############################
    
###############################
# PopSTR call level filters
###############################

class PopSTRCallRequireSupport(Reason):
    name = "PopSTRCallRequireSupport"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        read_support = sample["AD"]
        gts = sample.gt_alleles
        for gt in gts:
            if read_support[int(gt)] < self.threshold: return self.threshold
        return None
