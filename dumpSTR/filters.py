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
else: # pragma: no cover
    import strtools.utils.common as common # pragma: no cover
    import strtools.utils.tr_harmonizer as trh # pragma: no cover
    import strtools.utils.utils as utils # pragma: no cover

###################################
# Locus level filters
###################################

class Filter_MinLocusCallrate(vcf.filters.Base):
    """
    Class to filter VCF records by call rate

    This class extends vcf.filters.Base

    Parameters
    ----------
    min_locus_callrate : float
       Filters calls with lower than this fraction called

    Attributes
    ----------
    threshold : float
       Call rate threshold
    name : str
       Set to the name of the filter ('CALLRATE')
    """

    name = 'CALLRATE'
    def __init__(self, min_locus_callrate):
        self.threshold = min_locus_callrate

    def __call__(self, record):
        if record.call_rate < self.threshold: return record.call_rate
        else: return None

class Filter_MinLocusHWEP(vcf.filters.Base):
    """
    Class to filter VCF records by minimum HWE p-value

    This class extends vcf.filters.Base

    Parameters
    ----------
    min_locus_hwep : float
       Filters calls with HWE p-value lower than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    

    Attributes
    ----------
    threshold : float
       Filters calls with HWE p-value lower than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    
    name : str
       Set to the name of the filter ('HWE')
    """

    name = 'HWE'
    def __init__(self, min_locus_hwep, tr_harmonizer, uselength=False):
        self.threshold = min_locus_hwep
        self.tr_harmonizer = tr_harmonizer
        self.uselength = uselength

    def __call__(self, record):
        trrecord = self.tr_harmonizer.HarmonizeRecord(record)
        allele_freqs = trrecord.GetAlleleFreqs(uselength=self.uselength)
        genotype_counts = trrecord.GetGenotypeCounts(uselength=self.uselength)
        hwep = utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)
        if hwep < self.threshold: return hwep
        else: return None

class Filter_MinLocusHet(vcf.filters.Base):
    """
    Class to filter VCF records by minimum heterozygosity

    This class extends vcf.filters.Base

    Parameters
    ----------
    min_locus_het : float
       Filters calls with heterozygosity lower than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    

    Attributes
    ----------
    threshold : float
       Filters calls with heterozygosity lower than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    
    name : str
       Set to the name of the filter ('HETLOW')
    """

    name = 'HETLOW'
    def __init__(self, min_locus_het, tr_harmonizer, uselength=False):
        self.threshold = min_locus_het
        self.tr_harmonizer = tr_harmonizer
        self.uselength = uselength

    def __call__(self, record):
        trrecord = self.tr_harmonizer.HarmonizeRecord(record)
        het = utils.GetHeterozygosity(trrecord.GetAlleleFreqs(uselength=self.uselength))
        if het < self.threshold:
            return het
        return None

class Filter_MaxLocusHet(vcf.filters.Base):
    """
    Class to filter VCF records by maximum heterozygosity

    This class extends vcf.filters.Base

    Parameters
    ----------
    max_locus_het : float
       Filters calls with heterozygosity greater than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    

    Attributes
    ----------
    threshold : float
       Filters calls with heterozygosity greater than this
    tr_harmonizer : trh.TRRecordHarmonizer object
       Used to harmonize VCF records across genotyping tools
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same    
    name : str
       Set to the name of the filter ('HETHIGHT')
    """

    name = 'HETHIGH'
    def __init__(self, max_locus_het, tr_harmonizer, uselength=False):
        self.threshold = max_locus_het
        self.tr_harmonizer = tr_harmonizer
        self.uselength = uselength

    def __call__(self, record):
        trrecord = self.tr_harmonizer.HarmonizeRecord(record)
        het = utils.GetHeterozygosity(trrecord.GetAlleleFreqs(uselength=self.uselength))
        if het > self.threshold:
            return het
        return None

class Filter_LocusHrun(vcf.filters.Base):
    """
    Class to filter VCF records for penta- or hexanucleotide STRs with long homopolymer runs

    This is relevant only for HipSTR VCFs. 5-mer and 6-mer STRs with long T runs have been
    shown to be problematic. This filter removes 5-mers with stretches of T >= 5bp or 6-mers
    with stretches of T >= 6bp. If PERIOD is not in record.INFO, this class does nothing.

    This class extends vcf.filters.Base

    Attributes
    ----------
    name : str
       Set to the name of the filter ('HRUN')
    """
    
    name = 'HRUN'
    def __init__(self):
        self.threshold = 0 # e.g. pentamers with hruns of 5+threshold, hexa with 6+threshold

    def __call__(self, record):
        hrun = utils.GetHomopolymerRun(record.REF)
        if "PERIOD" not in record.INFO: return None # Don't apply if we don't know the period
        if record.INFO["PERIOD"] in [5,6] and hrun >= self.threshold+record.INFO["PERIOD"]:
            return hrun
        return None

def create_region_filter(name, filename):
    """Creates a locus-level filter based on a file of regions.

    Builds and returns a class extending vcf.filters.Base that
    can be used to filter any records overlapping intervals in the
    input BED file

    Parameters
    ----------
    name : str
       Name of the region filter to create.
       This will go in the FILTER field of the output VCF
    filename : str
       BED file containing the regions. Must be sorted by chrom, start
       If it's not bgzipped and tabixed, we'll attempt to do that.

    Returns
    filter_regions : vcf.filters.Base object.
       Returns None if we fail to load the regions
    -------

    """
    class Filter_Regions(vcf.filters.Base):
        'Filter regions from file'
        def __init__(self, name, filename):
            self.threshold = ""
            self.name = name
            self.pass_checks = True
            self.LoadRegions(filename)
        def LoadRegions(self, filename):
            if not os.path.exists(filename):
                self.regions = None
                common.WARNING("%s not found"%filename)
                self.pass_checks = False
                return
            self.regions = BedTool(filename)
            if not self.regions._tabixed():
                sys.stderr.write("Creating tabix index for %s\n"%filename)
                try:
                    self.regions.tabix(force=True)
                except:
                    self.pass_checks = False
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
    if not f.pass_checks: return None
    return f

###################################
# Call level filters - General
###################################

class Reason:
    """Base call-level filter class.

    Other classes extend this for each different call-level filter.
    Classes that extend this must implement a __call__ function that
    gets applied to each call. The __call__ function returns None for
    calls passing the filter, and something besides None for samples
    failing the filter.

    Attributes
    ----------
    name : str
        The name of the filter to put in the FORMAT:FILTER field of 
        filtered calls.
    """
    
    name = ""
    def GetReason(self):
        return self.name

class CallFilterMinValue(Reason):
    """Generic call-level filter based on minimum allowed value for a field.

    Extends Reason class.
    For any call-level value, such as DP, this class can be used to make
    a filter based on the minimum allowed value for that field.

    Parameters
    ----------
    name : str
        The name of the filter to put in the FORMAT:FILTER field of 
        filtered calls.
    field : str
        The FORMAT field to filter on
    threshold : float
        The minimum allowed value for the field.

    Attributes
    ----------
    name : str
        The name of the filter to put in the FORMAT:FILTER field of 
        filtered calls.
    field : str
        The FORMAT field to filter on
    threshold : float
        The minimum allowed value for the field.    

    Examples
    --------
    >>> min_dp_filt = CallFilterMinValue("LOWDP","DP",10)
    """
    
    def __init__(self, name, field, threshold):
        self.name = name
        self.field = field
        self.threshold = threshold
    def __call__(self, sample):
        if sample[self.field] < self.threshold: return sample[self.field]
        else: return None

class CallFilterMaxValue(Reason):
    """Generic call-level filter based on maximum allowed value for a field.

    Extends Reason class.
    For any call-level value, such as DP, this class can be used to make
    a filter based on the maximum allowed value for that field.

    Parameters
    ----------
    name : str
        The name of the filter to put in the FORMAT:FILTER field of 
        filtered calls.
    field : str
        The FORMAT field to filter on
    threshold : float
        The maximum allowed value for the field.

    Attributes
    ----------
    name : str
        The name of the filter to put in the FORMAT:FILTER field of 
        filtered calls.
    field : str
        The FORMAT field to filter on
    threshold : float
        The maximum allowed value for the field.    

    Examples
    --------
    >>> max_dp_filt = CallFilterMaxValue("HIGHDP","DP",1000)
    """

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
    """Filter HipSTR calls with many indels in flanks

    Extends Reason class.
    Filters on the percentage of reads with indels in flanks.
    Based on FORMAT:DFLANKINDEL and FORMAT:DP fields.

    Parameters
    ----------
    threshold : float
        Minimum percent of reads that can have indels in their flanks
    
    Attributes
    ----------
    threshold : float
        Minimum percent of reads that can have indels in their flanks
    name : str
        Set to "HipSTRCallFlankIndels"
    """
    
    name = "HipSTRCallFlankIndels"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if 1.0*sample['DFLANKINDEL']/sample['DP'] > self.threshold:
            return 1.0*sample['DFLANKINDEL']/sample['DP']
        else: return None

class HipSTRCallStutter(Reason):
    """Filter HipSTR calls with many stutter reads

    Extends Reason class.
    Filters on the percentage of reads with stutter errors
    Based on FORMAT:DSTUTTER and FORMAT:DP fields.

    Parameters
    ----------
    threshold : float
        Minimum percent of reads that can have stutter errors
    
    Attributes
    ----------
    threshold : float
        Minimum percent of reads that can have stutter errors
    name : str
        Set to "HipSTRCallStutter"
    """
    
    name = "HipSTRCallStutter"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if 1.0*sample['DSTUTTER']/sample['DP'] > self.threshold:
            return 1.0*sample['DSTUTTER']/sample['DP']
        else: return None

class HipSTRCallMinSuppReads(Reason):
    """Filter HipSTR calls for which alleles are supported by too few reads

    Extends Reason class.
    Filters on the number of reads supporting each called allele.
    Based on FORMAT:ALLREADS and FORMAT:GB fields.

    Parameters
    ----------
    threshold : int
        Minimum number of reads supporting each allele
    
    Attributes
    ----------
    threshold : int
        Minimum number of reads supporting each allele
    name : str
        Set to "HipSTRMinSuppReads"
    """

    name = "HipSTRMinSuppReads"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if sample["ALLREADS"] is None: return 0
        delim = "|"
        if "/" in sample["GB"]: delim = "/"
        gb = [int(item) for item in sample["GB"].split(delim)]
        allreads = sample["ALLREADS"].split(";")
        r1 = 0
        r2 = 0
        for item in allreads:
            allele, readcount = [int(item) for item in item.split("|")]
            if allele == gb[0]: r1 += readcount
            if allele == gb[1]: r2 += readcount
        min_read_count = min([r1, r2])
        if min_read_count < self.threshold: return min_read_count
        else: return None

###############################
# GangSTR filters
###############################

class GangSTRCallExpansionProbHom(Reason):
    """Filter GangSTR calls with low probability for homozygous expansion.

    Extends Reason class.
    Based on the QEXP field. Filter if QEXP[2] (prob hom expansion above INFO:THRESHOLD)
    is less than the threshold.

    Parameters
    ----------
    threshold : float
        Minimum homozygous expansion probability
    
    Attributes
    ----------
    threshold : float
        Minimum homozygous expansion probability
    name : str
        Set to "GangSTRCallExpansionProbHom"
    """

    name = "GangSTRCallExpansionProbHom"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        if sample["QEXP"][2] < self.threshold: return sample["QEXP"][2]
        else: return None

class GangSTRCallExpansionProbHet(Reason):
    """Filter GangSTR calls with low probability for heterozygous expansion.

    Extends Reason class.
    Based on the QEXP field. Filter if QEXP[1] (prob het expansion above INFO:THRESHOLD)
    is less than the threshold.

    Parameters
    ----------
    threshold : float
        Minimum heterozygous expansion probability
    
    Attributes
    ----------
    threshold : float
        Minimum heterozygous expansion probability
    name : str
        Set to "GangSTRCallExpansionProbHet"
    """

    name = "GangSTRCallExpansionProbHet"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob het expansion
        if sample["QEXP"][1] < self.threshold: return sample["QEXP"][1]
        else: return None

class GangSTRCallExpansionProbTotal(Reason):
    """Filter GangSTR calls with low probability for expansion (heterozygous or homozygous)

    Extends Reason class.
    Based on the QEXP field. Filter if QEXP[1]+QEXP[2] (prob het or hom expansion above INFO:THRESHOLD)
    is less than the threshold.

    Parameters
    ----------
    threshold : float
        Minimum expansion probability
    
    Attributes
    ----------
    threshold : float
        Minimum expansion probability
    name : str
        Set to "GangSTRCallExpansionProbTotal"
    """

    name = "GangSTRCallExpansionProbTotal"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        #### Prob het and hom expansion
        if (sample["QEXP"][1]+sample["QEXP"][2]) < self.threshold: return sample["QEXP"][1]+sample["QEXP"][2]
        else: return None

class GangSTRCallSpanOnly(Reason):
    """Filter GangSTR calls where only spanning reads were identified.

    Extends Reason class. Based on RC field.

    Attributes
    ----------
    name : str
        Set to "GangSTRCallSpanOnly"
    """
    
    name = "GangSTRCallSpanOnly"
    def __init__(self):
        pass
    def __call__(self, sample):
        #### Only spanning reads
        rcvals = [int(item) for item in sample["RC"].split(",")]
        if sample["DP"] == rcvals[1] : return rcvals[1]
        else: return None

class GangSTRCallSpanBoundOnly(Reason):
    """Filter GangSTR calls where only spanning or flanking reads were identified.

    Extends Reason class. Based on RC field.

    Attributes
    ----------
    name : str
        Set to "GangSTRCallSpanBoundOnly"
    """

    name = "GangSTRCallSpanBoundOnly"
    def __init__(self):
        pass
    def __call__(self, sample):
        #### Only spanning and bounded
        rcvals = [int(item) for item in sample["RC"].split(",")]
        if sample["DP"] == rcvals[1]+rcvals[3] : return rcvals[1]+rcvals[3]
        else: return None

class GangSTRCallBadCI(Reason):
    """Filter GangSTR calls where the ML genotype estimate is outside of CI.

    If 95% confidence interval does not include the maximum likelihood
    genotype call, the call is filtered.

    Extends Reason class. Based on REPCI and REPCN fields.

    Attributes
    ----------
    name ; str
        Set to "GangSTRCallBadCI"
    """
    
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
    """Fitler GangSTR calls where alleles not supported by enough reads.

    Extends Reason class. 
    This filter is still experimental and should be used with caution.
    For short alleles, require enclosing reads. If we see enough enclosing, the
    call is kept. If the allele was very short (20% of readlength) and we don't see
    enough enclosing, filter. If the call is longer, require that we at least see 
    supporting flanking reads.

    Parameters
    ----------
    threshold : int
       Minimum required number of supporting reads
    readlen : int
       Read length used for GangSTR calling.

    Attributes
    ----------
    threshold : int
       Minimum required number of supporting reads
    readlen : int
       Read length used for GangSTR calling. 
    name : str
       Set to "GangSTRCallRequireSupport"
    """
    
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
    """Filter PopSTR calls not supported by enough reads

    Extends Reason class.
    Relies on FORMAT:AD field.

    Parameters
    ----------
    threshold : int
       Require this many reads supporting each called allele.

    Attributes
    ----------
    threshold : int
       Require this many reads supporting each called allele.
    name : str
       Set to "PopSTRCallRequireSupport"
    """
    
    name = "PopSTRCallRequireSupport"
    def __init__(self, threshold):
        self.threshold = threshold
    def __call__(self, sample):
        read_support = sample["AD"]
        gts = sample.gt_alleles
        for gt in gts:
            if read_support[int(gt)] < self.threshold: return self.threshold
        return None

