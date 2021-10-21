"""
Locus-level and Call-level VCF filters
"""

import ast
import os

import numpy as np
from pybedtools import BedTool

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

class FilterBase:
    '''
    Base class for locus level filters.
    Just defines the interface
    '''
    name = 'NotYetImplemented'

    def __call__(self, record: trh.TRRecord):
        raise NotImplementedError

    def filter_name(self):
        raise NotImplementedError

    def description(self):
        return ''

###################################
# Locus level filters
###################################

class Filter_MinLocusCallrate(FilterBase):
    """
    Class to filter VCF records by call rate

    This class extends Base

    Parameters
    ----------
    min_locus_callrate : float
       Filters calls with lower than this fraction called

    Attributes
    ----------
    threshold : float
       Filters calls with lower than this fraction called
       Derived from the input min_locus_callrate
    """

    name = 'CALLRATE'
    """The name of the filter"""

    def __init__(self, min_locus_callrate):
        self.threshold = min_locus_callrate

    def __call__(self, record: trh.TRRecord):
        if record.GetCallRate() < self.threshold: return record.GetCallRate()
        else: return None

    def filter_name(self):
        return self.name + str(self.threshold)

class Filter_MinLocusHWEP(FilterBase):
    """
    Class to filter VCF records by minimum HWE p-value

    This class extends Base

    Parameters
    ----------
    min_locus_hwep : float
       Filters calls with HWE p-value lower than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same

    Attributes
    ----------
    threshold : float
       Filters calls with HWE p-value lower than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same
    """

    name = 'HWE'
    """The name of the filter"""

    def __init__(self, min_locus_hwep, uselength=False):
        self.threshold = min_locus_hwep
        self.uselength = uselength

    def __call__(self, record: trh.TRRecord):
        allele_freqs = record.GetAlleleFreqs(uselength=self.uselength)
        genotype_counts = record.GetGenotypeCounts(uselength=self.uselength)
        hwep = utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)
        if hwep < self.threshold: return hwep
        else: return None

    def filter_name(self):
        return self.name + str(self.threshold)

class Filter_MinLocusHet(FilterBase):
    """
    Class to filter VCF records by minimum heterozygosity

    This class extends Base

    Parameters
    ----------
    min_locus_het : float
       Filters calls with heterozygosity lower than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same

    Attributes
    ----------
    threshold : float
       Filters calls with heterozygosity lower than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same
    """

    name = 'HETLOW'
    """The name of the filter"""

    def __init__(self, min_locus_het, uselength=False):
        self.threshold = min_locus_het
        self.uselength = uselength

    def __call__(self, record: trh.TRRecord):
        het = utils.GetHeterozygosity(record.GetAlleleFreqs(uselength=self.uselength))
        if het < self.threshold:
            return het
        return None

    def filter_name(self):
        return self.name + str(self.threshold)

class Filter_MaxLocusHet(FilterBase):
    """
    Class to filter VCF records by maximum heterozygosity

    This class extends Base

    Parameters
    ----------
    max_locus_het : float
       Filters calls with heterozygosity greater than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same

    Attributes
    ----------
    threshold : float
       Filters calls with heterozygosity greater than this
    vcftype: trh.VCFTYPES
        the type of the VCF we're working with
    uselength : bool, optional
       If set to true, consider all alleles with the same length as the same
    """

    name = 'HETHIGH'
    """The name of the filter"""

    def __init__(self, max_locus_het, uselength=False):
        self.threshold = max_locus_het
        self.uselength = uselength

    def __call__(self, record: trh.TRRecord):
        het = utils.GetHeterozygosity(record.GetAlleleFreqs(uselength=self.uselength))
        if het > self.threshold:
            return het
        return None

    def filter_name(self):
        return self.name + str(self.threshold)

class Filter_LocusHrun(FilterBase):
    """
    Class to filter VCF records for penta- or hexanucleotide STRs with long homopolymer runs

    This only works on HipSTR VCFs. STRs with long homopolymer runs have been
    shown to be difficult for HipSTR to call.
    This filter removes 5-mers with homopolymer runs >= len 5
    and 6-mers with homopolymer runs >= len 6
    """

    name = 'HRUN'
    """The name of the filter"""

    def __init__(self):
        pass

    def __call__(self, record: trh.TRRecord):
        if record.HasFullStringGenotypes():
            hrun = utils.GetHomopolymerRun(record.full_alleles[0])
        else:
            hrun = utils.GetHomopolymerRun(record.ref_allele)
        if "PERIOD" not in record.info: return None # Don't apply if we don't know the period
        if record.info["PERIOD"] in [5, 6] and hrun >= record.info["PERIOD"]:
            return hrun
        return None

    def filter_name(self):
        return self.name

def create_region_filter(name, filename):
    """Creates a locus-level filter based on a file of regions.

    Builds and returns a class extending Base that
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
    filter_regions : Base object.
       Returns None if we fail to load the regions
    -------

    """
    class Filter_Regions(FilterBase):
        def __init__(self, name, filename):
            self.threshold = ""
            self.name = name
            self.pass_checks = True
            self.LoadRegions(filename)
        def LoadRegions(self, filename):
            if not filename.endswith(".bed.gz") and not filename.endswith(".bed.bgz"):
                #raise ValueError("Make sure %s is bgzipped and indexed"%filename)
                self.regions = None
                common.WARNING("Make sure %s is bgzipped and indexed"%filename)
                self.pass_checks = False
                return
            if not os.path.isfile(filename):
                #raise ValueError("Could not find regions BED file %s"%filename)
                self.regions = None
                common.WARNING("Could not find regions BED file %s"%filename)
                self.pass_checks = False
                return
            if not os.path.isfile(filename+".tbi"):
                #raise ValueError("Could not find tabix index %s.tbi"%filename)
                self.regions = None
                common.WARNING("Could not find tabix index %s.tbi"%filename)
                self.pass_checks = False
                return
            self.regions = BedTool(filename)
        def __call__(self, record: trh.TRRecord):
            interval = "%s:%s-%s"%(record.chrom, record.pos,
                                   record.pos + record.ref_allele_length)
            if self.regions is None: return None
            if "chr" in interval:
                interval2 = interval.replace("chr","")
            else: interval2 = "chr"+interval
            # Try with and without chr
            tb1 = self.regions.tabix_intervals(interval)
            if tb1.count() > 0: return self.name
            tb2 = self.regions.tabix_intervals(interval2)
            if tb2.count() > 0: return self.name
            return None
        def filter_name(self):
            return self.name
        def description(self):
            return 'Filter TRs overlapping this region'
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
    gets applied to each call. The __call__ function returns a 1D array of
    values, one per sample. For numeric arrays, nan values indicate the sample
    wasn't filtered, and any other value indicates that the sample was
    filtered (where the value indicates why).
    """

    name = ""
    """
    The name of the filter to put in the FORMAT:FILTER field of
    filtered calls.
    """

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
        self.name = name + str(threshold)
        self.field = field
        self.threshold = threshold
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        fieldvals = record.format[self.field][:, 0]
        sample_filter[fieldvals < self.threshold] = fieldvals[fieldvals < self.threshold]
        return sample_filter

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
        self.name = name + str(threshold)
        self.field = field
        self.threshold = threshold
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        fieldvals = record.format[self.field][:, 0]
        sample_filter[fieldvals > self.threshold] = fieldvals[fieldvals > self.threshold]
        return sample_filter

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
    """

    name = "HipSTRCallFlankIndels"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        ratio = record.format['DFLANKINDEL'][:, 0]/record.format['DP'][:, 0]
        sample_filter[ratio <= self.threshold] = np.nan
        sample_filter[ratio > self.threshold] = ratio[ratio > self.threshold]
        return sample_filter

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
    """

    name = "HipSTRCallStutter"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        ratio = record.format['DSTUTTER'][:, 0]/record.format['DP'][:, 0]
        sample_filter[ratio <= self.threshold] = np.nan
        sample_filter[ratio > self.threshold] = ratio[ratio > self.threshold]
        return sample_filter

class HipSTRCallMinSuppReads(Reason):
    """Filter HipSTR calls for which alleles are supported by too few reads

    Extends Reason class.
    Filters on the number of reads supporting each called allele.
    Based on FORMAT:ALLREADS and FORMAT:GB fields.
    Assumes that number of supporting reads is zero if:
    * that read length is not present in ALLREADS
    * or ALLREADS is unset for that sample ('.')
    * or ALLREADS is not present at that locus

    Parameters
    ----------
    threshold : int
        Minimum number of reads supporting each allele

    Attributes
    ----------
    threshold : int
        Minimum number of reads supporting each allele
    """

    name = "HipSTRMinSuppReads"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return np.full((record.GetNumSamples()), np.nan)

        if "ALLREADS" not in record.format:
            return np.zeros((record.GetNumSamples()), dtype=float)
        samples_to_check = (called_samples &
                          (record.format["ALLREADS"] != '') &
                          (record.format["ALLREADS"] != '.'))

        if not np.any(samples_to_check):
            # all samples were either not called or were missing ALLREADS
            # say that we filtered all the samples which were missing ALLREADS
            sample_filter = np.full((record.GetNumSamples()), np.nan)
            sample_filter[called_samples] = 0
            return sample_filter

        # Going to assume that either all samples are phased or none are
        first_gb = record.format["GB"][samples_to_check][0]
        if "/" in first_gb:
            delim = "/"
        elif "|" in first_gb:
            delim = '|'
        else:
            raise ValueError("Cant't identify phasing char ('|' or '/') in GB field")

        gb = np.char.split(record.format["GB"][samples_to_check], delim)
        gb = np.stack(gb).astype(int)
        # Format allreads like a python dict literal so we can interpret it
        # like that
        allreads = np.char.replace(record.format["ALLREADS"][samples_to_check], ";", ',')
        allreads = np.char.replace(allreads, '|', ':')
        allreads = np.char.add('{', np.char.add(allreads, '}'))
        min_counts = np.full((record.GetNumSamples()), np.nan)
        for idx, single_allreads in enumerate(allreads):
            reads_dict = ast.literal_eval(single_allreads)
            min_count = np.inf
            for gt in gb[idx, :]:
                gt = int(gt)
                if gt not in reads_dict:
                    min_count = 0
                else:
                    min_count = min(min_count, reads_dict[gt])
            min_counts[np.nonzero(samples_to_check)[0][idx]] = min_count
        min_counts[min_counts >= self.threshold] = np.nan
        # if ALLREADS is missing but the sample is present, filter
        min_counts[called_samples & ~samples_to_check] = 0
        return min_counts

###############################
# GangSTR filters
###############################

class GangSTRCallExpansionProbHom(Reason):
    """Filter GangSTR calls with low probability for homozygous expansion.

    Extends Reason class.
    Based on the QEXP field. Filter if QEXP[:, 2] (prob hom expansion above INFO:THRESHOLD)
    is less than the threshold.

    Parameters
    ----------
    threshold : float
        Minimum homozygous expansion probability

    Attributes
    ----------
    threshold : float
        Minimum homozygous expansion probability
    """

    name = "GangSTRCallExpansionProbHom"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        hom_expansion_prob = record.format["QEXP"][called_samples, 2]
        sample_filter[np.nonzero(called_samples)[0][hom_expansion_prob < self.threshold]] = \
            hom_expansion_prob[hom_expansion_prob < self.threshold]
        return sample_filter

class GangSTRCallExpansionProbHet(Reason):
    """Filter GangSTR calls with low probability for heterozygous expansion.

    Extends Reason class.
    Based on the QEXP field. Filter if QEXP[:, 1] (prob het expansion above INFO:THRESHOLD)
    is less than the threshold.

    Parameters
    ----------
    threshold : float
        Minimum heterozygous expansion probability

    Attributes
    ----------
    threshold : float
        Minimum heterozygous expansion probability
    """

    name = "GangSTRCallExpansionProbHet"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        het_expansion_prob = record.format["QEXP"][called_samples, 1]
        sample_filter[np.nonzero(called_samples)[0][het_expansion_prob < self.threshold]] = \
            het_expansion_prob[het_expansion_prob < self.threshold]
        return sample_filter

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
    """

    name = "GangSTRCallExpansionProbTotal"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        expansion_prob = record.format["QEXP"][called_samples, 1] + \
                record.format["QEXP"][called_samples, 2]
        sample_filter[np.nonzero(called_samples)[0][expansion_prob < self.threshold]] = \
                expansion_prob[expansion_prob < self.threshold]
        return sample_filter

class GangSTRCallSpanOnly(Reason):
    """Filter GangSTR calls where only spanning reads were identified.

    Extends Reason class. Based on RC field.

    """

    name = "GangSTRCallSpanOnly"
    def __init__(self):
        pass
    def __call__(self, record: trh.TRRecord):
        #### Only spanning reads
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        rcvals = np.char.split(record.format['RC'][called_samples], ',')
        rcvals = np.stack(rcvals, axis=0).astype(int)
        filter_indicies = rcvals[:, 1] == record.format['DP'][called_samples, 0]
        sample_filter[np.nonzero(called_samples)[0][filter_indicies]] = \
            rcvals[:, 1][filter_indicies]
        return sample_filter

class GangSTRCallSpanBoundOnly(Reason):
    """Filter GangSTR calls where only spanning or flanking reads were identified.

    Extends Reason class. Based on RC field.

    """

    name = "GangSTRCallSpanBoundOnly"
    """The name of the filter"""

    def __init__(self):
        pass
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        rcvals = np.char.split(record.format['RC'][called_samples], ',')
        rcvals = np.stack(rcvals, axis=0).astype(int)
        span_bound = rcvals[:, 1] + rcvals[:, 3]
        filter_indicies = span_bound == record.format['DP'][called_samples, 0]
        sample_filter[np.nonzero(called_samples)[0][filter_indicies]] = \
            span_bound[filter_indicies]
        return sample_filter

class GangSTRCallBadCI(Reason):
    """Filter GangSTR calls where the ML genotype estimate is outside of CI.

    If 95% confidence interval does not include the maximum likelihood
    genotype call, the call is filtered.

    Extends Reason class. Based on REPCI and REPCN fields.

    """

    name = "GangSTRCallBadCI"
    """The name of the filter"""

    def __init__(self):
        pass
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        called_samples = record.GetCalledSamples()
        if not np.any(called_samples):
            return sample_filter
        ml = record.format["REPCN"][called_samples]
        ci = np.char.split(record.format["REPCI"][called_samples], ",")
        ci = np.stack(ci)
        ci = np.char.split(ci, '-')
        ci = np.array(ci.tolist(), dtype=int) # now sample x ploidy x 2 (min, max)
        filter_per_gt = np.logical_or(ml < ci[:, :, 0], ci[:, :, 1] < ml)
        filter_indicies = np.any(filter_per_gt, axis=1)
        if not np.any(filter_indicies):
            return sample_filter
        problem_gt_indicies = np.argmax(filter_per_gt[filter_indicies, :],
                                        axis=1)
        sample_filter[np.nonzero(called_samples)[0][filter_indicies]] = \
                ml[filter_indicies, problem_gt_indicies]
        return sample_filter

'''
Since this class is experimental, I'm just not converting it from pyvcf to
cyvcf2 for now
class GangSTRCallRequireSupport(Reason)
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
    """

    name = "GangSTRCallRequireSupport"
    """The name of the filter"""

    def __init__(self, threshold, readlen):
        self.threshold = threshold
        self.readlen = readlen
    def __call__(self, record: trh.TRRecord):
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
'''

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
    """

    name = "PopSTRCallRequireSupport"
    """The name of the filter"""

    def __init__(self, threshold):
        self.threshold = threshold
        self.name += str(threshold)
    def __call__(self, record: trh.TRRecord):
        sample_filter = np.full((record.GetNumSamples()), np.nan)
        sample_list = np.arange(record.GetNumSamples())
        read_support = record.format["AD"] # samples x n_alleles
        gt_indicies = record.GetGenotypeIndicies()[:, :-1]
        for ploid in range(gt_indicies.shape[1]):
            # iterate over ploidy
            new_filters = read_support[sample_list, gt_indicies[:, ploid]] < self.threshold
            sample_filter[new_filters] = read_support[new_filters, gt_indicies[:, ploid]]
        return sample_filter

