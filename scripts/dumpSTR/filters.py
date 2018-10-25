"""
Locus-level VCF filters
"""

import utils
import vcf.filters

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
        hwep = utils.GetSTRHWE(record, uselength=self.uselength)
        if hwep < self.threshold: return hwep
        else: return None

class Filter_MinLocusHet(vcf.filters.Base):
    'Filter sites with low heterozygosity'

    name = 'HETLOW'

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

class Filter_MaxLocusHet(vcf.filters.Base):
    'Filter sites with high heterozygosity'

    name = 'HETHIGH'

    def __init__(self, max_locus_het, uselength=False):
        self.threshold = max_locus_het
        self.uselength = uselength

    def __call__(self, record):
        if self.uselength:
            het = utils.GetLengthHet(record)
        else: het = record.heterozygosity
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
            self.filename = filename
        def __call__(self, record):
            return None # TODO
    f = Filter_Regions(name, filename)
    return Filter_Regions(name, filename)
