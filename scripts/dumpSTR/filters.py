"""
Locus-level VCF filters
"""

import vcf.filters

class Filter_MinLocusCallrate(vcf.filters.Base):
    'Filter site by call rate'

    name = 'CALLRATE'

    def __init__(self, min_locus_callrate):
        self.threshold = min_locus_callrate

    def __call__(self, record):
        return None # TODO

class Filter_MinLocusHWEP(vcf.filters.Base):
    'Filter site by HWE'

    name = 'HWE'

    def __init__(self, min_locus_hwep, uselength=False):
        self.threshold = min_locus_hwep
        self.uselength = uselength

    def __call__(self, record):
        return None # TODO

class Filter_MinLocusHet(vcf.filters.Base):
    'Filter sites with low heterozygosity'

    name = 'HETLOW'

    def __init__(self, min_locus_het, uselength=False):
        self.threshold = min_locus_het
        self.uselength = uselength

    def __call__(self, record):
        return None # TODO

class Filter_MaxLocusHet(vcf.filters.Base):
    'Filter sites with high heterozygosity'

    name = 'HETHIGH'

    def __init__(self, max_locus_het, uselength=False):
        self.threshold = max_locus_het
        self.uselength = uselength

    def __call__(self, record):
        return None # TODO

class Filter_LocusHrun(vcf.filters.Base):
    'Filter penta/hexa sites with long homopolymer runs'
    
    name = 'HRUN'

    def __init__(self):
        self.threshold = 1 # e.g. pentamers with hruns of 5+1, hexa with 6+1

    def __call__(self, record):
        return None

def create_region_filter(name):
    return None # TODO for list of region files/names
