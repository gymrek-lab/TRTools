"""
Utilities for harmonizing tandem repeat VCF records
across various formats output by TR genotyping tools
"""

class TRRecordHarmonizer:
    def __init__(self, vcffile, vcftype="auto"):
        self.vcftype = self.InferVCFType(vcffile, vcftype)

    def InferVCFType(self, vcffile, vcftype):
        return "gangstr" # TODO fill this in! and make these names not strings

    def HarmonizeRecord(self, vcfrecord):
        if self.vcftype == "gangstr":
            return TRRecord(vcfrecord, vcfrecord.REF, vcfrecord.ALT, vcfrecord.INFO["RU"])
        else:
            return None # TODO implement other types

# TODO add __iter__ to go through samples?
# So we could do for sample in trrecord
class TRRecord:
    def __init__(self, vcfrecord, ref_allele, alt_alleles, motif):
        self.vcfrecord = vcfrecord
        self.ref_allele = ref_allele
        self.alt_alleles = alt_alleles
        self.motif = motif

    def GetLengthGenotype(self, vcfsample):
        gts = vcfsample.gt_alleles
        gts_bases = [([self.ref_allele]+self.alt_alleles)[int(gt)] for gt in gts]
        return [float(len(item))/len(self.motif) for item in gts_bases]
