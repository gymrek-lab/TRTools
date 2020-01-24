"""
Utilities for harmonizing tandem repeat VCF records
across various formats output by TR genotyping tools
"""

import enum
import math

types =  enum.Enum('Types', ['gangstr', 'advntr', 'hipstr', 'eh', 'popstr']) #TODO add Beagle

class TRRecordHarmonizer:
    #vcffile - a pyvcf reader instance
    #vcftype - the type of the program that produced this vcf file
    #it defaults to auto and is inferred from the vcf metadata
    #however, if the name of the program is not obvious from the command metadata
    #line in the vcf, then it will need to be manually specified
    #see the above list for possibilities
    #peculiarities of particular types:
    #gangstr
    #  does not include impurities, TODO partial repeats?
    #hipstr
    #  both ref and alt alleles may contain impurities (including partial repeats)
    #  hipstr sometimes includes non-repeat flanks in the original genotypes,
    #  which are trimmed off during harmonization. If the alt alleles
    #  have differently sized flanks than the ref allele then those
    #  alt alleles will be improperly trimmed
    #  hipstr reports the motif length but no the motif itself
    #  so the motif is inferred from the sequence. This normally results in an accurate
    #  motif, but in theory could be wrong on occasion.
    #advntr
    #  uses its own peculiar ID convention that may change from run to run
    #  may contain impurities or partial repeats?
    #popstr
    #  includes written out reference alleles (with impurities)
    #  but only repeat counts for alternate alleles. So the alternate alleles we output
    #  start at the beginning of the motif (even if its in the wrong part of the cycle),
    #  never contain impurities, and in the case of non-integer repeat counts, we
    #  include extra characters from the motif until we would exceed the count.
    #  since repeat counts are rounded to the nearest tenth, this means we could
    #  miss/add extra base pairs that weren't there
    def __init__(self, vcffile, vcftype="auto"):
        self.vcffile = vcffile
        if vcftype == "auto":
            self.vcftype = self.InferVCFType(vcffile)
        else:
            if vcftype not in types.__members__:
                raise ValueError("{} is not an excepted TR vcf type. Expected one of {}".format(vcftype, list(types.__members__)))
            self.vcftype = types[vcftype]

    def InferVCFType(self, vcffile):
        possible_vcf_types = set()
        for key, value in vcffile.metadata.items():
            #sometimes value is a list of values, other times it is the value itself
            #no documentation as to why this is the case, so consider both as possible
            if type(value) == list:
                values = value
            else:
                values = {value}
            for value in values:
                if key.upper() == 'COMMAND' and 'GANGSTR' in value.upper():
                    possible_vcf_types.add(types.gangstr)
                if key.upper() == 'COMMAND' and 'HIPSTR' in value.upper():
                    possible_vcf_types.add(types.hipstr)
                if key.upper() == 'SOURCE' and 'ADVNTR' in value.upper():
                    possible_vcf_types.add(types.advntr)
                if key.upper() == 'SOURCE' and 'POPSTR' in value.upper():
                    possible_vcf_types.add(types.popstr)
                #TODO expansion hunter - it doesn't document itself well

        if len(possible_vcf_types) == 0:
            raise ValueError('Could not identify the type of this vcf')

        if len(possible_vcf_types) > 1:
            raise ValueError(f'Confused - this vcf looks like it could have been any of the types: {possible_vcf_types}')

        return next(iter(possible_vcf_types))

    def __iter__(self):
        for record in self.vcffile:
            yield self.HarmonizeRecord(record)

    #harmonize to the gangstr standard of 
    #catcatcat catcatcatcat,catcat
    def HarmonizeRecord(self, vcfrecord):
        ref_allele = None
        alt_alleles = None
        motif = None

        if self.vcftype == types.gangstr:
            ref_allele = vcfrecord.REF.upper()
            if vcfrecord.ALT[0] is not None:
                    alt_alleles = []
                    for alt in vcfrecord.ALT:
                        alt_alleles.append(str(alt).upper())
            else:
                alt_alleles = [None]
            motif = vcfrecord.INFO["RU"].upper()
            id = vcfrecord.ID

        elif self.vcftype == types.hipstr:
            #remove the flanking variants
            #note that if there are both flanking variants
            #and some of the alterante alleles have insertions/deletions in the flanks
            #then this might inadvertantly elongate/truncate those alternate alleles
            pos = int(vcfrecord.POS)
            start_offset = int(vcfrecord.INFO['START']) - pos
            pos_end_offset = int(vcfrecord.INFO['END']) - pos
            neg_end_offset = pos_end_offset + 1 - len(vcfrecord.REF)

            # neg_end_offset is the number of flanking non repeat bp to remove from the
            # end of each allele
            # e.g. 'AAAT'[0:-1] == 'AAA'
            # however, if neg_end_offset == 0, then we would get
            # 'AAAA'[0:0] == '' which is not the intent
            #so we need an if statement to instead write 'AAAA'[0:]
            #which gives us 'AAAA'
            if neg_end_offset == 0:
                ref_allele = vcfrecord.REF[start_offset:].upper()
                if vcfrecord.ALT[0] is not None:
                    alt_alleles = []
                    for alt in vcfrecord.ALT:
                        alt_alleles.append(str(alt)[start_offset:].upper())
                else:
                    alt_alleles = [None]
            else:
                ref_allele = vcfrecord.REF[start_offset:neg_end_offset].upper()
                if vcfrecord.ALT[0] is not None:
                    alt_alleles = []
                    for alt in vcfrecord.ALT:
                        alt_alleles.append(str(alt)[start_offset:neg_end_offset].upper())
                else:
                    alt_alleles = [None]

            #get the motif. Hipstr doesn't tell us this explicitly, so figure it out
            #from the kmers in the sequence
            #do this by looking for the kmer of the given length which appears
            #most frequently at the same offset
            motif_len = int(vcfrecord.INFO["PERIOD"])
            best_kmer = None
            best_copies = 0
            for offset in range(0, motif_len):
                kmers = {}
                start_idx = start_offset
                while start_idx + motif_len <= len(ref_allele):
                    kmer = ref_allele[start_idx:(start_idx + motif_len)]
                    if kmer not in kmers:
                        kmers[kmer] = 1
                    else:
                        kmers[kmer] += 1
                    start_idx += motif_len
                current_best_kmer = max(kmers, key = lambda k: kmers[k])
                current_best_copies = kmers[current_best_kmer]
                if current_best_copies > best_copies:
                    best_kmer = current_best_kmer
                    best_copies = current_best_copies
                        
            motif = best_kmer
            id = vcfrecord.ID

        elif self.vcftype == types.advntr:
            ref_allele = vcfrecord.REF.upper()
            if vcfrecord.ALT[0] is not None:
                    alt_alleles = []
                    for alt in vcfrecord.ALT:
                        alt_alleles.append(str(alt).upper())
            else:
                alt_alleles = [None]

            motif = vcfrecord.INFO["RU"].upper()
            id = vcfrecord.INFO["VID"]
        
        elif self.vcftype == types.eh:
            raise NotImplementedError("Haven't implemented expansion hunter harmonization yet")
        elif self.vcftype == types.popstr:
            ref_allele = vcfrecord.REF.upper()
            motif = vcfrecord.INFO["Motif"].upper()

            if vcfrecord.ALT[0] is not None:
                    alt_alleles = []
                    for alt in vcfrecord.ALT:
                        n_repeats = float(str(alt)[1:-1])
                        int_repeats = math.floor(n_repeats)
                        extra = n_repeats - int_repeats
                        alt_allele = int_repeats*motif
                        idx = 0
                        while extra - 1/len(motif) >= 0:
                            alt_allele += motif[idx]
                            extra -= 1/len(motif)
                            idx += 1
                        alt_alleles.append(alt_allele)
            else:
                alt_alleles = [None]

            id = vcfrecord.ID
        else:
            raise ValueError("self.vcftype is the unexpected type {}".format(self.vcftype))

        return TRRecord(vcfrecord, ref_allele, alt_alleles, motif, id)

# TODO add __iter__ to go through samples?
# So we could do for sample in trrecord
class TRRecord:
    def __init__(self, vcfrecord, ref_allele, alt_alleles, motif, id):
        self.vcfrecord = vcfrecord
        self.id = id

        #alleles are stored as upper case strings with all the repeats written out
        #i.e. ACGACGACG
        #alleles may partial repeat copies or impurities in them
        #but will attempt to not contain any extra base pairs to either side of the repeat
        #but none of the statistics using this class should reflect that
        #because how these impurities are represented is not guaranteed 

        self.ref_allele = ref_allele

        #alt alleles is a list of alleles
        #if there are no alt alleles, it is the list [None]
        #as that's what pyvcf does
        self.alt_alleles = alt_alleles

        #motifs are as documented by the original 
        #or as derived from the original
        #they are not required to be canonical - so they may be 
        #cycles and/or reverse complements of the canonical motif
        self.motif = motif

    def GetStringGenotype(self, vcfsample):
        gts = vcfsample.gt_alleles
        gts_bases = [([self.ref_allele]+self.alt_alleles)[int(gt)] for gt in gts]
        return gts_bases

    def GetLengthGenotype(self, vcfsample):
        gts_bases = self.GetStringGenotype(vcfsample)
        return [float(len(item))/len(self.motif) for item in gts_bases]

    def GetAlleleCounts(self, uselength=True):
        allele_counts = {}
        for sample in self.vcfrecord:
            if sample.called:
                if uselength:
                    alleles = self.GetLengthGenotype(sample)
                else:
                    alleles = self.GetStringGenotype(sample)
                for a in alleles: allele_counts[a] += 1
        return allele_counts

    def __str__(self):
        string = "{} {} {} ".format(self.id, self.motif, self.ref_allele)
        if self.alt_alleles[0] is None:
            string += '.'
        else:
            string += ','.join(self.alt_alleles)
        return string

