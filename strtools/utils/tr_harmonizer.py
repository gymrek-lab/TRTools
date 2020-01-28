"""
Utilities for harmonizing tandem repeat VCF records
across various formats output by TR genotyping tools
"""

import enum
import math
import numpy as np
import os
import sys

if __name__ == "tr_harmonizer":
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import utils
else:
    import strtools.utils.utils as utils

# List of supported VCF types
# TODO: add Beagle
# TODO: add support for tool version numbers
VCFTYPES =  enum.Enum('Types', ['gangstr', 'advntr', 'hipstr', 'eh', 'popstr'])

class TRRecordHarmonizer:
    """
    Class for harmonizing TR VCF records across tools.

    The main purpose of this class is to infer which tool
    a VCF came from, and appropriately convert its records
    to TRRecord objects.

    Parameters
    ----------
    vcffile : PyVCF reader instance
    vcftype : {'auto', 'gangstr', 'advntr', 'hipstr', 'eh', 'popstr'}, optional
       Type of the VCF file. Default='auto'.
       If vcftype=='auto', attempts to infer the type.

    Attributes
    ----------
    vcffile : PyVCF reader instance
    vcftype : enum
       Type of the VCF file. Must be included in trh.VCFTYPES
    """

    def __init__(self, vcffile, vcftype="auto"):
        self.vcffile = vcffile
        if vcftype == "auto":
            self.vcftype = self.InferVCFType()
        else:
            if vcftype not in VCFTYPES.__members__:
                raise ValueError("{} is not an excepted TR vcf type. Expected one of {}".format(vcftype, list(VCFTYPES.__members__)))
            self.vcftype = VCFTYPES[vcftype]

    def InferVCFType(self):
        """
        Infer the genotyping tool used to create the VCF

        When we can, infer from header metadata.
        Otherwise, try to infer the type from the ALT field.

        Returns
        -------
        vcftype : enum
           Type of the VCF file. Must be included in trh.VCFTYPES

        Notes
        -----
        Some notes on the pecularities of each VCF type:
        GangSTR:
           Does not include sequence impurities or partial repeats.
           Full REF and ALT strings given
        HipSTR:
           Full REF and ALT strings given
           May contain sequence impurities and partial repeats
           Sometimes includes non-repeat flanks in the original genotypes,
           which are trimmed off during harmonization.
           If the alt alleles have differently sized flanks than the ref allele
           then those alt alleles will be improperly trimmed.
           HipSTR reports the motif length but not the motif itself,
           so the motif is inferred from the sequences. This is usually
           but not always the correct motif.
        adVNTR:
           Full REF and ALT strings given.
           Uses its own ID convention that can change from run to run.
           May contain impurities or partial repeats?
        popSTR:
           Includes full REF string, which can contain impurities.
           ALT field gives only repeat counts.
           Inferred ALT alleles start at the beginning of the motif
           (even if its in the wrong part of the cycle), never contain
           impurities, and in the case of non-integer repeat counts,
           we include extra characters from the motif until we would exceed
           the reported count. Since repeat counts are rounded to the nearest
           tenth, we could miss/add extra base pairs that weren't there.
        ExpansionHunter:
           Not implemented yet
        """
        possible_vcf_types = set()
        for key, value in self.vcffile.metadata.items():
            # Sometimes value is a list of values, other times it is the value itself
            if type(value) == list:
                values = value
            else:
                values = {value}
            for value in values:
                if key.upper() == 'COMMAND' and 'GANGSTR' in value.upper():
                    possible_vcf_types.add(VCFTYPES.gangstr)
                if key.upper() == 'COMMAND' and 'HIPSTR' in value.upper():
                    possible_vcf_types.add(VCFTYPES.hipstr)
                if key.upper() == 'SOURCE' and 'ADVNTR' in value.upper():
                    possible_vcf_types.add(VCFTYPES.advntr)
                if key.upper() == 'SOURCE' and 'POPSTR' in value.upper():
                    possible_vcf_types.add(VCFTYPES.popstr)
                # TODO expansion hunter - it doesn't document itself well

        if len(possible_vcf_types) == 0:
            raise ValueError('Could not identify the type of this vcf')

        if len(possible_vcf_types) > 1:
            raise ValueError('Confused - this vcf looks like it could have been any of the types: {}'.format(possible_vcf_types))

        return next(iter(possible_vcf_types))

    def __iter__(self):
        for record in self.vcffile:
            yield self.HarmonizeRecord(record)

    def HarmonizeRecord(self, vcfrecord):
        """
        Harmonize VCF record to the allele string representation

        This is the representation used by GangSTR, AdVNTR, HipSTR.

        Parameters
        ----------
        vcfrecord : pyvcf.Record
            A PyVCF record object

        Returns
        -------
        trrecord : TRRecord
            A harmonized TRRecord object
        """
        ref_allele = None
        alt_alleles = None
        motif = None

        if self.vcftype == VCFTYPES.gangstr:
            ref_allele = vcfrecord.REF.upper()
            if vcfrecord.ALT[0] is not None:
                alt_alleles = []
                for alt in vcfrecord.ALT:
                    alt_alleles.append(str(alt).upper())
            else:
                alt_alleles = [None]
            motif = vcfrecord.INFO["RU"].upper()
            record_id = vcfrecord.ID

        elif self.vcftype == VCFTYPES.hipstr:
            # Remove the flanking variants
            # Note that if there are both flanking variants
            # and some of the alterante alleles have insertions/deletions in the flanks
            # then this might inadvertantly elongate/truncate those alternate alleles
            pos = int(vcfrecord.POS)
            start_offset = int(vcfrecord.INFO['START']) - pos
            pos_end_offset = int(vcfrecord.INFO['END']) - pos
            neg_end_offset = pos_end_offset + 1 - len(vcfrecord.REF)

            # neg_end_offset is the number of flanking non repeat bp to remove from the
            # end of each allele
            # e.g. 'AAAT'[0:-1] == 'AAA'
            # however, if neg_end_offset == 0, then we would get
            # 'AAAA'[0:0] == '' which is not the intent
            # so we need an if statement to instead write 'AAAA'[0:]
            # which gives us 'AAAA'
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

            # Get the motif. Hipstr doesn't tell us this explicitly, so figure it out
            motif = utils.InferRepeatSequence(ref_allele[start_offset:], vcfrecord.INFO["PERIOD"])
            record_id = vcfrecord.ID

        elif self.vcftype == VCFTYPES.advntr:
            ref_allele = vcfrecord.REF.upper()
            if vcfrecord.ALT[0] is not None:
                alt_alleles = []
                for alt in vcfrecord.ALT:
                    alt_alleles.append(str(alt).upper())
            else:
                alt_alleles = [None]

            motif = vcfrecord.INFO["RU"].upper()
            record_id = vcfrecord.INFO["VID"]
        
        elif self.vcftype == VCFTYPES.eh:
            raise NotImplementedError("Haven't implemented expansion hunter harmonization yet")

        elif self.vcftype == VCFTYPES.popstr:
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

            record_id = vcfrecord.ID
        else:
            raise ValueError("self.vcftype is the unexpected type {}".format(self.vcftype))

        return TRRecord(vcfrecord, ref_allele, alt_alleles, motif, record_id)

class TRRecord:
    """
    Analagous to VCF Record, with TR fields harmonized

    Allows downstream functions to be agnostic to the
    genotyping tool used to create the record.

    Parameters
    ----------
    vcfrecord : Pyvcf.Record
       PyVCF Record object used to generate the metadata
    ref_allele : str
       Reference allele string
    alt_alleles : list of str
       List of alternate allele strings
    motif : str
       Repeat unit
    record_id : str
       Identifier for the record

    Attributes
    ----------
    vcfrecord : Pyvcf.Record
       PyVCF Record object used to generate the metadata
    ref_allele : str
       Reference allele string. Gets converted to uppercase
       e.g. ACGACGACG
    alt_allele : list of str
       List of alternate allele strings
    motif : str
       Repeat unit
    record_id : str
       Identifier for the record
    
    Notes
    -----
    Alleles are stored as upper case strings with all the repeats written out.
    Alleles may be partial repeat copies or impurities in them.
    But will attempt to not contain any extra base pairs to either side of the repeat
    None of the statistics using this class should reflect that
    Because how these impurities are represented is not guaranteed 
    """

    def __init__(self, vcfrecord, ref_allele, alt_alleles, motif, record_id):
        self.vcfrecord = vcfrecord
        self.ref_allele = ref_allele
        self.alt_alleles = alt_alleles
        self.motif = motif
        self.record_id = record_id
        if not self.CheckRecord():
            raise ValueError("Invalid TRRecord.")

    def CheckRecord(self):
        """
        Check that valid allele lists are valid given sample genotypes

        Returns
        -------
        is_valid : bool
            Returns True if the record is valid. Otherwise return False.
        """
        if len(self.alt_alleles) == 0 or None in self.alt_alleles:
            num_alleles = 1
        else: num_alleles = 1 + len(self.alt_alleles)
        for sample in self.vcfrecord:
            if not sample.called: continue
            gts = sample.gt_alleles
            for al in gts:
                if int(al) > num_alleles-1: return False
        return True

    def GetStringGenotype(self, vcfsample):
        """
        Get the genotype of a VCF sample

        Parameters
        ----------
        vcfsample : VCF.Call object
            The VCF.Call object for the sample

        Returns
        -------
        genotype : list of str
            The string representation of the call genotype
        """
        gts = vcfsample.gt_alleles
        gts_bases = [str(([self.ref_allele]+self.alt_alleles)[int(gt)]) for gt in gts]
        return gts_bases

    def GetLengthGenotype(self, vcfsample):
        """
        Get the genotype of a VCF sample

        Parameters
        ----------
        vcfsample : VCF.Call object
            The VCF.Call object for the sample

        Returns
        -------
        genotype : list of float
            The float representation of the call genotype
            Each item gives the repeat copy number of each allele
        """
        gts_bases = self.GetStringGenotype(vcfsample)
        return [float(len(item))/len(self.motif) for item in gts_bases]

    def GetGenotypeCounts(self, samplelist=[], uselength=True):
        """
        Get the counts of each genotype for a record

        Parameters
        ----------
        samplelist : list of str, optional
            List of samples to include when computing counts
        uselength : bool, optional
            If True, represent alleles a lengths
            else represent as strings

        Returns
        -------
        genotype_counts: dict of (tuple: int)
            Gives the count of each genotype.
            Genotypes are represented as tuples of alleles.
        """
        genotype_counts = {}
        for sample in self.vcfrecord:
            if len(samplelist) > 0 and sample.sample not in samplelist: continue
            if sample.called:
                if uselength:
                    alleles = tuple(self.GetLengthGenotype(sample))
                else:
                    alleles = tuple(self.GetStringGenotype(sample))
                genotype_counts[alleles] = genotype_counts.get(alleles, 0) + 1
        return genotype_counts

    def GetAlleleCounts(self, samplelist=[], uselength=True):
        """
        Get the counts of each allele for a record

        Parameters
        ----------
        samplelist : list of str, optional
            List of samples to include when computing counts
        uselength : bool, optional
            If True, represent alleles a lengths
            else represent as strings

        Returns
        -------
        allele_counts: dict of (object: int)
            Gives the count of each allele.
            Alleles may be represented as floats or strings
        """
        allele_counts = {}
        for sample in self.vcfrecord:
            if len(samplelist) > 0 and sample.sample not in samplelist: continue
            if sample.called:
                if uselength:
                    alleles = self.GetLengthGenotype(sample)
                else:
                    alleles = self.GetStringGenotype(sample)
                for a in alleles:
                    allele_counts[a] = allele_counts.get(a, 0) + 1
        return allele_counts

    def GetAlleleFreqs(self, samplelist=[], uselength=True):
        """
        Get the frequencies of each allele for a record

        Parameters
        ----------
        samplelist : list of str, optional
            List of samples to include when computing frequencies
        uselength : bool, optional
            If True, represent alleles a lengths
            else represent as strings

        Returns
        -------
        allele_freqs: dict of (object: float)
            Gives the frequency of each allele.
            Alleles may be represented as floats or strings
        """
        allele_counts = self.GetAlleleCounts(uselength=uselength, samplelist=samplelist)
        total = float(sum(allele_counts.values()))
        allele_freqs = {}
        for key in allele_counts.keys():
            allele_freqs[key] = allele_counts[key]/total
        return allele_freqs

    def GetMaxAllele(self, samplelist=[]):
        """
        Get the maximum allele length called in a record

        Parameters
        ----------
        samplelist : list of str, optional
            List of samples to include when computing frequencies

        Returns
        -------
        maxallele : float
            The maximum allele length called (in number of repeat units)
        """
        alleles = self.GetAlleleCounts(uselength=True, samplelist=samplelist).keys()
        if len(alleles) == 0: return np.nan
        return max(alleles)

    def __str__(self):
        string = "{} {} {} ".format(self.record_id, self.motif, self.ref_allele)
        if self.alt_alleles[0] is None:
            string += '.'
        else:
            string += ','.join(self.alt_alleles)
        return string

