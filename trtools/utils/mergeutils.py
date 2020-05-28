"""
Utilities for reading multiple VCFs simulataneously
and keeping them in sync.
"""

import numpy as np
import os
import sys
import vcf

if __name__ == "mergeutils" or __package__ is None:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import common
    import utils
    import tr_harmonizer as trh
else:
    import trtools.utils.common as common # pragma: no cover
    import trtools.utils.utils as utils # pragma: no cover
    import trtools.utils.tr_harmonizer as trh # pragma: no cover

def LoadReaders(vcffiles, region=None):
    r"""Return list of VCF readers

    Parameters
    ----------
    vcffiles : list of str
        List of VCF files to merge
    region : str, optional
        Chrom:start-end to restrict to
    
    Returns
    -------
    readers : list of vcf.Reader
        VCF readers list for all files to merge
    """
    for f in vcffiles:
        if not f.endswith(".vcf.gz") and not f.endswith(".vcf.bgz"):
            raise ValueError("Make sure %s is bgzipped and indexed"%f)
        if not os.path.isfile(f):
            raise ValueError("Could not find VCF file %s"%f)
        if not os.path.isfile(f+".tbi"):
            raise ValueError("Could not find VCF index %s.tbi"%f)
    readers = [vcf.Reader(open(f, "rb")) for f in vcffiles]
    if region is None:
        return readers
    else: return [r.fetch(region) for r in readers]

def GetSharedSamples(readers, usefilenames=False):
    r"""Get list of samples used in all files being merged

    Parameters
    ----------
    readers : list of vcf.Reader objects

    Returns
    -------
    samples : list of str
        Samples present in all readers
    """
    if len(readers) == 0: return set()
    samples = set(readers[0].samples)
    if len(readers) == 1: return samples
    for r in readers[1:]:
        samples = samples.intersection(set(r.samples))
    return samples

def GetSamples(readers, usefilenames=False):
    r"""Get list of samples used in all files being merged

    Parameters
    ----------
    readers : list of vcf.Reader objects
    usefilenames : bool, optional
       If True, add filename to sample names.
       Useful if sample names overlap across files

    Returns
    -------
    samples : list of str
       List of samples in merged list
    """
    samples = []
    for r in readers:
        if usefilenames:
            samples = samples + [r.filename.strip(".vcf.gz")+":"+ s for s in r.samples]
        else: samples = samples + r.samples
    if len(set(samples))!=len(samples):
        common.WARNING("Duplicate samples found.")
        return []
    return samples

def InferAndCheckVCFType(vcfreaders, vcftype):
    """Infer vcf type of readers

    If vcftype is "auto", try to infer types of each reader.
    If it is not auto, just return that string.
    If they are all the same, return that type
    If not, return error

    Parameters
    ----------
    vcfreaders : list of vcf.Reader
      Readers being merged
    vcftype : str
      Type of VCF

    Returns
    -------
    vcftype : str
      Inferred VCF type
    """
    types = []
    for reader in vcfreaders:
        reader_type = utils.InferVCFType(reader, vcftype)
        if reader_type == None:
            raise ValueError("Error in detecting vcf type. Please make sure the vcf file is from a supported variant calling software.")
        types.append(reader_type)
    if len(set(types)) == 1: 
        return types[0]
    else: raise ValueError("VCF files are of mixed types.")

def GetChromOrder(r, chroms):
    r"""Get the chromosome order of a record

    Parameters
    ----------
    r : vcf.Record
    chroms : list of str
       Ordered list of chromosomes

    Returns
    -------
    order : int
       Index of r.CHROM in chroms
       Return np.inf if can't find r.CHROM
    """
    if r is None: return np.inf
    else: return chroms.index(r.CHROM)

def GetChromOrderEqual(chrom_order, min_chrom):
    r"""Check chrom order

    Parameters
    ----------
    chrom_order : int
       Chromosome order
    min_chrom : int
       Current chromosome order

    Returns
    -------
    equal : bool
       Return True if chrom_order==min_chrom and chrom_order != np.inf
    """
    if chrom_order == np.inf: return False
    return chrom_order == min_chrom

def GetPos(r):
    r"""Get the position of a record

    Parameters
    ----------
    r : vcf.Record

    Returns
    -------
    pos : int
       If r is None, returns np.inf
    """
    if r is None: return np.inf
    else: return r.POS

def CheckPos(record, chrom, pos):
    r"""Check a record is at the specified position
    
    Parameters
    ----------
    r : vcf.Record
       VCF Record being checked
    chrom : str
       Chromosome name
    pos : int
       Chromosome position

    Returns
    -------
    check : bool
       Return True if the current record is at this position
    """
    if record is None: return False
    return record.CHROM==chrom and record.POS==pos

def GetMinRecords(record_list, chroms):
    r"""Check if each record is next up in sort order

    Return a vector of boolean set to true if
    the record is in lowest sort order of all the records
    Use order in chroms to determine sort order of chromosomes

    Parameters
    ----------
    record_list : list of vcf.Record
       list of current records from each file being merged
    chroms : list of str
       Ordered list of all chromosomes
    
    Returns
    -------
    checks : list of bool
       Set to True for records that are first in sort order
    """
    chrom_order = [GetChromOrder(r, chroms) for r in record_list]
    pos = [GetPos(r) for r in record_list]
    min_chrom = min(chrom_order)
    allpos = [pos[i] for i in range(len(pos)) if GetChromOrderEqual(chrom_order[i], min_chrom)]
    if len(allpos) > 0:
        min_pos = min(allpos)
    else:
        return [False]*len(record_list)
    return [CheckPos(r, chroms[min_chrom], min_pos) for r in record_list]

def DoneReading(records):
    r"""Check if all records are at the end of the file

    Parameters
    ----------
    records : list of vcf.Record
       List of records from files to merge
 
    Returns
    -------
    check : list of bool
       Set to True if all record is None
       indicating we're done reading the file
    """
    return all([item is None for item in records])

def PrintCurrentRecords(current_records, is_min):
    r"""Debug function to print current records for each file

    Parameters
    ----------
    current_records : list of vcf.Record
       List of current records from merged files
    is_min : list of bool
       List of check for if record is first in sort order
    """
    info = []
    for i in range(len(is_min)):
        try:
            chrom = current_records[i].CHROM
            pos = current_records[i].POS
        except:
            chrom = None
            pos = None
            common.WARNING("Missing CHROM and POS in record.")
        info.append("%s:%s:%s"%(chrom, pos, is_min[i]))
    common.MSG("\t".join(info)+"\n")

def CheckMin(is_min):
    r"""Check if we're progressing through VCFs

    Parameters
    ----------
    is_min : list of bool
        List indicating if each record is first in sort order

    Returns
    -------
    check : bool
        Set to True if something went wrong
    """
    if sum(is_min)==0:
        raise ValueError("Unexpected error. Stuck in infinite loop and exiting.")
    return False

def GetNextRecords(readers, current_records, increment):
    r"""Increment readers of each file

    Increment readers[i] if increment[i] set to true
    Else keep current_records[i]

    Parameters
    ----------
    readers : list of vcf.Reader
       List of readers for all files being merged
    current_records : list of vcf.Record
       List of current records for all readers
    increment : list of bool
       List indicating if each file should be incremented

    Returns
    -------
    new_records : list of vcf.Record
       List of next records for each file
    """
    new_records = []
    for i in range(len(readers)):
        if increment[i]:
            try:
                new_records.append(next(readers[i]))
            except: new_records.append(None)
        else: new_records.append(current_records[i])
    return new_records

