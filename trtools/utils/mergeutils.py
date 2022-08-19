"""
Utilities for reading multiple VCFs simulataneously
and keeping them in sync.
"""

import numpy as np
import os
import sys

from cyvcf2 import cyvcf2

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh

from typing import List, Union, Any, Optional, Callable, Tuple

CYVCF_RECORD = cyvcf2.Variant
CYVCF_READER = cyvcf2.VCF
COMPARABILITY_CALLBACK = Callable[[List[Optional[trh.TRRecord]], List[int], int], Union[bool, List[bool]]]


def LoadReaders(vcffiles: List[str], region: Optional[str] = None) -> List[CYVCF_READER]:
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
            raise ValueError("Make sure %s is bgzipped and indexed" % f)
        if not os.path.isfile(f):
            raise ValueError("Could not find VCF file %s" % f)
        if not os.path.isfile(f + ".tbi"):
            raise ValueError("Could not find VCF index %s.tbi" % f)
    readers = [vcf.Reader(open(f, "rb")) for f in vcffiles]
    if region is None:
        return readers
    else:
        return [r.fetch(region) for r in readers]


def GetSharedSamples(readers: List[CYVCF_READER]) -> List[str]:
    r"""Get list of samples used in all files being merged

    Parameters
    ----------
    readers : list of cyvcf.VCF objects

    Returns
    -------
    samples : list of str
        Samples present in all readers
    """
    if len(readers) == 0: return list()
    samples = set(readers[0].samples)
    if len(readers) == 1: return list(samples)
    for r in readers[1:]:
        samples = samples.intersection(set(r.samples))
    return list(samples)


def GetSamples(readers: List[CYVCF_READER], filenames: Optional[str] = None) -> List[str]:
    r"""Get list of samples used in all files being merged

    Parameters
    ----------
    readers : list of cyvcf2.VCF
    usefilenames : optional list of filenames
       If present, add filename to sample names.
       Useful if sample names overlap across files
       Must be the same length as readers

    Returns
    -------
    samples : list of str
       List of samples in merged list
    """
    samples = []
    if filenames:
        if len(readers) != len(filenames):
            raise ValueError("Must have same number of VCFs and VCF filenames")
        for r, filename in zip(readers, filenames):
            filename = filename.strip(".vcf.gz")
            samples += [filename + ":" + s for s in r.samples]
    else:
        for r in readers:
            if set(samples).intersection(set(r.samples)):
                raise ValueError("Found the same sample ID(s) in multiple input VCFs")
            samples += r.samples
    return samples


def GetAndCheckVCFType(vcfs: List[CYVCF_READER], vcftype: str) -> str:
    """Infer type of multiple VCFs

    If they are all the same, return that type
    If not, return error

    Parameters
    ----------
    vcfs: list of cyvcf2.VCF
      Multiple VCFs
    vcftype : str
      If it is unclear which of a few VCF callers produced the underlying
      VCFs (because the output markings of those VCF callers are similar)
      this string can be supplied by the user to choose from among
      those callers.

    Returns
    -------
    vcftype : str
      Inferred VCF type

    Raises
    ------
    TypeError
      If one of the VCFs does not look like it was produced by any supported TR
      caller, or if one of the VCFs looks like it could have been produced by
	  more than one supported TR caller and vcftype == 'auto', or if, for one
      of the VCFs, vcftype doesn't match any of the callers that could have
	  produced that VCF, or if the types of the VCFs don't match
    """
    types = []
    for vcf in vcfs:
        vcf_type = trh.InferVCFType(vcf, vcftype)
        types.append(vcf_type)
    if len(set(types)) == 1:
        return types[0]
    else:
        raise ValueError("VCF files are of mixed types.")


def GetChromOrder(r: CYVCF_RECORD, chroms: List[str]) -> Union[int, float]:
    r"""Get the chromosome order of a record

    Parameters
    ----------
    r : vcf.Record
    chroms : list of str
       Ordered list of chromosomes

    Returns
    -------
    order : int or float
       Index of r.CHROM in chroms, int
       Return np.inf if can't find r.CHROM, float
    """
    if r is None:
        return np.inf
    else:
        return chroms.index(r.CHROM)


def GetChromOrderEqual(chrom_order: Union[int, float], min_chrom: int) -> bool:
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


def GetPos(r: CYVCF_RECORD) -> Union[int, float]:
    r"""Get the position of a record

    Parameters
    ----------
    r : vcf.Record

    Returns
    -------
    pos : int
       If r is None, returns np.inf, which is a float
    """
    if r is None:
        return np.inf
    else:
        return r.POS


def CheckPos(record: CYVCF_RECORD, chrom: str, pos: int) -> bool:
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
    return record.CHROM == chrom and record.POS == pos


def GetMinRecords(record_list: List[Optional[trh.TRRecord]], chroms: List[str]) -> List[bool]:
    r"""Check if each record is next up in sort order

    Return a vector of boolean set to true if
    the record is in lowest sort order of all the records
    Use order in chroms to determine sort order of chromosomes

    Parameters
    ----------
    record_list : list of CYVCF_RECORD
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
        return [False] * len(record_list)
    return [CheckPos(r, chroms[min_chrom], min_pos) for r in record_list]


def default_callback(records: List[trh.TRRecord], chrom_order: List[int], min_chrom_index: int) -> bool:
    return True


def GetIncrementAndComparability(record_list: List[Optional[trh.TRRecord]],
                                 chroms: List[str],
                                 overlap_callback: COMPARABILITY_CALLBACK = default_callback) \
        -> Tuple[List[bool], Union[bool, List[bool]]]:

    r"""Get list that says which records should be skipped in the next
     iteration (increment), and whether they are all comparable / mergable
     The value of increment elements is determined by the (harmonized) position of corresponding records


    Parameters
    ----------
    record_list : trh.TRRecord
       list of current records from each file being merged

    chroms : list of str
       Ordered list of all chromosomes

    overlap_callback: Callable[[List[Optional[trh.TRRecord]], List[int], int], Union[bool, List[bool]]
        Function that calculates whether the records are comparable

    Returns
    -------
    increment : list of bool
       List or bools, where items are set to True when the record at the index of the item should be
       skipped during VCF file comparison.
    comparable: bool or list of bool
        Value, that determines whether current records are comparable / mergable, depending on the callback
    """
    chrom_order = [np.inf if r is None else chroms.index(r.chrom) for r in record_list]
    pos = [np.inf if r is None else r.pos for r in record_list]
    min_chrom_index = min(chrom_order)
    curr_pos=[pos[i] for i in range(len(chrom_order)) if chrom_order[i]==min_chrom_index]
    min_pos = min(curr_pos)
    increment = \
        [chrom_order[i] == min_chrom_index and pos[i] == min_pos and record_list[i] is not None
         for i in range(len(chrom_order))]
    comparable = overlap_callback(record_list, chrom_order, min_chrom_index)

    return increment, comparable


def DoneReading(records: List[Union[CYVCF_RECORD, trh.TRRecord]]) -> bool:
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


def DebugPrintRecordLocations(current_records: List[CYVCF_RECORD], is_min: List[bool]) -> None:
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
        chrom = current_records[i].CHROM if current_records[i] else None
        pos = current_records[i].POS if current_records[i] else None
        info.append("%s:%s:%s" % (chrom, pos, is_min[i]))
    common.MSG("\t".join(info) + "\n", debug=True)


def CheckMin(is_min: List[bool]) -> bool:
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
    if sum(is_min) == 0:
        raise ValueError("Unexpected error. Stuck in infinite loop and exiting.")
    return False


def GetNextRecords(readers: List[CYVCF_READER], current_records: List[CYVCF_RECORD], increment: List[bool]) \
        -> List[CYVCF_RECORD]:
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
            except StopIteration:
                new_records.append(None)
        else:
            new_records.append(current_records[i])
    return new_records


def InitReaders(readers: List[CYVCF_READER]) -> List[CYVCF_RECORD]:
    r"""Increment readers of each file

        Returns list of first records from list of readers.

        Parameters
        ----------
        readers : list of cyvcf2.VCF
           List of readers for all files being merged

        Returns
        -------
        list of vcf.Record
           List of next records for each file
        """
    return [next(reader) for reader in readers]
