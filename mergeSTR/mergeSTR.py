#!/usr/bin/env python3
"""
Tool for merging STR VCF files from GangSTR
"""

# Load external libraries
import argparse
import os
import numpy as np
import sys
import vcf

# Load local libraries
if __name__ == "mergeSTR" or __name__ == '__main__' or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import common
    import utils
else: # pragma: no cover
    import strtools.utils.common as common  # pragma: no cover
    import strtools.utils.utils as utils  # pragma: no cover

NOCALLSTRING = "."
FORMATFIELDS = ["GT","DP","Q","REPCN","REPCI","RC","ML","INS","STDERR","ENCLREADS","FLNKREADS","QEXP"]

def LoadReaders(vcffiles):
    r"""Return list of VCF readers

    Parameters
    ----------
    vcffiles : list of str
        List of VCF files to merge
    
    Returns
    -------
    readers : list of vcf.Reader
        VCF readers list for all files to merge
    """
    for f in vcffiles:
        if not f.endswith(".vcf.gz"):
            raise ValueError("Make sure %s is bgzipped and indexed"%f)
        if not os.path.isfile(f):
            raise ValueError("Could not find VCF file %s"%f)
        if not os.path.isfile(f+".tbi"):
            raise ValueError("Could not find VCF index %s.tbi"%f)
    return [vcf.Reader(open(f, "rb")) for f in vcffiles]

def GetInfoString(info):
    r"""Create a VCF INFO header string

    Parameters
    ----------
    info : PyVCF info field

    Returns
    -------
    infostring : str
       Formatted info string for header
    """
    return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">'%(info.id, info.num, info.type, info.desc)

def GetFormatString(fmt):
    r"""Create a VCF FORMAT header string

    Parameters
    ----------
    fmt : PyVCF format field
    
    Returns
    -------
    formatstring : str
       Formatted format string for header
    """
    return '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">'%(fmt.id, fmt.num, fmt.type, fmt.desc)

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
        raise ValueError("Duplicate samples found.")
    return samples

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

# TODO needs editing
def WriteMergedHeader(vcfw, args, readers, cmd):
    r"""Write merged header for VCFs in args.vcfs

    Also do some checks on the VCFs to make sure merging
    is appropriate

    Parameters
    ----------
    vcfw : file object
       Writer to write the merged VCF
    args : argparse namespace
       Contains user options
    readers : list of vcf.Reader
       List of readers to merge
    cmd : str
       Command used to call this program    
    """
    # Check contigs the same for all readers
    contigs = readers[0].contigs
    for i in range(1, len(readers)):
        if readers[i].contigs != contigs:
            raise ValueError("Different contigs found across VCF files. Make sure all files used the same reference")
    # Write VCF format, commands, and contigs
    vcfw.write("##fileformat=VCFv4.1\n")
    for r in readers: vcfw.write("##command="+r.metadata["command"][0]+"\n")
    vcfw.write("##command="+cmd+"\n")
    for key,val in contigs.items():
        vcfw.write("##contig=<ID=%s,length=%s>\n"%(val.id, val.length))
    ######## TODO change this part based on inferred VCF type
    # Write GangSTR specific INFO fields - TODO automatically infer from each tool. Avoid fields not implemented yet
    for field in ["END", "PERIOD", "RU", "REF","STUTTERUP","STUTTERDOWN","STUTTERP","EXPTHRESH"]:
        vcfw.write(GetInfoString(readers[0].infos[field])+"\n")
    # Write GangSTR specific FORMAT fields
    rmfields = []
    for field in FORMATFIELDS:
        if field not in readers[0].formats:
            rmfields.append(field)
        else: vcfw.write(GetFormatString(readers[0].formats[field])+"\n")
    for field in rmfields: FORMATFIELDS.remove(field)
    ############## end TODO
    # Write sample list
    samples = GetSamples(readers, usefilenames=args.update_sample_from_file)
    header_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcfw.write("#"+"\t".join(header_fields+samples)+"\n")

def GetRefAllele(current_records, mergelist):
    r"""Get reference allele for a set of records

    Parameters
    ----------
    current_records : list of vcf.Record
       List of records being merged
    mergelist : list of bool
       Indicates whether each record is included in merge

    Returns
    -------
    ref : str
       Reference allele string
    """
    refs = []
    chrom = ""
    pos = -1
    for i in range(len(mergelist)):
        if mergelist[i]:
            chrom = current_records[i].CHROM
            pos = current_records[i].POS
            refs.append(current_records[i].REF.upper())
    if len(set(refs)) != 1:
        raise ValueError("Conflicting refs found at %s:%s"%(chrom, pos))
    return refs[0]

def GetAltAlleles(current_records, mergelist):
    r"""Get list of alt alleles
    
    Parameters
    ----------
    current_records : list of vcf.Record
       List of records being merged
    mergelist : list of bool
       Indicates whether each record is included in merge

    Returns
    -------
    alts : list of str
       List of alternate allele strings
    """
    alts = set()
    for i in range(len(mergelist)):
        if mergelist[i]:
            ralts = current_records[i].ALT
            for item in ralts:
                if item is not None and item:
                    alts.add(item.sequence.upper())
    return sorted(list(alts), key=len)

def GetID(idval):
    r"""Get the ID for a a record

    If not set, output "."

    Parameters
    ----------
    idval : str
       ID of the record

    Returns
    -------
    idval : str
       Return ID. if None, return "."
    """
    if idval is None: return "."
    else: return idval

def GetInfoItem(current_records, mergelist, info_field, fail=True):
    """Get INFO item for a group of records

    Make sure it's the same across merged records
    if fail=True, die if items not the same.

    Parameters
    ----------
    current_records : list of vcf.Record
       List of records being merged
    mergelist : list of bool
       List of indicators of whether to merge each record
    info_field : str
       INFO field being merged
    fail : bool
       If True, throw error if fields don't have same value

    Returns
    -------
    infostring : str
       INFO string to add (key=value)
    """
    vals = set()
    for i in range(len(mergelist)):
        if mergelist[i]:
            vals.add(current_records[i].INFO[info_field])
    if len(vals)==1: return "%s=%s"%(info_field, vals.pop())
    else:
        raise ValueError("WARNING more than one value found for %s"%info_field)

def GetGT(gt_alleles, alleles):
    r"""Update GT field based on ref/alt alleles

    Parameters
    ----------
    gt_alleles : list of str
       List of GT allele strings
    alleles : list of str
       List of REF + ALT alleles
    
    Returns
    -------
    newgt : list of str
       List of new GT field based on updated allele list
    """
    newgt = [alleles.index(gta.upper()) for gta in gt_alleles]
    return "/".join([str(item) for item in newgt])

# TODO needs editing
def GetSampleInfo(record, alleles, formats, args):
    r"""Output sample FORMAT info

    Parameters
    ----------
    record : vcf.Record
       VCF record being summarized
    alleles : list of str
       List of REF + ALT alleles
    formats : list of str
       List of VCF FORMAT items
    args : argparse namespace
       User options

    Returns
    -------
    sampleinfo : str
       FORMAT fields for the sample
    """
    easy_fmts = FORMATFIELDS
    record_items = []
    for sample in record:
        sample_items = []
        if not sample.called:
            record_items.append(".")
            continue
        for fmt in formats:
            if fmt == "GT":
                sample_items.append(GetGT(sample.gt_bases.split(sample.gt_phase_char()), alleles))
            elif fmt in easy_fmts:
                try:
                    val = sample[fmt]
                    if type(val)==list: val = ",".join([str(item) for item in val])
                except:
                    val = NOCALLSTRING
                sample_items.append(val)
        record_items.append(":".join([str(item) for item in sample_items]))
    return record_items

# TODO update 
def MergeRecords(readers, current_records, mergelist, vcfw, args):
    r"""Merge records from different files

    Merge all records with indicator set to True in mergelist
    Output merged record to vcfw

    Parameters
    ----------
    readers : list of vcf.Reader
       List of readers being merged
    current_records : list of vcf.Record
       List of current records for each reader
    mergelist : list of bool
       Indicates whether to include each reader in merge
    vcfw : file
       File to write output to
    args : argparse namespace
       Contains user options
    """
    output_items = []
    use_ind = [i for i in range(len(mergelist)) if mergelist[i]]
    if len(use_ind)==0: return
    alt_alleles = GetAltAlleles(current_records, mergelist)
    ref_allele = GetRefAllele(current_records, mergelist)
    # Set common fields
    output_items.append(current_records[use_ind[0]].CHROM) # CHROM
    output_items.append(str(current_records[use_ind[0]].POS)) # POS
    output_items.append(GetID(current_records[use_ind[0]].ID)) # ID
    output_items.append(ref_allele) # REF
    if len(alt_alleles) == 0:
        output_items.append(".")
    else: output_items.append(",".join(alt_alleles)) # ALT
    output_items.append(".") # QUAL
    output_items.append(".") # FILTER
    # Set INFO
    info_items = []
    for info_field in ["END", "PERIOD", "RU", "REF"]:
        info_items.append(GetInfoItem(current_records, mergelist, info_field))
    for info_field in ["STUTTERUP","STUTTERDOWN","STUTTERP","EXPTHRESH"]:
        info_items.append(GetInfoItem(current_records, mergelist, info_field, fail=False))
    info_items = [item for item in info_items if item is not None]
    output_items.append(";".join(info_items))
    # Set FORMAT
    formats = FORMATFIELDS[:]
    output_items.append(":".join(formats))
    # Set sample info
    alleles = [ref_allele]+alt_alleles
    for i in range(len(mergelist)):
        if mergelist[i]:
            output_items.extend(GetSampleInfo(current_records[i], alleles, formats, args))
        else:
            output_items.extend([NOCALLSTRING]*len(readers[i].samples)) # NOCALL
    vcfw.write("\t".join(output_items)+"\n")

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
        info.append("%s:%s:%s"%(chrom, pos, is_min[i]))
    sys.stderr.write("\t".join(info)+"\n")

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
        sys.stderr.write("Unexpected error. Stuck in infinite loop and exiting.")
        return True
    return False

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcfs", help="Comma-separated list of VCF files to merge (must be sorted, bgzipped and indexed)", type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    req_group.add_argument("--vcftype", help="Options=%s"%trh.VCFTYPES.__members__, type=str, default="auto")
    ### Special merge options ###
    spec_group = parser.add_argument_group("Special merge options")
    spec_group.add_argument("--update-sample-from-file", help="Use file names, rather than sample header names, when merging", action="store_true")
    ### Optional arguments ###
    opt_group = parser.add_argument_group("Optional arguments")
    opt_group.add_argument("--verbose", help="Print out extra info", action="store_true")
    opt_group.add_argument("--quiet", help="Don't print out anything", action="store_true")
    ### Parse args ###
    args = parser.parse_args()
    return args

def main(args):
    ### Check VCF files ###
    vcfs = args.vcfs.split(",")
    if not os.path.exists(vcfs[0]):
        common.WARNING("%s does not exist"%args.vcfs)
        return 1
    if not os.path.exists(vcfs[1]):
        common.WARNING("%s does not exist"%args.vcfs)
        return 1
    ### Load readers ###
    try:
        vcfreaders = LoadReaders(args.vcfs.split(","))
    except ValueError: return 1
    if len(vcfreaders) == 0: return 1
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    ### Set up VCF writer ###
    vcfw = open(args.out + ".vcf", "w")
    try:
        WriteMergedHeader(vcfw, args, vcfreaders, " ".join(sys.argv))
    except ValueError: return 1

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = GetMinRecords(current_records, chroms)
    done = DoneReading(current_records)
    while not done:
        if args.verbose: PrintCurrentRecords(current_records, is_min)
        if CheckMin(is_min): return 1
        MergeRecords(vcfreaders, current_records, is_min, vcfw, args)
        current_records = GetNextRecords(vcfreaders, current_records, is_min)
        is_min = GetMinRecords(current_records, chroms)
        done = DoneReading(current_records)
    return 0 

if __name__ == "__main__":  # pragma: no cover
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)
