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

import trtools.utils.common as common
import trtools.utils.mergeutils as mergeutils
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from trtools import __version__


NOCALLSTRING = "."

# Tool-specific fields to merge. (FIELDNAME, req'd). req'd is True if merge should
# fail when all records don't have identical values for that field
INFOFIELDS = {
    "gangstr": [("END", True), ("RU", True), ("PERIOD", True), ("REF", True), \
                ("EXPTHRESH", True), ("STUTTERUP", False), \
                ("STUTTERDOWN", False), ("STUTTERP", False)],
    "hipstr": [("INFRAME_PGEOM", False), ("INFRAME_UP", False), ("INFRAME_DOWN", False), \
               ("OUTFRAME_PGEOM", False), ("OUTFRAME_UP", False), ("OUTFRAME_DOWN", False), \
               ("BPDIFFS", False), ("START", True), ("END", True), ("PERIOD", True), \
               ("AN", False), ("REFAC", False), ("AC", False), ("NSKIP", False), \
               ("NFILT", False), ("DP", False), ("DSNP", False), ("DSTUTTER", False), \
               ("DFLANKINDEL", False)],
    "eh": [("END", True), ("REF", True), ("REPID", True), ("RL", True), \
           ("RU", True), ("SVTYPE", True), ("VARID", False)],
    "popstr": [("Motif", True)], # TODO ("RefLen", True) omitted. since it is marked as "A" incorrectly
    "advntr": [("END", True), ("VID", False), ("RU", True), ("RC", True)]
}

# Tool-specific format fields to merge
# Not all fields currently handled
# If not listed here, it is ignored
FORMATFIELDS = {
    "gangstr": ["DP","Q","REPCN","REPCI","RC","ENCLREADS","FLNKREADS","ML","INS","STDERR","QEXP"],
    "hipstr": ["GB","Q","PQ","DP","DSNP","PSNP","PDP","GLDIFF","DSTUTTER","DFLANKINDEL","AB","FS","DAB","ALLREADS","MALLREADS"],
    "eh": ["ADFL","ADIR","ADSP","LC","REPCI","REPCN","SO"],
    "popstr": ["AD","DP","PL"],
    "advntr": ["DP","SR","FR","ML"]
}


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

def WriteMergedHeader(vcfw, args, readers, cmd, vcftype):
    r"""Write merged header for VCFs in args.vcfs

    Also do some checks on the VCFs to make sure merging
    is appropriate.
    Return info and format fields to use

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
    vcftype : str
       Type of VCF files being merged

    Returns
    -------
    useinfo : list of (str, bool)
       List of (info field, required) to use downstream
    useformat: list of str
       List of format field strings to use downstream
    """
    # Check contigs the same for all readers
    contigs = readers[0].contigs
    for i in range(1, len(readers)):
        if readers[i].contigs != contigs:
            raise ValueError(
                "Different contigs found across VCF files. Make sure all "
                "files used the same reference. Consider using this "
                "command:\n\t"
                "bcftools reheader -f ref.fa.fai file.vcf.gz -o file_rh.vcf.gz")
    # Write VCF format, commands, and contigs
    vcfw.write("##fileformat=VCFv4.1\n")

    for r in readers:
        if "command" in r.metadata:
            vcfw.write("##command="+r.metadata["command"][0]+"\n")
    vcfw.write("##command="+cmd+"\n")
    for key,val in contigs.items():
        vcfw.write("##contig=<ID=%s,length=%s>\n"%(val.id, val.length))
    # Write INFO fields, different for each tool
    useinfo = []
    for (field, reqd) in INFOFIELDS[vcftype]:
        if field not in readers[0].infos:
            common.WARNING("Expected info field %s not found. Skipping"%field)
        else:
            vcfw.write(GetInfoString(readers[0].infos[field])+"\n")
            useinfo.append((field, reqd))
    # Write FORMAT fields, different for each tool
    useformat = []
    for field in FORMATFIELDS[vcftype]:
        if field not in readers[0].formats:
            common.WARNING("Expected format field %s not found. Skipping"%field)
        else:
            vcfw.write(GetFormatString(readers[0].formats[field])+"\n")
            useformat.append(field)
    # Write sample list
    samples = mergeutils.GetSamples(readers, usefilenames=args.update_sample_from_file)
    if len(samples) == 0:
        return None, None
    header_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcfw.write("#"+"\t".join(header_fields+samples)+"\n")
    return useinfo, useformat

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
                    alts.add(GetAlleleString(item))
    return sorted(list(alts), key=len)

def GetAlleleString(allele):
    """Get string representation of allele
    
    If it is a sequence, return upper case sequence
    If _SV type, return string representation

    Parameters
    ----------
    allele : ALT allele from vcf.Record
    
    Returns
    -------
    str_allele : str
       String representation of the allele
    """
    try:
        return allele.sequence.upper()
    except:
        return str(allele)

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
    if fail=False, only do something if we have a rule on how to handle that field

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
    if not fail: return None # TODO in future implement smart merging of select fields
    vals = set()
    for i in range(len(mergelist)):
        if mergelist[i]:
            if info_field in current_records[i].INFO:
                vals.add(current_records[i].INFO[info_field])
            else:
                raise ValueError("Missing info field %s"%info_field)
    if len(vals)==1:
        return "%s=%s"%(info_field, vals.pop())
    else:
        common.WARNING("Incompatible info field value %s"%info_field)
        return None

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
    # TODO check this is correct - need to do GetAlleleString?
    newgt = [alleles.index(gta.upper()) for gta in gt_alleles]
    return "/".join([str(item) for item in newgt])

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
    assert "GT" not in formats # since we will add that
    record_items = []
    for sample in record:
        sample_items = []
        if not sample.called:
            record_items.append(".")
            continue
        # Add GT
        sample_items.append(GetGT(sample.gt_bases.split(sample.gt_phase_char()), alleles))
        # Add rest of formats
        for fmt in formats:
            val = sample[fmt]
            if type(val)==list: val = ",".join([str(item) for item in val])
            sample_items.append(val)
        record_items.append(":".join([str(item) for item in sample_items]))
    return record_items

def MergeRecords(readers, current_records, mergelist, vcfw, args, useinfo, useformat):
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
    useinfo : list of (str, bool)
       List of (info field, required) to use downstream
    useformat: list of str
       List of format field strings to use downstream
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
    for (field, reqd) in useinfo:
        inf = GetInfoItem(current_records, mergelist, field, fail=reqd)
        if inf is not None:
            info_items.append(inf)
    info_items = [item for item in info_items if item is not None]
    output_items.append(";".join(info_items))
    # Set FORMAT - add GT to front
    output_items.append(":".join(["GT"]+useformat))
    # Set sample info
    alleles = [ref_allele]+alt_alleles
    for i in range(len(mergelist)):
        if mergelist[i]:
            output_items.extend(GetSampleInfo(current_records[i], alleles, useformat, args))
        else:
            output_items.extend([NOCALLSTRING]*len(readers[i].samples)) # NOCALL
    vcfw.write("\t".join(output_items)+"\n")

def getargs():  # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcfs", help="Comma-separated list of VCF files to merge (must be sorted, bgzipped and indexed)", type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    req_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VcfTypes.__members__], type=str, default="auto")
    ### Special merge options ###
    spec_group = parser.add_argument_group("Special merge options")
    spec_group.add_argument("--update-sample-from-file", help="Use file names, rather than sample header names, when merging", action="store_true")
    ### Optional arguments ###
    opt_group = parser.add_argument_group("Optional arguments")
    opt_group.add_argument("--verbose", help="Print out extra info", action="store_true")
    opt_group.add_argument("--quiet", help="Don't print out anything", action="store_true")
    ## Version argument ##
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    ### Parse args ###
    args = parser.parse_args()
    return args

def main(args):    
    ### Check and Load VCF files ###
    vcfs = args.vcfs.split(",")
    vcfreaders = utils.LoadReaders(args.vcfs.split(","), checkgz = True)
    if vcfreaders is None:
        return 1
    if len(vcfreaders) == 0: return 1
    contigs = vcfreaders[0].contigs
    # WriteMergedHeader will confirm that the list of contigs is the same for
    # each vcf, so just pulling it from one here is fine
    chroms = list(contigs)

    ### Check inferred type of each is the same
    vcftype = mergeutils.GetAndCheckVCFType(vcfreaders, args.vcftype)

    ### Set up VCF writer ###
    vcfw = open(args.out + ".vcf", "w")
    useinfo, useformat = WriteMergedHeader(vcfw, args, vcfreaders, " ".join(sys.argv), vcftype)
    if useinfo is None or useformat is None: return 1

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    # Check if contig ID is set in VCF header for all records
    done = mergeutils.DoneReading(current_records)
    while not done:
        if not all([r.CHROM in chroms for r in current_records if r is not None]):
            common.WARNING("Error: An input VCF file has a record with a CHROM not in its contig list.")
            return 1
        is_min = mergeutils.GetMinRecords(current_records, chroms)
        if args.verbose: mergeutils.DebugPrintRecordLocations(current_records, is_min)
        if mergeutils.CheckMin(is_min): return 1
        MergeRecords(vcfreaders, current_records, is_min, vcfw, args, useinfo, useformat)
        current_records = mergeutils.GetNextRecords(vcfreaders, current_records, is_min)
        done = mergeutils.DoneReading(current_records)
    return 0 

def run(): # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
