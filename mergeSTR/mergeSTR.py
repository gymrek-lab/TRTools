#!/usr/bin/env python3

"""
Tool for merging STR VCF files from GangSTR

"""

# Test command:
# ./mergeSTR/mergeSTR.py --vcfs ERR1955393.17.sorted.vcf.gz,ERR1955394.17.sorted.vcf.gz --out test

### GangSTR merging ###
# By default, merge the following FORMAT fields in the logical way:
# GT, DP, Q, REPCN, REPCI, RC, ML, INS, STDERR

# By default, require the following INFO fields to be the same across files:
# END, PERIOD, RU, REF

# Stutter model (INFO:STUTTERUP, STUTTERDOWN, and STUTTERP) included if it is the same across files.
# FORMAT:QEXP included if EXPTHRESH is the same across files

# Special options required for:
#- Merging GGL (affects INFO:GRID and FORMAT:GGL)

# Ignores filters. For now this should be done before filtering with DumpSTR

# Load external libraries
import argparse
import os
import numpy as np
import sys
import vcf

# Load local libraries
import strtools.utils.common as common
import strtools.utils.utils as utils

NOCALLSTRING = "."

def GetInfoString(info):
    return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">'%(info.id, info.num, info.type, info.desc)

def GetFormatString(fmt):
    return '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">'%(fmt.id, fmt.num, fmt.type, fmt.desc)

def GetSamples(readers, usefilenames=False):
    samples = []
    for r in readers:
        if usefilenames:
            samples = samples + [r.filename.strip(".vcf.gz")+":"+ s for s in r.samples]
        else: samples = samples + r.samples
    if len(set(samples))!=len(samples):
        common.ERROR("Duplicate samples found. Quitting")
    return samples

def MergeGRID(current_records, mergelist):
    """
    Merge the INFO/GRID field
    """
    return (-1, -1) # TODO
    
def WriteMergedHeader(vcfw, args, readers, cmd):
    """
    Write merged header for VCFs in args.vcfs
    Also do some checks on the VCFs to make sure merging
    is appropriate
    """
    # Check contigs the same for all readers
    contigs = readers[0].contigs
    for i in range(1, len(readers)):
        if readers[i].contigs != contigs:
            common.ERROR("Different contigs found across VCF files. Make sure all files used the same reference")
    # Write VCF format, commands, and contigs
    vcfw.write("##fileformat=VCFv4.1\n")
    for r in readers: vcfw.write("##command="+r.metadata["command"][0]+"\n")
    vcfw.write("##command="+cmd+"\n")
    for key,val in contigs.items():
        vcfw.write("##contig=<ID=%s,length=%s>\n"%(val.id, val.length))
    # Write GangSTR specific INFO fields
    for field in ["END", "PERIOD", "RU", "REF","STUTTERP","STUTTERDOWN","STUTTERP","EXPTHRESH"]:
        vcfw.write(GetInfoString(readers[0].infos[field])+"\n")
    if args.merge_ggl: vcfw.write(GetInfoString(readers[0].infos["GRID"])+"\n")
    # Write GangSTR specific FORMAT fields
    for field in ["GT", "DP", "Q", "REPCN", "REPCI", "RC", "ML", "INS", "STDERR", "QEXP"]:
        vcfw.write(GetFormatString(readers[0].formats[field])+"\n")
    if args.merge_ggl: vcfw.write(GetFormatString(readers[0].formats["GGL"])+"\n")
    # Write sample list
    samples=GetSamples(readers, usefilenames=args.update_sample_from_file)
    header_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcfw.write("#"+"\t".join(header_fields+samples)+"\n")

def GetChromOrder(r, chroms):
    if r is None: return np.inf
    else: return chroms.index(r.CHROM)

def GetPos(r):
    if r is None: return np.inf
    else: return r.POS

def CheckPos(record, chrom, pos):
    if record is None: return False
    return record.CHROM==chrom and record.POS==pos

def GetChromOrderEqual(chrom_order, min_chrom):
    if chrom_order == np.inf: return False
    return chrom_order == min_chrom

def GetMinRecords(record_list, chroms, debug=True):
    """
    Return a vector of boolean set to true if 
    the record is in lowest sort order of all the records
    Use order in chroms to determine sort order of chromosomes
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

def GetRefAllele(current_records, mergelist):
    refs = []
    chrom = ""
    pos = -1
    for i in range(len(mergelist)):
        if mergelist[i]:
            chrom = current_records[i].CHROM
            pos = current_records[i].POS
            refs.append(current_records[i].REF.upper())
    if len(set(refs)) != 1:
        common.ERROR("Conflicting refs found at %s:%s"%(chrom, pos))
    return refs[0]

def GetAltAlleles(current_records, mergelist):
    alts = set()
    for i in range(len(mergelist)):
        if mergelist[i]:
            ralts = current_records[i].ALT
            for item in ralts:
                if item is not None and item:
                    alts.add(item.sequence.upper())
    return sorted(list(alts), key=len)

def GetID(idval):
    if idval is None: return "."
    else: return idval

def GetInfoItem(current_records, mergelist, info_field, fail=True):
    """
    Get info item. Make sure it's the same across merged records
    if fail=True, die if items not the same
    """
    vals = set()
    for i in range(len(mergelist)):
        if mergelist[i]:
            vals.add(current_records[i].INFO[info_field])
    if len(vals)==1: return "%s=%s"%(info_field, vals.pop())
    else:
        if fail: common.ERROR("More than one value found for %s"%info_field)
        sys.stderr.write("WARNING more than one value found for %s"%info_field)
        return None

def GetGT(gt_alleles, alleles):
    """
    Update GT field based on ref/alt alleles
    """
    newgt = [alleles.index(gta.upper()) for gta in gt_alleles]
    return "/".join([str(item) for item in newgt])

def GetSampleInfo(record, alleles, formats, use_qexp, args):
    """
    Output sample info
    """
    easy_fmts = ["DP","Q","REPCN","REPCI","RC","ML","INS","STDERR"]
    if use_qexp: easy_fmts.append("QEXP")
    record_items = []
    for sample in record:
        sample_items = []
        if not sample.called:
            record_items.append(".")
            continue
        for fmt in formats:
            if fmt == "GT":
                sample_items.append(GetGT(sample.gt_bases.split(sample.gt_phase_char()), alleles))
            if fmt == "GGL":
                if args.merge_ggl:
                    sample_items.append(".") # TODO merge this for real
            if fmt in easy_fmts:
                val = sample[fmt]
                if type(val)==list: val = ",".join([str(item) for item in val])
                sample_items.append(val)
        record_items.append(":".join([str(item) for item in sample_items]))
    return record_items

def MergeRecords(readers, current_records, mergelist, vcfw, args):
    """
    Merge all records with indicator set to True in mergelist
    Output merged record to vcfw
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
    if args.merge_ggl:
        info_items.append("GRID=%s,%s"%MergeGRID(current_records, mergelist))
    info_items = [item for item in info_items if item is not None]
    output_items.append(";".join(info_items))
    # should we use qexp?
    use_qexp = ("EXPTHRESH" in [item.split("=")[0] for item in info_items])
    # Set FORMAT
    formats = ["GT","DP","Q","REPCN","REPCI","RC","ML","INS","STDERR"]
    if use_qexp: formats.append("QEXP")
    if args.merge_ggl: formats.append("GGL")
    output_items.append(":".join(formats))
    # Set sample info
    alleles = [ref_allele]+alt_alleles
    for i in range(len(mergelist)):
        if mergelist[i]:
            output_items.extend(GetSampleInfo(current_records[i], alleles, formats, use_qexp, args))
        else:
            output_items.extend([NOCALLSTRING]*len(readers[i].samples)) # NOCALL
    vcfw.write("\t".join(output_items)+"\n")

def LoadReaders(vcffiles):
    """
    Return list of VCF readers
    """
    if len(vcffiles) == 0:
        common.ERROR("No VCF files found")
    for f in vcffiles:
        if not f.endswith(".vcf.gz"):
            common.ERROR("Make sure %s is bgzipped and indexed"%f)
        if not os.path.isfile(f):
            common.ERROR("Could not find VCF file %s"%f)
        if not os.path.isfile(f+".tbi"):
            common.ERROR("Could not find VCF index %s.tbi"%f)
    return [vcf.Reader(open(f, "rb")) for f in vcffiles]

def DoneReading(records):
    """
    Check if all records are at the end of the file
    """
    return all([item is None for item in records])

def GetNextRecords(readers, current_records, increment):
    """
    Increment readers[i] if increment[i] set to true
    Else keep current_records[i]
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
    if sum(is_min)==0:
        common.ERROR("Unexpected error. Stuck in infinite loop and exiting.")

def main():
    parser = argparse.ArgumentParser(__doc__)
    ### Required arguments ###
    req_group = parser.add_argument_group("Required arguments")
    req_group.add_argument("--vcfs", help="Comma-separated list of VCF files to merge (must be sorted, bgzipped and indexed)", type=str, required=True)
    req_group.add_argument("--out", help="Prefix to name output files", type=str, required=True)
    ### Special merge options ###
    spec_group = parser.add_argument_group("Special merge options")
    spec_group.add_argument("--update-sample-from-file", help="Use file names, rather than sample header names, when merging", action="store_true")
    spec_group.add_argument("--merge-ggl", help="Merge GGL fields", action="store_true")
    ### Optional arguments ###
    opt_group = parser.add_argument_group("Optional arguments")
    opt_group.add_argument("--verbose", help="Print out extra info", action="store_true")
    opt_group.add_argument("--quiet", help="Don't print out anything", action="store_true")
    ### Parse args ###
    args = parser.parse_args()
    if args.merge_ggl: common.ERROR("--merge-ggl not implemented yet") # TODO remove

    ### Load readers ###
    vcfreaders = LoadReaders(args.vcfs.split(","))
    contigs = vcfreaders[0].contigs
    chroms = list(contigs)

    ### Set up VCF writer ###
    vcfw = open(args.out + ".vcf", "w")
    WriteMergedHeader(vcfw, args, vcfreaders, " ".join(sys.argv))

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = GetMinRecords(current_records, chroms, debug=args.verbose)
    done = DoneReading(current_records)
    while not done:
        if args.verbose: PrintCurrentRecords(current_records, is_min)
        CheckMin(is_min)
        MergeRecords(vcfreaders, current_records, is_min, vcfw, args)
        current_records = GetNextRecords(vcfreaders, current_records, is_min)
        is_min = GetMinRecords(current_records, chroms)
        done = DoneReading(current_records)

if __name__ == "__main__":
    main()
