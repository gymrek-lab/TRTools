#!/usr/bin/env python3

"""
Tool for merging STR VCF files from GangSTR

"""

# Test command:
# ./mergeSTR/mergeSTR.py --vcfs test1.vcf.gz,test2.vcf.gz --out test

### GangSTR merging ###
# By default, merge the following FORMAT fields in the logical way:
# GT, DP, Q, REPCN, REPCI, RC, ML, INS, STDERR

# By default, require the following INFO fields to be the same across files:
# END, PERIOD, RU, REF

# Stutter model (INFO:STUTTERUP, STUTTERDOWN, and STUTTERP) included if it is the same across files.
# FORMAT:QEXP included if EXPTHRESH is the same across files

# Special options required for:
#- Merging GGL (affects INFO:GRID and FORMAT:GGL)

# Load external libraries
import argparse
import vcf

# Load local libraries
import strtools.utils.common as common
import strtools.utils.utils as utils

def WriteMergedHeader(vcfw, args):
    """
    Write merged header for VCFs in args.vcfs
    Also do some checks on the VCFs to make sure merging
    is appropriate
    """
    # TODO
    pass

def GetMinRecords(record_list):
    """
    Return a vector of boolean set to true if 
    the record is in lowest sort order of all the records
    """
    return [False]*len(record_list) # TODO

def MergeRecords(current_records, mergelist, vcfw, args):
    """
    Merge all records with indicator set to True in mergelist
    Output merged record to vcfw
    """
    pass # TODO

def LoadReaders(vcffiles):
    """
    Return list of VCF readers
    """
    return []
    # TODO

def DoneReading(records):
    """
    Check if all records are at the end of the file
    """
    return True # TODO

def GetNextRecords(readers, current_records, increment):
    """
    Increment readers[i] if increment[i] set to true
    Else keep current_records[i]
    """
    new_records = []
    for i in range(len(readers)):
        if increment[i]:
            new_records.append(next(readers[i]))
        else: new_records.append(current_records[i])
    return new_records

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
    ### Parser args ###
    args = parser.parse_args()

    ### Set up VCF writer ###
    vcfw = open(args.out + ".vcf", "w")
    WriteMergedHeader(vcfw, args)

    ### Load readers ###
    vcfreaders = LoadReaders(args.vcfs.split(","))

    ### Walk through sorted readers, merging records as we go ###
    current_records = [next(reader) for reader in vcfreaders]
    is_min = GetMinRecords(current_records)
    done = DoneReading(current_records)
    while not done:
        MergeRecords(current_records, is_min, vcfw, args)
        current_records = GetNextRecords(vcfreaders, current_records, is_min)
        is_min = GetMinRecords(current_records)
        done = DoneReading(current_records)

if __name__ == "__main__":
    main()
