#!/usr/bin/env python3

"""
Tool for filtering and QC of STR genotypes

"""

# Imports
import sys
import os

# Load local libraries
import strtools.utils.common as common
import strtools.utils.utils as utils
import statSTR.stats as stats

# Load external libraries
import argparse
import inspect
import sys
import vcf
from vcf.parser import _Filter
from vcf.parser import _Format
from vcf.parser import _Info

def GetSamples(readers, usefilenames=False):
    samples = []
    for r in readers:
        if usefilenames:
            samples = samples + [r.filename.strip(".vcf.gz")+":"+ s for s in r.samples]
        else: samples = samples + r.samples
    if len(set(samples))!=len(samples):
        common.ERROR("Duplicate samples found. Quitting")
    return samples


def WriteMergedHeader(vcfw, args, readers, cmd):
    # Write sample list
    samples=GetSamples(readers, usefilenames=args.update_sample_from_file)
    header_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    vcfw.write("#"+"\t".join(header_fields+samples)+"\n")


def MakeWriter(outfile, invcf, command):
    invcf.metadata["command-statSTR"] = [command]
    writer = vcf.Writer(open(outfile, "w"), invcf)
    return writer

def CheckStats(invcf, args):
    """
    Perform checks on user input for filters

    Input:
    - invcf (vcf.Reader)
    - args (argparse namespace)

    Exit program if checks fail
    """
    if args.stats is not None:
        statSet = set(args.stats)
        invalidStat = [x for x in ['thresh','mean','sd'] if x not in statSet]
        if len(invalidStat) > 0: 
            common.ERROR("Invalid stat " + invalidStat )

def WriteLocLog(loc_info, fname):
    """
    Write locus-level features to log file

    Input:
    - loc_info (dict): dictionary with locus-level stats
    - fname (str): output filename
    """
    f = open(fname, "w")
    keys = list(loc_info.keys())
    keys.remove("totalcalls")
    if loc_info["PASS"] == 0: callrate = 0
    else: callrate = loc_info["totalcalls"]*1.0/loc_info["PASS"]
    f.write("MeanSamplesPerPassingSTR\t%s\n"%callrate)
    for k in keys:
        f.write("Stats:%s\t%s\n"%(k, loc_info[k]))
    f.close()

def WriteSampLog(sample_info, reasons, fname):
    """
    Write sample-level features to log file

    Input:
    - sample_info (dict): dictionary of stats for each sample
    - reasons (list<str>): list of possible feature reasons
    - fname (str): output filename
    """
    header = ["sample", "numcalls","meanDP"] + reasons
    f = open(fname, "w")
    f.write("\t".join(header)+"\n")
    for s in sample_info:
        numcalls = sample_info[s]["numcalls"]
        if numcalls > 0:
            meancov = sample_info[s]["totaldp"]*1.0/numcalls
        else: meancov = 0
        items = [s, numcalls, meancov]
        for r in reasons: items.append(sample_info[s][r])
        f.write("\t".join([str(item) for item in items])+"\n")
    f.close()

def GetAllCallFilters():
    """
    List all possible call filters by
    extracting from filters module

    Output:
    - reasons (list<str>): list of call-level filter reasons
    """
    reasons = []
    for name, obj in inspect.getmembers(stats):
        if inspect.isclass(obj) and issubclass(obj, stat.Reason) and not obj is stat.Reason:
            reasons.append(obj.name)
    return reasons


def BuildCallStats(args):
    """
    Build list of locus-level filters to include

    Input:
    - args (namespace from parser.parse_args)
    
    Output:
    - cdict (list<filters.Filter>): list of call-level filters
    """
    cdict = []
    if args.thresh is not None:
        cdict.append(stats.thresh(args.thresh))
    if args.mean is not None:
        cdict.append(stats.mean(args.mean))
    if args.sd is not None:
        cdict.append(stats.sd(args.sd))
    return cdict

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)

    #### GangSTR
    gangstr_call_group = parser.add_argument_group("Call-level filters specific to GangSTR output")
    gangstr_call_group.add_argument("--expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    gangstr_call_group.add_argument("--expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)

    debug_group = parser.add_argument_group("Debugging parameters")
    debug_group.add_argument("--num-records", help="Only process this many records", type=int)

    args = parser.parse_args()
    # Load VCF file
    invcf = vcf.Reader(open(args.vcf, "rb"))

    # Set up filter list
    CheckStats(invcf, args)
    invcf.stats = {}
    stat_list = BuildLocusStats(args)
    for f in filter_list:
        short_doc = f.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        invcf.filters[f.filter_name()] = _Filter(f.filter_name(), short_doc)
    call_filters = BuildCallFilters(args)
    # Add new FORMAT fields
    if "FILTER" not in invcf.formats:
        invcf.formats["FILTER"] = _Format("FILTER", 1, "String", "Call-level filter")

    # Add new INFO fields
    invcf.infos["AC"] = _Info("AC", -1, "Integer", "Alternate allele counts", source=None, version=None)
    invcf.infos["REFAC"] = _Info("REFAC", 1, "Integer", "Reference allele count", source=None, version=None)
    invcf.infos["HET"] = _Info("HET", 1, "Float", "Heterozygosity", source=None, version=None)
    invcf.infos["HWEP"] = _Info("HWEP", 1, "Float", "HWE p-value for obs. vs. exp het rate", source=None, version=None)
    invcf.infos["HRUN"] = _Info("HRUN", 1, "Integer", "Length of longest homopolymer run", source=None, version=None)

    # Set up output files
    outvcf = MakeWriter(args.out + ".vcf", invcf, " ".join(sys.argv))

    # Set up sample info
    all_reasons = GetAllCallFilters()
    sample_info = {}
    for s in invcf.samples:
        sample_info[s] = {"numcalls": 0, "totaldp": 0}
        for r in all_reasons: sample_info[s][r]  = 0

    # Set up locus info
    loc_info = {"totalcalls": 0, "PASS": 0} 
    for filt in filter_list: loc_info[filt.filter_name()] = 0

    # Go through each record
    record_counter = 0
    while True:
        try:
            record = next(invcf)
        except IndexError:
            common.WARNING("Skipping TR that couldn't be parsed by PyVCF. Check VCF format")
            continue
        except StopIteration: break
        record_counter += 1
        if args.num_records is not None and record_counter > args.num_records: break
        # Call-level filters
        record = ApplyCallFilters(record, invcf, call_filters, sample_info)

        # Locus-level filters
        record.FILTER = None
        output_record = True
        for filt in filter_list:
            if filt(record) == None: continue
            if args.drop_filtered:
                output_record = False
                break
            record.add_filter(filt.filter_name())
            loc_info[filt.filter_name()] += 1
        if args.drop_filtered:
            if record.call_rate == 0: output_record = False
        if output_record:
            # Recalculate locus-level INFO fields
            record.INFO["HRUN"] = utils.GetHomopolymerRun(record.REF)
            if record.num_called > 0:
                if args.use_length:
                    record.INFO["HET"] = utils.GetLengthHet(record)
                else: record.INFO["HET"] = record.heterozygosity
                record.INFO["HWEP"] = utils.GetSTRHWE(record, uselength=args.use_length)
                record.INFO["AC"] = [int(item*(2*record.num_called)) for item in record.aaf]
                record.INFO["REFAC"] = int((1-sum(record.aaf))*(2*record.num_called))
            else:
                record.INFO["HET"] = -1
                record.INFO["HWEP"] = -1
                record.INFO["AC"] = [0]*len(record.ALT)
                record.INFO["REFAC"] = 0
            # Recalc filter
            if record.FILTER is None and not args.drop_filtered:
                record.FILTER = "PASS"
                loc_info["PASS"] += 1
                loc_info["totalcalls"] += record.num_called
            # Output the record
            outvcf.write_record(record)

    # Output log info
    WriteSampLog(sample_info, all_reasons, args.out + ".samplog.tab")
    WriteLocLog(loc_info, args.out+".loclog.tab")

if __name__ == "__main__":
    main()

