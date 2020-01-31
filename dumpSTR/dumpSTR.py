#!/usr/bin/env python3
"""
Tool for filtering and QC of STR genotypes
"""

# Imports
import sys
import os

# Handle STRTools imports differently depending on where we're calling this from
if __name__ == "dumpSTR" or __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "strtools", "utils"))
    import filters # If running from source code
    import common
    import tr_harmonizer as trh
    import utils
else:
    import dumpSTR.filters as filters # If running as a package
    import strtools.utils.common as common
    import strtools.utils.tr_harmonizer as trh
    import strtools.utils.utils as utils
    
# Load external libraries
import argparse
import inspect
import sys
import vcf
from vcf.parser import _Filter
from vcf.parser import _Format
from vcf.parser import _Info

def MakeWriter(outfile, invcf, command):
    r"""Create a VCF writer with a dumpSTR header

    Adds a header line with the dumpSTR command used

    Parameters
    ----------
    outfile : str
       Name of the output file
    invcf : vcf.Reader object
       Input VCF. Used to grab header info
    command : str
       String command used to run dumpSTR

    Returns
    -------
    writer : vcf.Writer object
       VCF writer initialized with header of input VCF
    """
    invcf.metadata["command-DumpSTR"] = [command]
    writer = vcf.Writer(open(outfile, "w"), invcf)
    return writer

# TODO update with vcf type info
# TODO add fields for other tools
def CheckFilters(invcf, args, vcftype="auto"):
    r""" Perform checks on user input for filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments
    vcftype : str
        Specifies which tool this VCF came from

    Returns
    -------
    checks : bool
        Set to True if all filters look ok. 
        Set to False if filters are invalid 
    """
    if args.min_call_DP is not None:
        if args.min_call_DP < 0:
            common.WARNING("Invalid min_call_DP <0")
            return False
        if "DP" not in invcf.formats:
            common.WARNING("No DP FORMAT found")
            return False
    if args.max_call_DP is not None:
        if args.max_call_DP < 0:
            common.WARNING("Invalid min_call_DP <0")
            return False
        if args.min_call_DP is not None and args.max_call_DP <= args.min_call_DP:
            common.WARNING("--max-call-DP must be > --min-call-DP")
            return False
        if "DP" not in invcf.formats:
            common.WARNING("No DP FORMAT found")
            return False
    if args.min_call_Q is not None:
        if args.min_call_Q < 0 or args.min_call_Q > 1:
            common.WARNING("--min-call-Q must be between 0 and 1")
            return False
        if "Q" not in invcf.formats:
            common.WARNING("No Q FORMAT found")
            return False
    if args.max_call_flank_indel is not None:
        if args.max_call_flank_indel < 0 or args.max_call_flank_indel > 1:
            common.WARNING("--max-call-flank-indel must be between 0 and 1")
            return False
        if "DP" not in invcf.formats or "DFLANKINDEL" not in invcf.formats:
            common.WARNING("No DP or DFLANKINDEL FORMAT found")
            return False
    if args.max_call_stutter is not None:
        if args.max_call_stutter < 0 or args.max_call_stutter > 1:
            common.WARNING("--max-call-stutter must be between 0 and 1")
            return False
        if "DP" not in invcf.formats or "DSTUTTER" not in invcf.formats:
            common.WARNING("No DP or DSTUTTER FORMAT found")        
            return False
    if args.expansion_prob_het is not None:
        if args.expansion_prob_het < 0 or args.expansion_prob_het > 1:
            common.WARNING("--expansion-prob-het must be between 0 and 1")
            return False
    if args.expansion_prob_hom is not None:
        if args.expansion_prob_hom < 0 or args.expansion_prob_hom > 1:
            common.WARNING("--expansion-prob-hom must be between 0 and 1")
            return False
    if args.expansion_prob_total is not None:
        if args.expansion_prob_total < 0 or args.expansion_prob_total > 1:
            common.WARNING("--expansion-prob-total must be between 0 and 1")
            return False
    if args.require_support < 0:
        common.WARNING("--require-support must be >= 0")
        return False
    if args.require_support > 0 and args.readlen is None:
        common.WARNING("Using --require-support requires setting --readlen")
        return False
    if args.readlen is not None and args.readlen < 20:
        common.WARNING("--readlen must be an integer value >= 20")
        return False
    return True

def WriteLocLog(loc_info, fname):
    r"""Write locus-level features to log file

    Parameters
    ----------
    loc_info : dict of str->value
       Dictionary containing locus-level stats.
       Must have at least keys: 'totalcalls', 'PASS'
    fname : str
       Output log filename

    Returns
    -------
    success : bool
       Set to true if outputting the log was successful
    """
    if not os.path.exists(os.path.dirname(fname)):
        common.WARNING("Could not write log to %s"%fname)
        return False
    f = open(fname, "w")
    keys = list(loc_info.keys())
    if 'totalcalls' not in keys or 'PASS' not in keys:
        common.WARNING("Unexpected error occurred when writing log file. Missing required stats")
        return False
    keys.remove("totalcalls")
    if loc_info["PASS"] == 0:
        callrate = 0
    else:
        callrate = float(loc_info["totalcalls"])/loc_info["PASS"]
    f.write("MeanSamplesPerPassingSTR\t%s\n"%callrate)
    for k in keys:
        f.write("FILTER:%s\t%s\n"%(k, loc_info[k]))
    f.close()
    return True

def WriteSampLog(sample_info, reasons, fname):
    r"""Write sample-level features to log file
    
    Parameters
    ----------
    sample_info : dict of str->value
        Dictionary of statistics for each sample
    reasons: list of str
        List of possible call filter reasons
    fname : str
        Output filename

    Returns
    -------
    success : bool
       Set to true if outputting the log was successful
    """
    if not os.path.exists(os.path.dirname(fname)):
        common.WARNING("Could not write log to %s"%fname)
        return False
    if 'numcalls' not in sample_info or not 'totaldp' in sample_info:
        common.WARNING("Unexpected error occurred when writing sample log file. Missing required stats")
        return False        
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
    return True

def GetAllCallFilters():
    r"""List all possible call filters by extracting from filters module
    
    Returns
    -------
    reasons : list of str
        A list of call-level filter reasons
    """
    reasons = []
    for name, obj in inspect.getmembers(filters):
        if inspect.isclass(obj) and issubclass(obj, filters.Reason) and not obj is filters.Reason:
            reasons.append(obj.name)
    return reasons

def FilterCall(sample, call_filters):
    r"""Apply call-level filters and return filter reason.

    Parameters
    ----------
    sample : vcf._Call
       The call to be filtered
    call_filters : list of filters.Reason
       List of call-level filters

    Returns
    -------
    reasons : list of str
       List of string description of fitler reasons
    """
    reasons = []
    for cfilt in call_filters:
        if cfilt(sample) is not None: reasons.append(cfilt.GetReason())
    return reasons

def ApplyCallFilters(record, reader, call_filters, sample_info):
    r"""Apply call-level filters to a record

    Returns a new record with FILTER populated for each sample.
    Also updates sample_info with sample level stats
    
    Parameters
    ----------
    record : vcf._Record
       The record to apply filters to
    reader : vcf.Reader
       The vcf.Reader object we're reading from
    call_filters : list of filters.Reason
       List of call filters to apply
    sample_info : dict of str->value
       Dictionary of sample stats to keep updated

    Returns
    -------
    new_record : vcf._Record
       Modified record object with FILTER field set
    """
    if "FILTER" in record.FORMAT:
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
    else:
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':')+["FILTER"])
        record.FORMAT = record.FORMAT+":FILTER"
    for fmt in samp_fmt._fields:
        if fmt == "FILTER" and "FILTER" not in record.FORMAT:
            samp_fmt._types.append("String")
            samp_fmt._nums.append(1)
        else:
            entry_type = reader.formats[fmt].type
            entry_num  = reader.formats[fmt].num
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
    # Get data
    new_samples = []
    for sample in record:
        sampdat = []
        if sample['GT'] is None or sample['GT'] == "./." or sample['GT'] == ".":
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "FILTER":
                    sampdat.append("NOCALL")
                else: sampdat.append(sample[key])
            call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
            new_samples.append(call)
            continue
        filter_reasons = FilterCall(sample, call_filters)
        if len(filter_reasons) > 0:
            for r in filter_reasons:
                sample_info[sample.sample][r] += 1
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "GT":
                    sampdat.append("./.")
                else:
                    if key == "FILTER": sampdat.append(",".join(filter_reasons))
                    else: sampdat.append(None)
        else:
            sample_info[sample.sample]["numcalls"] += 1
            sample_info[sample.sample]["totaldp"] += sample["DP"]
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "FILTER": sampdat.append("PASS")
                else: sampdat.append(sample[key])
        call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
        new_samples.append(call)
    record.samples = new_samples
    return record

# TODO update with new filters for different VCFs
def BuildCallFilters(args):
    r"""Build list of locus-leve filters to include

    Parameters
    ----------
    args : argparse namespace
       User input arguments used to decide on filters
    
    Returns
    -------
    filter_list : list of filters.Filter
       List of call-level filters to apply
    """
    filter_list = []
    if args.min_call_DP is not None:
        filter_list.append(filters.LowCallDepth(args.min_call_DP))
    if args.max_call_DP is not None:
        filter_list.append(filters.HighCallDepth(args.max_call_DP))
    if args.min_call_Q is not None:
        filter_list.append(filters.LowCallQ(args.min_call_Q))
    if args.max_call_flank_indel is not None:
        filter_list.append(filters.CallFlankIndels(args.max_call_flank_indel))
    if args.max_call_stutter is not None:
        filter_list.append(filters.CallStutter(args.max_call_stutter))
    if args.min_supp_reads > 0:
        filter_list.append(filters.CallMinSuppReads(args.min_supp_reads))
    if args.expansion_prob_het is not None:
        filter_list.append(filters.ProbHet(args.expansion_prob_het))
    if args.expansion_prob_hom is not None:
        filter_list.append(filters.ProbHom(args.expansion_prob_hom))
    if args.expansion_prob_total is not None:
        filter_list.append(filters.ProbTotal(args.expansion_prob_total))
    if args.filter_span_only:
        filter_list.append(filters.SpanOnly())
    if args.filter_spanbound_only:
        filter_list.append(filters.SpanBoundOnly())
    if args.filter_badCI:
        filter_list.append(filters.BadCI())
    if args.require_support > 0:
        filter_list.append(filters.RequireSupport(args.require_support, args.readlen))
    return filter_list

def BuildLocusFilters(args):
    r"""Build list of locus-level fitlers to include.

    These filters should in general not be tool specific

    Paramters
    ---------
    args : argparse namespace
       User input arguments used to decide on filters

    Returns
    -------
    filter_list : list of filters.Filter
       List of locus-level filters
    success : bool
       Return False if we had an issue building filters
    """
    filter_list = []
    if args.min_locus_callrate is not None:
        filter_list.append(filters.Filter_MinLocusCallrate(args.min_locus_callrate))
    if args.min_locus_hwep is not None:
        filter_list.append(filters.Filter_MinLocusHWEP(args.min_locus_hwep, args.use_length))
    if args.min_locus_het is not None:
        filter_list.append(filters.Filter_MinLocusHet(args.min_locus_het, args.use_length))
    if args.max_locus_het is not None:
        filter_list.append(filters.Filter_MaxLocusHet(args.max_locus_het, args.use_length))
    if args.filter_hrun is not None:
        filter_list.append(filters.Filter_LocusHrun())
    if args.filter_regions is not None:
        filter_region_files = args.filter_regions.split(",")
        if args.filter_regions_names is not None:
            filter_region_names = args.filter_regions_names.split(",")
            if len(filter_region_names) != len(filter_region_files):
                common.WARNING("Length of --filter-regions-names must match --filter-regions\n")
                return filter_list, False
        else: filter_region_names = [str(item) for item in list(range(len(filter_region_files)))]
        for i in range(len(filter_region_names)):
            filter_list.append(filters.create_region_filter(filter_region_names[i], filter_region_files[i]))
    return filter_list, True

def getargs():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)

    locus_group = parser.add_argument_group("Locus-level filters")
    locus_group.add_argument("--min-locus-callrate", help="Minimum locus call rate", type=float)
    locus_group.add_argument("--min-locus-hwep", help="Filter loci failing HWE at this p-value threshold", type=float)
    locus_group.add_argument("--min-locus-het", help="Minimum locus heterozygosity", type=float)
    locus_group.add_argument("--max-locus-het", help="Maximum locus heterozygosity", type=float)
    locus_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    locus_group.add_argument("--filter-regions", help="Comma-separated list of BED files of regions to filter", type=str)
    locus_group.add_argument("--filter-regions-names", help="Comma-separated list of filter names for each BED filter file", type=str)
    locus_group.add_argument("--filter-hrun", help="Filter STRs with long homopolymer runs.", action="store_true")
    locus_group.add_argument("--drop-filtered", help="Drop filtered records from output", action="store_true")

    #### General call-level filters (HipSTR or GangSTR)
    call_group = parser.add_argument_group("General call-level filters")
    call_group.add_argument("--min-call-DP", help="Minimum call coverage", type=int)
    call_group.add_argument("--max-call-DP", help="Maximum call coverage", type=int)
    call_group.add_argument("--min-call-Q", help="Minimum call quality score", type=float)    

    #### HipSTR
    hipstr_call_group = parser.add_argument_group("Call-level filters specific to HipSTR output")
    hipstr_call_group.add_argument("--max-call-flank-indel", help="Maximum call flank indel rate", type=float)
    hipstr_call_group.add_argument("--max-call-stutter", help="Maximum call stutter rate", type=float)
    hipstr_call_group.add_argument("--min-supp-reads", help="Minimum supporting reads for each allele", default=0, type=int)
    
    #### GangSTR
    gangstr_call_group = parser.add_argument_group("Call-level filters specific to GangSTR output")
    gangstr_call_group.add_argument("--expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    gangstr_call_group.add_argument("--expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    gangstr_call_group.add_argument("--expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)
    gangstr_call_group.add_argument("--filter-span-only", help="Filter out all calls that only have spanning read support", action="store_true")
    gangstr_call_group.add_argument("--filter-spanbound-only", help="Filter out all reads except spanning and bounding", action="store_true")
    gangstr_call_group.add_argument("--filter-badCI", help="Filter regions where the ML estimate is not in the CI", action="store_true")
    gangstr_call_group.add_argument("--require-support", help="Require each allele call to have at least n supporting reads", type=int, default=0)
    gangstr_call_group.add_argument("--readlen", help="Read length used (bp). Required if using --require-support", type=int)

    debug_group = parser.add_argument_group("Debugging parameters")
    debug_group.add_argument("--num-records", help="Only process this many records", type=int)
    debug_group.add_argument("--die-on-warning", help="Quit if a record can't be parsed", action="store_true")
    debug_group.add_argument("--verbose", help="Print out extra info", action="store_true")

    args = parser.parse_args()
    return args

def main(args=None):
    if args is None:
        args = getargs()
    # Load VCF file
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    invcf = vcf.Reader(filename=args.vcf)

    # Set up filter list and check for errors
    if not CheckFilters(invcf, args): return 1
    invcf.filters = {}
    filter_list, build_success = BuildLocusFilters(args)
    if not build_success: return 1

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
            if args.die_on_warning: return 1
        except StopIteration: break
        if args.verbose:
            common.MSG("Processing %s:%s"%(record.CHROM, record.POS))
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
            trrecord = trh.TRRecord(record, record.REF, record.ALT, "", "") # TODO handle different VCF types
            # Recalculate locus-level INFO fields
            record.INFO["HRUN"] = utils.GetHomopolymerRun(record.REF)
            if record.num_called > 0:
                allele_freqs = trrecord.GetAlleleFreqs(uselength=args.use_length)
                genotype_counts = trrecord.GetGenotypeCounts(uselength=args.use_length)
                record.INFO["HET"] = utils.GetHeterozygosity(allele_freqs)
                record.INFO["HWEP"] = utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)
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

    return 0

if __name__ == "__main__":
    # Set up args
    args = getargs()
    # Run main function
    retcode = main(args)
    sys.exit(retcode)

