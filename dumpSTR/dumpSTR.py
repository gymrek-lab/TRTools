#!/usr/bin/env python3
"""
Tool for filtering and QC of STR genotypes
"""

# Load external libraries
import argparse
import inspect
import os
# Imports
import sys

import vcf
from vcf.parser import _Filter, _Format, _Info

# Handle TRTools imports differently depending on where we're calling this from
if __name__ == "dumpSTR" or __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trtools", "utils"))
    import filters # If running from source code
    import common
    import tr_harmonizer as trh
    import utils
    import version
else:  # pragma: no cover
    import dumpSTR.filters as filters  # pragma: no cover
    import trtools.utils.common as common # pragma: no cover
    import trtools.utils.tr_harmonizer as trh # pragma: no cover
    import trtools.utils.utils as utils # pragma: no cover
    import trtools.utils.version as version
__version__ = version.__version__


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
       Set to None if we had a problem writing the file
    """
    invcf.metadata["command-DumpSTR"] = [command]
    try:
        writer = vcf.Writer(open(outfile, "w"), invcf)
    except OSError:
        writer = None
    return writer

def CheckLocusFilters(args, vcftype):
    r"""Perform checks on user inputs for locus-level filters

    Parameters
    ----------
    args : argparse namespace
        Contains user arguments
    vcftype : enum.
        Specifies which tool this VCF came from.
        Must be included in trh.VCFTYPES

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.min_locus_hwep is not None:
        if args.min_locus_hwep < 0 or args.min_locus_hwep > 1:
            common.WARNING("Invalid --min-locus-hwep. Must be between 0 and 1")
            return False
    if args.min_locus_het is not None:
        if args.min_locus_het < 0 or args.min_locus_het > 1:
            common.WARNING("Invalid --min-locus-het. Must be between 0 and 1")
            return False
    if args.max_locus_het is not None:
        if args.max_locus_het < 0 or args.max_locus_het > 1:
            common.WARNING("Invalid --max-locus-het. Must be between 0 and 1")
            return False
    if args.min_locus_het is not None and args.max_locus_het is not None:
        if args.max_locus_het < args.min_locus_het:
            common.WARNING("Cannot have --max-locus-het less than --min-locus-het")
            return False
    if args.use_length and vcftype not in [trh.VCFTYPES["hipstr"]]:
        common.WARNING("--use-length is only meaningful for HipSTR, which reports sequence level differences.")
    if args.filter_hrun and vcftype not in [trh.VCFTYPES["hipstr"]]:
        common.WARNING("--filter-run only relevant to HipSTR files. This filter will have no effect.")
    if args.filter_regions is not None:
        if args.filter_regions_names is not None:
            filter_region_files = args.filter_regions.split(",")
            filter_region_names = args.filter_regions_names.split(",")
            if len(filter_region_names) != len(filter_region_files):
                common.WARNING("Length of --filter-regions-names must match --filter-regions.")
                return False
    return True

def CheckHipSTRFilters(invcf, args):
    r"""Check HipSTR call-level filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.hipstr_max_call_flank_indel is not None:
        if args.hipstr_max_call_flank_indel < 0 or args.hipstr_max_call_flank_indel > 1:
            common.WARNING("--hipstr-max-call-flank-indel must be between 0 and 1")
            return False
        assert "DP" in invcf.formats and "DFLANKINDEL" in invcf.formats # should always be true
    if args.hipstr_max_call_stutter is not None:
        if args.hipstr_max_call_stutter < 0 or args.hipstr_max_call_stutter > 1:
            common.WARNING("--hipstr-max-call-stutter must be between 0 and 1")
            return False
        assert "DP" in invcf.formats and "DSTUTTER" in invcf.formats # should always be true
    if args.hipstr_min_supp_reads is not None:
        if args.hipstr_min_supp_reads < 0:
            common.WARNING("--hipstr-min-supp-reads must be >= 0")
            return False
        assert "ALLREADS" in invcf.formats and "GB" in invcf.formats
    if args.hipstr_min_call_DP is not None:
        if args.hipstr_min_call_DP < 0:
            common.WARNING("--hipstr-min-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.hipstr_max_call_DP is not None:
        if args.hipstr_max_call_DP < 0:
            common.WARNING("--hipstr-max-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.hipstr_min_call_DP is not None and args.hipstr_max_call_DP is not None:
        if args.hipstr_max_call_DP < args.hipstr_min_call_DP:
            common.WARNING("--hipstr-max-call-DP must be >= --hipstr-min-call-DP")
            return False
    if args.hipstr_min_call_Q is not None:
        if args.hipstr_min_call_Q < 0 or args.hipstr_min_call_Q > 1:
            common.WARNING("--hipstr-min-call-Q must be between 0 and 1")
            return False
        assert "Q" in invcf.formats
    return True

def CheckGangSTRFilters(invcf, args):
    r"""Check GangSTR call-level filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.gangstr_min_call_DP is not None:
        if args.gangstr_min_call_DP < 0:
            common.WARNING("--gangstr-min-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.gangstr_max_call_DP is not None:
        if args.gangstr_max_call_DP < 0:
            common.WARNING("--gangstr-max-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.gangstr_min_call_DP is not None and args.gangstr_max_call_DP is not None:
        if args.gangstr_max_call_DP < args.gangstr_min_call_DP:
            common.WARNING("--gangstr-max-call-DP must be >= --gangstr-min-call-DP")
            return False
    if args.gangstr_min_call_Q is not None:
        if args.gangstr_min_call_Q < 0 or args.gangstr_min_call_Q > 1:
            common.WARNING("--gangstr-min-call-Q must be between 0 and 1")
            return False
        assert "Q" in invcf.formats
    if args.gangstr_expansion_prob_het is not None:
        if args.gangstr_expansion_prob_het < 0 or args.gangstr_expansion_prob_het > 1:
            common.WARNING("--gangstr-expansion-prob-het must be between 0 and 1")
            return False
        assert "QEXP" in invcf.formats
    if args.gangstr_expansion_prob_hom is not None:
        if args.gangstr_expansion_prob_hom < 0 or args.gangstr_expansion_prob_hom > 1:
            common.WARNING("--gangstr-expansion-prob-hom must be between 0 and 1")
            return False
        assert "QEXP" in invcf.formats
    if args.gangstr_expansion_prob_total is not None:
        if args.gangstr_expansion_prob_total < 0 or args.gangstr_expansion_prob_total > 1:
            common.WARNING("--gangstr-expansion-prob-total must be between 0 and 1")
            return False
        assert "QEXP" in invcf.formats
    if args.gangstr_require_support is not None:
        if args.gangstr_require_support < 0:
            common.WARNING("--gangstr-require-support must be >= 0")
            return False
        if args.gangstr_require_support > 0 and args.gangstr_readlen is None:
            common.WARNING("Using --gangstr-require-support requires setting --gangstr-readlen")
            return False
        if args.gangstr_readlen is not None and args.gangstr_readlen < 20:
            common.WARNING("--gangstr-readlen must be an integer value >= 20")
            return False
        assert "ENCLREADS" in invcf.formats and "FLNKREADS" in invcf.formats and "RC" in invcf.formats
    return True

def CheckAdVNTRFilters(invcf, args):
    r"""Check adVNTR call-level filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.advntr_min_call_DP is not None:
        if args.advntr_min_call_DP < 0:
            common.WARNING("--advntr-min-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.advntr_max_call_DP is not None:
        if args.advntr_max_call_DP < 0:
            common.WARNING("--advntr-max-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.advntr_min_call_DP is not None and args.advntr_max_call_DP is not None:
        if args.advntr_max_call_DP < args.advntr_min_call_DP:
            common.WARNING("--advntr-max-call-DP must be >= --advntr-min-call-DP")
            return False
    if args.advntr_min_spanning is not None:
        if args.advntr_min_spanning < 0:
            common.WARNING("--advntr-min-spanning must be >=0")
            return False
        assert "SR" in invcf.formats
    if args.advntr_min_flanking is not None:
        if args.advntr_min_flanking < 0:
            common.WARNING("--advntr-min-flanking must be >=0")
            return False
        assert "FR" in invcf.formats
    if args.advntr_min_ML is not None:
        if args.advntr_min_ML < 0:
            common.WARNING("--advntr-min-ML must be >= 0")
            return False
        assert "ML" in invcf.formats
    return True

# TODO remove pragma line when EH implemented
def CheckEHFilters(invcf, args): # pragma: no cover
    r"""Check ExpansionHunter call-level filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.eh_min_ADFL is not None:
        if args.eh_min_ADFL < 0:
            common.WARNING("--eh-min-ADFL must be >= 0")
            return False
        assert "ADFL" in invcf.formats
    if args.eh_min_ADIR is not None:
        if args.eh_min_ADIR < 0:
            common.WARNING("--eh-min-ADIR must be >= 0")
            return False
        assert "ADIR" in invcf.formats
    if args.eh_min_ADSP is not None:
        if args.eh_min_ADSP < 0:
            common.WARNING("--eh-min-ADSP must be >= 0")
            return False
        assert "ADSP" in invcf.formats
    if args.eh_min_call_LC is not None:
        if args.eh_min_call_LC < 0:
            common.WARNING("--eh-min-call-LC must be >= 0")
            return False
        assert "LC" in invcf.formats
    if args.eh_max_call_LC is not None:
        if args.eh_max_call_LC < 0:
            common.WARNING("--eh-max-call-LC must be >= 0")
            return False
        assert "LC" in invcf.formats
    if args.eh_min_call_LC is not None and args.eh_max_call_LC is not None:
        if args.eh_max_call_LC < args.eh_min_call_LC:
            common.WARNING("--eh-max-call-LC must be >= --eh-min-call-LC")
            return False
    return True

def CheckPopSTRFilters(invcf, args):
    r"""Check PopSTR call-level filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if args.popstr_min_call_DP is not None:
        if args.popstr_min_call_DP < 0:
            common.WARNING("--popstr-min-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.popstr_max_call_DP is not None:
        if args.popstr_max_call_DP < 0:
            common.WARNING("--popstr-max-call-DP must be >= 0")
            return False
        assert "DP" in invcf.formats
    if args.popstr_min_call_DP is not None and args.popstr_max_call_DP is not None:
        if args.popstr_max_call_DP < args.popstr_min_call_DP:
            common.WARNING("--popstr-max-call-DP must be >= --popstr-min-call-DP")
            return False
    if args.popstr_require_support is not None:
        if args.popstr_require_support < 0:
            common.WARNING("--popstr-require-support must be >= 0")
            return False
        assert "AD" in invcf.formats
    return True

def CheckFilters(invcf, args, vcftype):
    r"""Perform checks on user input for filters

    Parameters
    ----------
    invcf : str
        vcf.Reader object
    args : argparse namespace
        Contains user arguments
    vcftype : enum.
        Specifies which tool this VCF came from.
        Must be included in trh.VCFTYPES

    Returns
    -------
    checks : bool
        Set to True if all filters look ok.
        Set to False if filters are invalid
    """
    if not CheckLocusFilters(args, vcftype):
        return False

    # Check HipSTR specific filters
    if args.hipstr_max_call_flank_indel is not None or \
       args.hipstr_max_call_stutter is not None or \
       args.hipstr_min_supp_reads is not None or \
       args.hipstr_min_call_DP is not None or \
       args.hipstr_max_call_DP is not None or \
       args.hipstr_min_call_Q is not None:
        if vcftype != trh.VCFTYPES["hipstr"]:
            common.WARNING("HipSTR options can only be applied to HipSTR VCFs")
            return False
        else:
            if not CheckHipSTRFilters(invcf, args):
                return False

    # Check GangSTR specific filters
    if args.gangstr_min_call_DP is not None or \
       args.gangstr_max_call_DP is not None or \
       args.gangstr_min_call_Q is not None or \
       args.gangstr_expansion_prob_het is not None or \
       args.gangstr_expansion_prob_hom is not None or \
       args.gangstr_expansion_prob_total is not None or \
       args.gangstr_filter_span_only or \
       args.gangstr_filter_spanbound_only or \
       args.gangstr_filter_badCI or \
       args.gangstr_require_support is not None or \
       args.gangstr_readlen is not None:
        if vcftype != trh.VCFTYPES["gangstr"]:
            common.WARNING("GangSTR options can only be applied to GangSTR VCFs")
            return False
        else:
            if not CheckGangSTRFilters(invcf, args):
                return False

    # Check adVNTR specific filters
    if args.advntr_min_call_DP is not None or \
       args.advntr_max_call_DP is not None or \
       args.advntr_min_spanning is not None or \
       args.advntr_min_flanking is not None or \
       args.advntr_min_ML is not None:
        if vcftype != trh.VCFTYPES["advntr"]:
            common.WARNING("adVNTR options can only be applied to adVNTR VCFs")
            return False
        else:
            if not CheckAdVNTRFilters(invcf, args):
                return False

    # Check EH specific filters
    if args.eh_min_ADFL is not None or \
       args.eh_min_ADIR is not None or \
       args.eh_min_ADSP is not None or \
       args.eh_min_call_LC is not None or \
       args.eh_max_call_LC is not None:
        if vcftype != trh.VCFTYPES["eh"]:
            common.WARNING("ExpansionHunter options can only be applied to ExpansionHunter VCFs")
            return False
        else:  # pragma: no cover
            if not CheckEHFilters(invcf, args):  # pragma: no cover
                return False  # pragma: no cover

    # Check popSTR specific filters
    if args.popstr_min_call_DP is not None or \
       args.popstr_max_call_DP is not None or \
       args.popstr_require_support is not None:
        if vcftype != trh.VCFTYPES["popstr"]:
            common.WARNING("popSTR options can only be applied to popSTR VCFs")
            return False
        else:
            if not CheckPopSTRFilters(invcf, args):
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
    f = open(fname, "w")
    keys = list(loc_info.keys())
    assert "totalcalls" in keys and "PASS" in keys
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
    header = ["sample", "numcalls","meanDP"] + reasons
    f = open(fname, "w")
    f.write("\t".join(header)+"\n")
    for s in sample_info:
        assert "numcalls" in sample_info[s] and "totaldp" in sample_info[s]
        numcalls = sample_info[s]["numcalls"]
        if numcalls > 0:
            meancov = sample_info[s]["totaldp"]*1.0/numcalls
        else: meancov = 0
        items = [s, numcalls, meancov]
        for r in reasons: items.append(sample_info[s][r])
        f.write("\t".join([str(item) for item in items])+"\n")
    f.close()
    return True

def GetAllCallFilters(call_filters):
    r"""List all possible call filters

    Parameters
    ----------
    call_filters : list of filters.Reason
        List of all call-level filters

    Returns
    -------
    reasons : list of str
        A list of call-level filter reason strings
    """
    reasons = []
    for filt in call_filters:
        reasons.append(filt.name)
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
    # Add FILTER to end of formats
    if "FILTER" in record.FORMAT:
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
    else:
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':')+["FILTER"])
        record.FORMAT = record.FORMAT+":FILTER"
    for fmt in samp_fmt._fields:
        if fmt == "FILTER":
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
        if not sample.called:
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
            try:
                sample_info[sample.sample]["totaldp"] += sample["DP"]
            except AttributeError:
                sample_info[sample.sample]["totaldp"] += sample["LC"] # EH has LC field
            except AttributeError: pass # If no DP, or LC, no coverage data available
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "FILTER": sampdat.append("PASS")
                else: sampdat.append(sample[key])
        call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
        new_samples.append(call)
    record.samples = new_samples
    return record

def BuildCallFilters(args):
    r"""Build list of locus-level filters to include

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

    # HipSTR call-level filters
    if args.hipstr_max_call_flank_indel is not None:
        filter_list.append(filters.HipSTRCallFlankIndels(args.hipstr_max_call_flank_indel))
    if args.hipstr_max_call_stutter is not None:
        filter_list.append(filters.HipSTRCallStutter(args.hipstr_max_call_stutter))
    if args.hipstr_min_supp_reads is not None:
        filter_list.append(filters.HipSTRCallMinSuppReads(args.hipstr_min_supp_reads))
    if args.hipstr_min_call_DP is not None:
        filter_list.append(filters.CallFilterMinValue("HipSTRCallMinDepth", "DP", args.hipstr_min_call_DP))
    if args.hipstr_max_call_DP is not None:
        filter_list.append(filters.CallFilterMaxValue("HipSTRCallMaxDepth", "DP", args.hipstr_max_call_DP))
    if args.hipstr_min_call_Q is not None:
        filter_list.append(filters.CallFilterMinValue("HipSTRCallMinQ", "Q", args.hipstr_min_call_Q))

    # GangSTR call-level filters
    if args.gangstr_min_call_DP is not None:
        filter_list.append(filters.CallFilterMinValue("GangSTRCallMinDepth", "DP", args.gangstr_min_call_DP))
    if args.gangstr_max_call_DP is not None:
        filter_list.append(filters.CallFilterMaxValue("GangSTRCallMaxDepth", "DP", args.gangstr_max_call_DP))
    if args.gangstr_min_call_Q is not None:
        filter_list.append(filters.CallFilterMinValue("GangSTRCallMinQ", "Q", args.gangstr_min_call_Q))
    if args.gangstr_expansion_prob_het is not None:
        filter_list.append(filters.GangSTRCallExpansionProbHet(args.gangstr_expansion_prob_het))
    if args.gangstr_expansion_prob_hom is not None:
        filter_list.append(filters.GangSTRCallExpansionProbHom(args.gangstr_expansion_prob_hom))
    if args.gangstr_expansion_prob_total is not None:
        filter_list.append(filters.GangSTRCallExpansionProbTotal(args.gangstr_expansion_prob_total))
    if args.gangstr_filter_span_only:
        filter_list.append(filters.GangSTRCallSpanOnly())
    if args.gangstr_filter_spanbound_only:
        filter_list.append(filters.GangSTRCallSpanBoundOnly())
    if args.gangstr_filter_badCI:
        filter_list.append(filters.GangSTRCallBadCI())
    if args.gangstr_require_support is not None:
        filter_list.append(filters.GangSTRCallRequireSupport(args.gangstr_require_support, args.gangstr_readlen))

    # adVNTR call-level filters
    if args.advntr_min_call_DP is not None:
        filter_list.append(filters.CallFilterMinValue("AdVNTRCallMinDepth", "DP", args.advntr_min_call_DP))
    if args.advntr_max_call_DP is not None:
        filter_list.append(filters.CallFilterMaxValue("AdVNTRCallMaxDepth", "DP", args.advntr_max_call_DP))
    if args.advntr_min_spanning is not None:
        filter_list.append(filters.CallFilterMinValue("AdVNTRCallMinSpanning", "SR", args.advntr_min_spanning))
    if args.advntr_min_flanking is not None:
        filter_list.append(filters.CallFilterMinValue("AdVNTRCallMinFlanking", "FR", args.advntr_min_flanking))
    if args.advntr_min_ML is not None:
        filter_list.append(filters.CallFilterMinValue("AdVNTRCallMinML", "ML", args.advntr_min_ML))

    # EH call-level filters
    if args.eh_min_call_LC is not None:
        filter_list.append(filters.CallFilterMinValue("EHCallMinDepth", "LC", args.eh_min_call_LC))  # pragma: no cover
    if args.eh_max_call_LC is not None:
        filter_list.append(filters.CallFilterMaxValue("EHCallMaxDepth", "LC", args.eh_max_call_LC))  # pragma: no cover
    if args.eh_min_ADFL is not None:
        filter_list.append(filters.CallFilterMinValue("EHCallMinADFL", "ADFL", args.eh_min_ADFL))  # pragma: no cover
    if args.eh_min_ADIR is not None:
        filter_list.append(filters.CallFilterMinValue("EHCallMinADFL", "ADIR", args.eh_min_ADIR))  # pragma: no cover
    if args.eh_min_ADSP is not None:
        filter_list.append(filters.CallFilterMinValue("EHCallMinADSP", "ADSP", args.eh_min_ADSP)) # pragma: no cover

    # popSTR call-level filters
    if args.popstr_min_call_DP is not None:
        filter_list.append(filters.CallFilterMinValue("PopSTRMinCallDepth", "DP", args.popstr_min_call_DP))
    if args.popstr_max_call_DP is not None:
        filter_list.append(filters.CallFilterMaxValue("PopSTRMaxCallDepth", "DP", args.popstr_max_call_DP))
    if args.popstr_require_support is not None:
        filter_list.append(filters.PopSTRCallRequireSupport(args.popstr_require_support))
    return filter_list

def BuildLocusFilters(args, vcftype: trh.VCFTYPES):
    r"""Build list of locus-level filters to include.

    These filters should in general not be tool specific

    Parameters
    ---------
    args : argparse namespace
       User input arguments used to decide on filters
    vcftype:
        the type of the vcf we're working with

    Returns
    -------
    filter_list : list of filters.Filter
       List of locus-level filters
    """
    filter_list = []
    if args.min_locus_callrate is not None:
        filter_list.append(filters.Filter_MinLocusCallrate(args.min_locus_callrate))
    if args.min_locus_hwep is not None:
        filter_list.append(filters.Filter_MinLocusHWEP(args.min_locus_hwep,
                                                       vcftype, args.use_length))
    if args.min_locus_het is not None:
        filter_list.append(filters.Filter_MinLocusHet(args.min_locus_het,
                                                      vcftype, args.use_length))
    if args.max_locus_het is not None:
        filter_list.append(filters.Filter_MaxLocusHet(args.max_locus_het,
                                                      vcftype, args.use_length))
    if args.filter_hrun:
        filter_list.append(filters.Filter_LocusHrun())
    if args.filter_regions is not None:
        filter_region_files = args.filter_regions.split(",")
        if args.filter_regions_names is not None:
            filter_region_names = args.filter_regions_names.split(",")
        else: filter_region_names = [str(item) for item in list(range(len(filter_region_files)))]
        for i in range(len(filter_region_names)):
            region_filter = filters.create_region_filter(filter_region_names[i], filter_region_files[i])
            if region_filter is not None:
                filter_list.append(region_filter)
            else:
                raise ValueError('Could not load regions file: {}'.format(filter_region_files[i]))
    return filter_list

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(__doc__)
    # In/out are always relevant
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)
    inout_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VCFTYPES.__members__], type=str, default="auto")

    # Locus-level filters are not specific to any tool
    locus_group = parser.add_argument_group("Locus-level filters (tool agnostic)")
    locus_group.add_argument("--min-locus-callrate", help="Minimum locus call rate", type=float)
    locus_group.add_argument("--min-locus-hwep", help="Filter loci failing HWE at this p-value threshold", type=float)
    locus_group.add_argument("--min-locus-het", help="Minimum locus heterozygosity", type=float)
    locus_group.add_argument("--max-locus-het", help="Maximum locus heterozygosity", type=float)
    locus_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    locus_group.add_argument("--filter-regions", help="Comma-separated list of BED files of regions to filter. Must be bgzipped and tabix indexed", type=str)
    locus_group.add_argument("--filter-regions-names", help="Comma-separated list of filter names for each BED filter file", type=str)
    locus_group.add_argument("--filter-hrun", help="Filter STRs with long homopolymer runs.", action="store_true")
    locus_group.add_argument("--drop-filtered", help="Drop filtered records from output", action="store_true")

    ###### Tool specific filters #####
    hipstr_call_group = parser.add_argument_group("Call-level filters specific to HipSTR output")
    hipstr_call_group.add_argument("--hipstr-max-call-flank-indel", help="Maximum call flank indel rate", type=float)
    hipstr_call_group.add_argument("--hipstr-max-call-stutter", help="Maximum call stutter rate", type=float)
    hipstr_call_group.add_argument("--hipstr-min-supp-reads", help="Minimum supporting reads for each allele", type=int)
    hipstr_call_group.add_argument("--hipstr-min-call-DP", help="Minimum call coverage", type=int)
    hipstr_call_group.add_argument("--hipstr-max-call-DP", help="Maximum call coverage", type=int)
    hipstr_call_group.add_argument("--hipstr-min-call-Q", help="Minimum call quality score", type=float)

    gangstr_call_group = parser.add_argument_group("Call-level filters specific to GangSTR output")
    gangstr_call_group.add_argument("--gangstr-min-call-DP", help="Minimum call coverage", type=int)
    gangstr_call_group.add_argument("--gangstr-max-call-DP", help="Maximum call coverage", type=int)
    gangstr_call_group.add_argument("--gangstr-min-call-Q", help="Minimum call quality score", type=float)
    gangstr_call_group.add_argument("--gangstr-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    gangstr_call_group.add_argument("--gangstr-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    gangstr_call_group.add_argument("--gangstr-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)
    gangstr_call_group.add_argument("--gangstr-filter-span-only", help="Filter out all calls that only have spanning read support", action="store_true")
    gangstr_call_group.add_argument("--gangstr-filter-spanbound-only", help="Filter out all reads except spanning and bounding", action="store_true")
    gangstr_call_group.add_argument("--gangstr-filter-badCI", help="Filter regions where the ML estimate is not in the CI", action="store_true")
    gangstr_call_group.add_argument("--gangstr-require-support", help="Require each allele call to have at least n supporting reads", type=int)
    gangstr_call_group.add_argument("--gangstr-readlen", help="Read length used (bp). Required if using --require-support", type=int)

    advntr_call_group = parser.add_argument_group("Call-level filters specific to adVNTR output")
    advntr_call_group.add_argument("--advntr-min-call-DP", help="Minimum call coverage", type=int)
    advntr_call_group.add_argument("--advntr-max-call-DP", help="Maximum call coverage", type=int)
    advntr_call_group.add_argument("--advntr-min-spanning", help="Minimum spanning read count (SR field)", type=int)
    advntr_call_group.add_argument("--advntr-min-flanking", help="Minimum flanking read count (FR field)", type=int)
    advntr_call_group.add_argument("--advntr-min-ML", help="Minimum value of maximum likelihood (ML field)", type=float)

    eh_call_group = parser.add_argument_group("Call-level filters specific to ExpansionHunter output")
    eh_call_group.add_argument("--eh-min-ADFL", help="Minimum number of flanking reads consistent with the allele", type=int)
    eh_call_group.add_argument("--eh-min-ADIR", help="Minimum number of in-repeat reads consistent with the allele", type=int)
    eh_call_group.add_argument("--eh-min-ADSP", help="Minimum number of spanning reads consistent with the allele", type=int)
    eh_call_group.add_argument("--eh-min-call-LC", help="Minimum call coverage", type=int)
    eh_call_group.add_argument("--eh-max-call-LC", help="Maximum call coverage", type=int)
    # TODO: add SO field filter. After clarifying possible formats it can take

    popstr_call_group = parser.add_argument_group("Call-level filters specific to PopSTR output")
    popstr_call_group.add_argument("--popstr-min-call-DP", help="Minimum call coverage", type=int)
    popstr_call_group.add_argument("--popstr-max-call-DP", help="Maximum call coverage", type=int)
    popstr_call_group.add_argument("--popstr-require-support", help="Require each allele call to have at least n supporting reads", type=int)

    # Debugging options
    debug_group = parser.add_argument_group("Debugging parameters")
    debug_group.add_argument("--num-records", help="Only process this many records", type=int)
    debug_group.add_argument("--die-on-warning", help="Quit if a record can't be parsed", action="store_true")
    debug_group.add_argument("--verbose", help="Print out extra info", action="store_true")
    # Version option
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    # Load VCF file
    if not os.path.exists(args.vcf):
        common.WARNING("%s does not exist"%args.vcf)
        return 1
    invcf = vcf.Reader(filename=args.vcf)

    # Set up record harmonizer and infer VCF type
    vcftype = trh.InferVCFType(invcf)

    # Check filters all make sense
    if not CheckFilters(invcf, args, vcftype): return 1

    # Set up locus-level filter list
    try:
        filter_list = BuildLocusFilters(args, vcftype)
    except ValueError:
        return 1
    invcf.filters = {}
    for f in filter_list:
        short_doc = f.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        invcf.filters[f.filter_name()] = _Filter(f.filter_name(), short_doc)

    # Set up call-level filters
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
    if not os.path.exists(os.path.dirname(os.path.abspath(args.out))):
        common.WARNING("Output directory does not exist")
        return 1
    outvcf = MakeWriter(args.out + ".vcf", invcf, " ".join(sys.argv))
    if outvcf is None: return 1

    # Set up sample info
    all_reasons = GetAllCallFilters(call_filters)
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
            trrecord = trh.HarmonizeRecord(vcftype, record)
            # Recalculate locus-level INFO fields
            record.INFO["HRUN"] = utils.GetHomopolymerRun(record.REF)
            if record.num_called > 0:
                allele_freqs = trrecord.GetAlleleFreqs(uselength=args.use_length)
                genotype_counts = trrecord.GetGenotypeCounts(uselength=args.use_length)
                record.INFO["HET"] = utils.GetHeterozygosity(allele_freqs)
                record.INFO["HWEP"] = utils.GetHardyWeinbergBinomialTest(allele_freqs, genotype_counts)
                record.INFO["AC"] = [int(item*(3*record.num_called)) for item in record.aaf]
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

def run(): # pragma: no cover
    args = getargs()
    retcode = main(args)
    sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
