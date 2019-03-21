#!/usr/bin/env python3

"""
Tool for post-processing and priotization of STR genotypes

"""

# Imports


# Load local libraries
import strtools.utils.common as common
import dumpSTR.filters as filters

# Load external libraries
import argparse
import inspect
import vcf
import sys

def MakeWriter(outfile, invcf, command):
    invcf.metadata["command-postmaSTR"] = [command]
    writer = vcf.Writer(open(outfile, "w"), invcf)
    return writer

def CheckFilters(args):
    """
    Perform checks on user input for filters
    Input:
    - invcf (vcf.Reader)
    - args (argparse namespace)
    Exit program if checks fail
    """
    if args.affec_max_expansion_prob_het is not None:
        if args.affec_max_expansion_prob_het < 0 or args.affec_max_expansion_prob_het > 1:
            common.ERROR("--affec-max-expansion-prob-het must be between 0 and 1")
    if args.affec_min_expansion_prob_het is not None:
        if args.affec_min_expansion_prob_het < 0 or args.affec_min_expansion_prob_het > 1:
            common.ERROR("--affec-min-expansion-prob-het must be between 0 and 1")
    if args.affec_min_expansion_prob_het is not None and args.affec_max_expansion_prob_het is not None:
        if args.affec_min_expansion_prob_het > args.affec_max_expansion_prob_het:
            common.ERROR("--affec-min-expansion-prob-het must be less than --affec-max-expansion-prob-het")
    if args.unaff_max_expansion_prob_het is not None:
        if args.unaff_max_expansion_prob_het < 0 or args.unaff_max_expansion_prob_het > 1:
            common.ERROR("--unaff-max-expansion-prob-het must be between 0 and 1")
    if args.unaff_min_expansion_prob_het is not None:
        if args.unaff_min_expansion_prob_het < 0 or args.unaff_min_expansion_prob_het > 1:
            common.ERROR("--unaff-min-expansion-prob-het must be between 0 and 1")
    if args.unaff_min_expansion_prob_het is not None and args.unaff_max_expansion_prob_het is not None:
        if args.unaff_min_expansion_prob_het > args.unaff_max_expansion_prob_het:
            common.ERROR("--unaff-min-expansion-prob-het must be less than --unaff-max-expansion-prob-het")
    if args.affec_max_expansion_prob_hom is not None:
        if args.affec_max_expansion_prob_hom < 0 or args.affec_max_expansion_prob_hom > 1:
            common.ERROR("--affec-max-expansion-prob-hom must be between 0 and 1")
    if args.affec_min_expansion_prob_hom is not None:
        if args.affec_min_expansion_prob_hom < 0 or args.affec_min_expansion_prob_hom > 1:
            common.ERROR("--affec-min-expansion-prob-hom must be between 0 and 1")
    if args.affec_min_expansion_prob_hom is not None and args.affec_max_expansion_prob_hom is not None:
        if args.affec_min_expansion_prob_hom < args.affec_max_expansion_prob_hom:
            common.ERROR("--affec-min-expansion-prob-hom must be less than --affec-max-expansion-prob-hom")
    if args.unaff_max_expansion_prob_hom is not None:
        if args.unaff_max_expansion_prob_hom < 0 or args.unaff_max_expansion_prob_hom > 1:
            common.ERROR("--unaff-max-expansion-prob-hom must be between 0 and 1")
    if args.unaff_min_expansion_prob_hom is not None:
        if args.unaff_min_expansion_prob_hom < 0 or args.unaff_min_expansion_prob_hom > 1:
            common.ERROR("--unaff-min-expansion-prob-hom must be between 0 and 1")
    if args.unaff_min_expansion_prob_hom is not None and args.unaff_max_expansion_prob_hom is not None:
        if args.unaff_min_expansion_prob_hom < args.unaff_max_expansion_prob_hom:
            common.ERROR("--unaff-min-expansion-prob-hom must be less than --unaff-max-expansion-prob-hom")
    if args.affec_max_expansion_prob_total is not None:
        if args.affec_max_expansion_prob_total < 0 or args.affec_max_expansion_prob_total > 1:
            common.ERROR("--affec-max-expansion-prob-total must be between 0 and 1")
    if args.affec_min_expansion_prob_total is not None:
        if args.affec_min_expansion_prob_total < 0 or args.affec_min_expansion_prob_total > 1:
            common.ERROR("--affec-min-expansion-prob-total must be between 0 and 1")
    if args.affec_min_expansion_prob_total is not None and args.affec_max_expansion_prob_total is not None:
        if args.affec_min_expansion_prob_total < args.affec_max_expansion_prob_total:
            common.ERROR("--affec-min-expansion-prob-total must be less than --affec-max-expansion-prob-total")
    if args.unaff_max_expansion_prob_total is not None:
        if args.unaff_max_expansion_prob_total < 0 or args.unaff_max_expansion_prob_total > 1:
            common.ERROR("--unaff-max-expansion-prob-total must be between 0 and 1")
    if args.unaff_min_expansion_prob_total is not None:
        if args.unaff_min_expansion_prob_total < 0 or args.unaff_min_expansion_prob_total > 1:
            common.ERROR("--unaff-min-expansion-prob-total must be between 0 and 1")
    if args.unaff_min_expansion_prob_total is not None and args.unaff_max_expansion_prob_total is not None:
        if args.unaff_min_expansion_prob_total < args.unaff_max_expansion_prob_total:
            common.ERROR("--unaff-min-expansion-prob-total must be less than --unaff-max-expansion-prob-total")
    if args.affec_min_call_count < 0:
        common.ERROR("Minimum number of affected calls (" + str(args.affec_min_call_count) + \
                     ") must be 0 or more")
    if args.unaff_min_call_count < 0:
        common.ERROR("Minimum number of unaffected calls (" + str(args.unaff_min_call_count) + \
                     ") must be 0 or more")

def GetAllCallFilters():
    """
    List all possible call filters by
    extracting from filters module
    Output:
    - reasons (list<str>): list of call-level filter reasons
    """
    reasons = []
    for name, obj in inspect.getmembers(filters):
        if inspect.isclass(obj) and issubclass(obj, filters.Reason) and not obj is filters.Reason:
            reasons.append(obj.name)
    return reasons

def FilterCall(sample, call_filters):
    """
    Apply call-level filters and return filter reason.
    Input:
    - sample (vcf._Call)
    - call_filters (list<filters.Reason>): list of call level filters
    Return:
    - reason (list<str>): list of string description of filter reasons
    """
    reasons = []
    for cfilt in call_filters:
        if cfilt(sample) is not None: reasons.append(cfilt.GetReason())
    return reasons
def PrintAcceptedRecord(record, isAffected):
    affec_CN = []
    unaff_CN = []
    for sample in record:
        if isAffected[sample.sample]:
            affec_CN.append(str(sample['REPCN']))
        else:
            unaff_CN.append(str(sample['REPCN']))
    print("\t".join([str(record.CHROM), str(record.POS), record.INFO['RU'], str(record.INFO['EXPTHRESH'])])\
          + "\t" + ", ".join(unaff_CN) + "\t" +  ", ".join(affec_CN))

def ApplyPostMaSTRCallFilters(record, reader, call_filters, sample_info, isAffected, args):
    """
    Apply call level filters to a record
    Return a new record with FILTER populated for each sample
    Update sample_info with sample level stats
    Input:
    - record (vcf._Record)
    - reader (vcf.Reader)
    - call_filters ({'affec': list<filters.Reason>, 'unaff': list<filters.Reason>}): list of call filters for affected and unaffected
    - sample_info (dict): dictionary of sample stats
    - isAffected (dict): dictionary of phenotypes
    Output:
    - modified record (vcf._Record). 
    """
    min_affec = args.affec_min_call_count
    min_unaff = args.unaff_min_call_count
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
    num_affec_calls = 0
    num_unaff_calls = 0
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
        
        if isAffected[sample.sample]:
            filter_reasons = FilterCall(sample, call_filters['affec'])
        else:
            filter_reasons = FilterCall(sample, call_filters['unaff'])

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
            if isAffected[sample.sample]:
                num_affec_calls = num_affec_calls + 1
            else:
                num_unaff_calls = num_unaff_calls + 1
            sample_info[sample.sample]["numcalls"] += 1
            sample_info[sample.sample]["totaldp"] += sample["DP"]
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "FILTER": sampdat.append("PASS")
                else: sampdat.append(sample[key])
        call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
        new_samples.append(call)
    record.samples = new_samples
    if num_unaff_calls >= min_unaff and num_affec_calls >= min_affec:
        PrintAcceptedRecord(record, isAffected)
        return record
    else:
        return None


def BuildPostMaSTRCallFilters(args):
    """
    Build list of locus-level filters to include for affected and unaffected samples.
    Input:
    - args (namespace from parser.parse_args)
    
    Output:
    - cdict ({'affec': list<filters.Filter>, 'unaff': list<filters.Filter>}): list of call-level filters
    """
    cdict = {'affec':[], 'unaff':[]}
    if args.affec_min_expansion_prob_het is not None:
        cdict['affec'].append(filters.ProbHet(args.affec_min_expansion_prob_het))
    if args.affec_max_expansion_prob_het is not None:
        cdict['affec'].append(filters.MaxProbHet(args.affec_max_expansion_prob_het))
    if args.affec_min_expansion_prob_hom is not None:
        cdict['affec'].append(filters.ProbHom(args.affec_min_expansion_prob_hom))
    if args.affec_max_expansion_prob_hom is not None:
        cdict['affec'].append(filters.MaxProbHom(args.affec_max_expansion_prob_hom))
    if args.affec_min_expansion_prob_total is not None:
        cdict['affec'].append(filters.ProbTotal(args.affec_min_expansion_prob_total))
    if args.affec_max_expansion_prob_total is not None:
        cdict['affec'].append(filters.MaxProbTotal(args.affec_max_expansion_prob_total))
    if args.unaff_min_expansion_prob_het is not None:
        cdict['unaff'].append(filters.ProbHet(args.unaff_min_expansion_prob_het))
    if args.unaff_max_expansion_prob_het is not None:
        cdict['unaff'].append(filters.MaxProbHet(args.unaff_max_expansion_prob_het))
    if args.unaff_min_expansion_prob_hom is not None:
        cdict['unaff'].append(filters.ProbHom(args.unaff_min_expansion_prob_hom))
    if args.unaff_max_expansion_prob_hom is not None:
        cdict['unaff'].append(filters.MaxProbHom(args.unaff_max_expansion_prob_hom))
    if args.unaff_min_expansion_prob_total is not None:
        cdict['unaff'].append(filters.ProbTotal(args.unaff_min_expansion_prob_total))
    if args.unaff_max_expansion_prob_total is not None:
        cdict['unaff'].append(filters.MaxProbTotal(args.unaff_max_expansion_prob_total))
    return cdict


def ParseFam(args):
    """
    Parse fam file and extract affected and unaffected sample IDs.
    Input:
    - args (namespace from parser.parse_args)

    Output:
    - isAffected ({str: bool}): dictionary for affected and unaffected sample status
    """
    filename = args.fam
    min_affec = args.affec_min_call_count
    min_unaff = args.unaff_min_call_count
    isAffected = {}
    with open(filename, 'r') as f:
        i = 0
        count_affec = 0
        count_unaff = 0
        for line in f:
            i = i + 1
            recs = line.strip().split('\t')
            if len(recs) < 6:
                common.ERROR("Insufficient number of columns in line " + str(i) + " of fam file: " + filename)
            sid = recs[1]
            phe = recs[5]
            if phe == '2':
                isAffected[sid] = True
                count_affec = count_affec + 1
            else:
                isAffected[sid] = False
                count_unaff = count_unaff + 1
        if count_affec < min_affec:
            common.ERROR("Minimum number of affected calls (" + str(min_affec) + \
                         ") larger than number of affected samples in fam file (" + str(count_affec) + ")")
        if count_unaff < min_unaff:
            common.ERROR("Minimum number of unaffected calls (" + str(min_unaff) + \
                         ") larger than number of unaffected samples in fam file (" + str(count_unaff) + ")")
    return isAffected
                

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)    
    inout_group.add_argument("--fam", help="Input fam file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)

    locus_group = parser.add_argument_group("Locus-level filters")
    locus_group.add_argument("--use-length", help="Calculate per-locus stats (het, HWE) collapsing alleles by length", action="store_true")
    
    #### Affected sample filters
    affec_group = parser.add_argument_group("Call-level filters specific to affected samples")
    affec_group.add_argument("--affec-max-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    affec_group.add_argument("--affec-min-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion more than this", type=float)
    affec_group.add_argument("--affec-max-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    affec_group.add_argument("--affec-min-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion more than this", type=float)
    affec_group.add_argument("--affec-max-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)
    affec_group.add_argument("--affec-min-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion more than this", type=float)
    affec_group.add_argument("--affec-min-call-count", help="Minimum number of affected calls per TR region", type=int, default=1)
    
    #### Unaffected sample filters
    unaff_group = parser.add_argument_group("Call-level filters specific to unaffected samples")
    unaff_group.add_argument("--unaff-max-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion less than this", type=float)
    unaff_group.add_argument("--unaff-min-expansion-prob-het", help="Expansion prob-value threshold. Filters calls with probability of heterozygous expansion more than this", type=float)
    unaff_group.add_argument("--unaff-max-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion less than this", type=float)
    unaff_group.add_argument("--unaff-min-expansion-prob-hom", help="Expansion prob-value threshold. Filters calls with probability of homozygous expansion more than this", type=float)
    unaff_group.add_argument("--unaff-max-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion less than this", type=float)
    unaff_group.add_argument("--unaff-min-expansion-prob-total", help="Expansion prob-value threshold. Filters calls with probability of total expansion more than this", type=float)
    unaff_group.add_argument("--unaff-min-call-count", help="Minimum number of unaffected calls per TR region", type=int, default=1)

    debug_group = parser.add_argument_group("Debugging parameters")
    debug_group.add_argument("--num-records", help="Only process this many records", type=int)
    debug_group.add_argument("--die-on-warning", help="Quit if a record can't be parsed", action="store_true")
    debug_group.add_argument("--verbose", help="Print out extra info", action="store_true")

    args = parser.parse_args()
    # Load VCF file
    invcf = vcf.Reader(filename=args.vcf)
    # Load FAM file
    isAffected = ParseFam(args)

    # Checking filters for validity
    CheckFilters(args)
    invcf.filters = {}
    
    # Two sets of parameters set for call filters (for affected and unaffected)
    call_filters = BuildPostMaSTRCallFilters(args)

    # Set up output files
    outvcf = MakeWriter(args.out + ".vcf", invcf, " ".join(sys.argv))

    all_reasons = GetAllCallFilters()
    sample_info = {}
    for s in invcf.samples:
        sample_info[s] = {"numcalls": 0, "totaldp": 0}
        for r in all_reasons: sample_info[s][r]  = 0

    print("CHROM\tPOS\tMOTIF\tEXPTHRESH\tUNAFF\tAFFEC")
    # Go through each record
    record_counter = 0
    while True:
        try:
            record = next(invcf)
        except IndexError:
            common.WARNING("Skipping TR that couldn't be parsed by PyVCF. Check VCF format")
            if args.die_on_warning: sys.exit(1)
        except StopIteration: break
        if args.verbose:
            common.MSG("Processing %s:%s"%(record.CHROM, record.POS))
        record_counter += 1
        if args.num_records is not None and record_counter > args.num_records: break
        # Call-level filters
        record = ApplyPostMaSTRCallFilters(record, invcf, call_filters, sample_info,\
                                           isAffected, args)
        output_record = True
        if record is None:
            output_record = False
        if output_record:
            outvcf.write_record(record)

if __name__ == "__main__":
    main()
