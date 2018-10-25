#!/usr/bin/env python

"""
Tool for filtering and QC of STR genotypes

Example command:

./dumpSTR.py \
--vcf /storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr22.allfilters.vcf.gz \
--out test \
--min-call-DP 10 \
--max-call-DP 1000 \
--min-call-Q 0.9 \
--max-call-flank-indel 0.15 \
--max-call-stutter 0.15 \
--min-locus-callrate 0.8 \
--min-locus-hwep 0.01 \
--min-locus-het 0 \
--max-locus-het 1 \
--use-length \
--filter-regions /storage/resources/dbase/human/hg19/hg19_segmentalduplications.bed \
--filter-regions-names SEGDUP \
--filter-hrun
"""

# TODO:
# - Implement locus-level fitlers in filters.py
# - Add new INFO fields with updated info (allele counts
# - Info for log files
# - better checking of user input

# Load external libraries
import argparse
import vcf
from vcf.parser import _Filter
from vcf.parser import _Format

# Load custom libraries
import filters

def FilterCall(sample, args):
    if args.min_call_DP is not None:
        try:
            if sample["DP"] < args.min_call_DP: return "LowCallDepth"
        except KeyError:
            sys.stderr.write("ERROR: No DP field found\n")
            sys.exit(1)
    if args.max_call_DP is not None:
        try:
            if sample["DP"] > args.max_call_DP: return "HighCallDepth"
        except KeyError:
            sys.stderr.write("ERROR: No DP field found\n")
            sys.exit(1)
    if args.min_call_Q is not None:
        try:
            if sample["Q"] < args.min_call_Q: return "LowCallQ"
        except KeyError:
            sys.stderr.write("ERROR: No Q field found\n")
            sys.exit(1)
    if args.max_call_flank_indel is not None:
        try:
            if 1.0*sample['DFLANKINDEL']/sample['DP'] > args.max_call_flank_indel:
                return "CallFlankIndels"
        except KeyError:
            sys.stderr.write("ERROR: No DP or DFLANKINDEL field found\n")
            sys.exit(1)
    if args.max_call_stutter is not None:
        try:
            if 1.0*sample['DSTUTTER']/sample['DP'] > args.max_call_stutter:
                return "CallStutter"
        except KeyError:
            sys.stderr.write("ERROR: No DP or DSTUTTER field found\n")
            sys.exit(1)
    return None

def ApplyCallFilters(record, reader, args):
    if "FILTER" in record.FORMAT:
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
    else: samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':')+["FILTER"])
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
        filter_reason = FilterCall(sample, args)
        if filter_reason is not None:
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "GT":
                    sampdat.append("./.")
                else:
                    if key == "FILTER": sampdat.append(filter_reason.replace(" ", "_").upper())
                    else: sampdat.append(None)
        else:
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key == "FILTER": sampdat.append("PASS")
                else: sampdat.append(sample[key])
        call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
        new_samples.append(call)
    record.samples = new_samples
    return record

def BuildFilters(args):
    fdict = []
    if args.min_locus_callrate is not None:
        fdict.append(filters.Filter_MinLocusCallrate(args.min_locus_callrate))
    if args.min_locus_hwep is not None:
        fdict.append(filters.Filter_MinLocusHWEP(args.min_locus_hwep, args.use_length))
    if args.min_locus_het is not None:
        fdict.append(filters.Filter_MinLocusHet(args.min_locus_het, args.use_length))
    if args.max_locus_het is not None:
        fdict.append(filters.Filter_MaxLocusHet(args.max_locus_het, args.use_length))
    if args.filter_hrun is not None:
        fdict.append(filters.Filter_LocusHrun())
    if args.filter_regions is not None:
        filter_region_files = args.filter_regions.split(",")
        if args.filter_regions_names is not None:
            filter_region_names = args.filter_regions_names.split(",")
            if len(filter_region_names) != len(filter_region_files):
                sys.stderr.write("ERROR: length of --filter-regions-names must match --filter-regions\n")
                sys.exit(1)
        else: filter_region_names = [str(item) for item in list(range(len(filter_region_files)))]
        for i in range(len(filter_region_names)):
            fdict.append(filters.create_region_filter(filter_region_names[i], filter_region_files[i]))
    return fdict

def main():
    parser = argparse.ArgumentParser(__doc__)
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file", type=str, required=True)
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)

    call_group = parser.add_argument_group("Call-level filters") # TODO add GangSTR filters
    call_group.add_argument("--min-call-DP", help="Minimum call coverage", type=int)
    call_group.add_argument("--max-call-DP", help="Maximum call coverage", type=int)
    call_group.add_argument("--min-call-Q", help="Minimum call quality score", type=float)
    call_group.add_argument("--max-call-flank-indel", help="Maximum call flank indel rate", type=float)
    call_group.add_argument("--max-call-stutter", help="Maximum call stutter rate", type=float)

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

    args = parser.parse_args()

    # Load VCF file
    invcf = vcf.Reader(open(args.vcf, "rb"))

    # Set up filter list
    invcf.filters = {}
    filter_list = BuildFilters(args)
    for f in filter_list:
        short_doc = f.__doc__ or ''
        short_doc = short_doc.split('\n')[0].lstrip()
        invcf.filters[f.filter_name()] = _Filter(f.filter_name(), short_doc)

    # Add new FORMAT fields
    if "FILTER" not in invcf.formats:
        invcf.formats["FILTER"] = _Format("FILTER", 1, "String", "Call-level filter")

    # Add new INFO fields - TODO

    # Set up output files
    outvcf = vcf.Writer(open(args.out + ".vcf", "w"), invcf)

    # Open again to get rid of any changes
    invcf = vcf.Reader(open(args.vcf, "rb"))
    loclog = open(args.out + ".loclog.tab", "w")
    samplog = open(args.out + ".samplog.tab", "w")

    # Go through each record
    for record in invcf:
        # Call-level filters
        record = ApplyCallFilters(record, invcf, args)

        # Locus-level filters
        record.FILTER = None
        output_record = True
        for filt in filter_list:
            if filt(record) == None: continue
            if args.drop_filtered:
                output_record = False
                break
            record.add_filter(filt.filter_name())
        if output_record:
            if record.FILTER is None and not args.drop_filtered: record.FILTER = "PASS"
            outvcf.write_record(record)

if __name__ == "__main__":
    main()
