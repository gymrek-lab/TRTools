"""
Tool for annotating TR VCF files

TODO:

* Add documentation to functions
* Add unit tests
* Add command line tests
* Add README and link to other docs
* Add tensorqtl dosage output
* force normalize for pgen output
* Update docs for CALLERS about Beagle
* test we give error if no annotation options specified

"""

import argparse
import cyvcf2
import enum
import numpy as np
import os
from pgenlib import PgenWriter
import sys

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils
from trtools import __version__

# Variables for PGEN output
DEFAULT_PGEN_BATCHSIZE = 1000
DUMMY_REF = "A"
DUMMY_ALT = "T"

# Info fields copied from reference panel for each tool
INFOFIELDS = {
    trh.VcfTypes.hipstr: ["START","END","PERIOD"],
    trh.VcfTypes.advntr: ["RU", "VID"],
    trh.VcfTypes.gangstr: ["RU"],
    trh.VcfTypes.eh: ["RU","VARID","RL"]
}

class OutputFileTypes(enum.Enum):
    """Different supported output file types."""
    vcf = "vcf"
    pgen = "pgen"
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)

def GetVCFWriter(reader, fname, command, vcftype, dosage_type=None, refreader=None):
    reader.add_to_header("##command-AnnotaTR=" + command)
    # Add dosage lines to header
    if dosage_type is not None:
        reader.add_format_to_header({
            'ID': 'TRDS',
            'Description': 'TR genotype dosage, method={method}'.format(method=str(dosage_type)),
            'Type': 'Float',
            'Number': 1
        })
        reader.add_info_to_header({
            'ID': 'DSLEN',
            'Description': 'Minimum and maximum dosages, used if normalization was applied',
            'Type': 'Float',
            'Number': '2'
        })
    # Copy over relevant header lines from refpanel
    if refreader is not None:
        refheader = refreader.raw_header.split("\n")
        for item in refheader:
            if item.startswith("##source"):
                reader.add_to_header("##preimputation_source" + item.strip()[8:])
            if item.startswith("##command"):
                reader.add_to_header("##preimputation_command" + item.strip()[9:])
            if item.startswith("##contig") or item.startswith("##ALT"):
                reader.add_to_header(item.strip())
        for infofield in INFOFIELDS[vcftype]:
            if refreader.contains(infofield):
                # Note can't pass headerinfo directly because extra
                # quotes in Description need to be removed...
                headerinfo = refreader.get_header_type(infofield)
                reader.add_info_to_header({
                    'ID': headerinfo["ID"],
                    'Description': headerinfo["Description"].replace('"',''),
                    'Type': headerinfo["Type"],
                    "Number": headerinfo["Number"]
                    })
            else:
                common.WARNING("Could not find required header field {field} in refpanel".format(field=infofield))
                return None
    writer = cyvcf2.Writer(fname, reader)
    return writer

def LoadMetadataFromRefPanel(refreader, vcftype):
    metadata = {} # chr:pos:ref->info
    for record in refreader:
        locdata = {}
        for infofield in INFOFIELDS[vcftype]:
            infodata = record.INFO.get(infofield, None)
            if infodata is not None:
                locdata[infofield] = infodata
        # Check if we have all required INFO fields
        # Otherwise assume it is not a TR
        if len(locdata.keys())!=len(INFOFIELDS[vcftype]):
            continue
        # Add to metadata
        locuskey = "{chrom}:{pos}:{ref}".format(
            chrom=record.CHROM,
            pos=record.POS, 
            ref=record.REF
        )
        # Quit if we found a duplicate TR locus
        if locuskey in metadata.keys():
            common.WARNING("Error: duplicate locus detected in refpanel: {locus}".format(locus=locuskey))
            return {}
        metadata[locuskey] = locdata
    return metadata

def GetPGenPvarWriter(reader, outprefix, variant_ct):
    # Write .psam 
    with open(outprefix+".psam", "w") as f:
        f.write("#IID\tSEX\n")
        for sample in reader.samples:
            f.write("{sample}\t0\n")
    # Get pvar writer
    pvar_writer = open(outprefix+".pvar", "w")
    pvar_writer.write("##INFO=<ID=DSLEN,Number=2,Type=Float,Description=\"Minimum and maximum dosages, used if normalization was applied\">\n")
    pvar_writer.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "INFO"])+"\n")
    # Get pgen writer we will continue writing to
    pgen_writer = PgenWriter(bytes(outprefix+".pgen", "utf8"),
        len(reader.samples), variant_ct=variant_ct, dosage_present=True)
    return pgen_writer, pvar_writer

def WritePvarVariant(pvar_writer, record, minlen, maxlen):
    out_items = [record.CHROM, str(record.POS), record.ID, DUMMY_REF, DUMMY_ALT,
        "%.2f,%.2f"%(minlen, maxlen)]
    pvar_writer.write("\t".join(out_items)+"\n")

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--vcf", help="Input STR VCF file. Must be bgzipped/indexed", type=str, required=True)
    inout_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VcfTypes.__members__],
        type=str, default="auto")
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)
    inout_group.add_argument("--outtype", help="Options=%s"%[str(item) for item in OutputFileTypes.__members__],
        type=str, nargs="+", default=["vcf"])
    annot_group = parser.add_argument_group("Annotations")
    annot_group.add_argument(
        "--dosages", 
        help="Compute genotype dosages. " 
             "Optionally specify how. Options=%s"%[str(item) for item in trh.TRDosageTypes.__members__],
        type=str)
    annot_group.add_argument(
        "--ref-panel", 
        help="Annotate Beagle-imputed VCF with TR metadata from the reference panel. "
             "The reference must be the same VCF used for imputation. ", 
        type=str)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
    print(args.out)
    ###### Check input options #######
    if not os.path.exists(args.vcf):
        common.WARNING("Error: %s does not exist"%args.vcf)
        return 1
    if not os.path.exists(os.path.dirname(os.path.abspath(args.out))):
        common.WARNING("Error: The directory which contains the output location {} does"
                       " not exist".format(args.out))
        return 1
    if os.path.isdir(args.out) and args.out.endswith(os.sep):
        common.WARNING("Error: The output location {} is a "
                       "directory".format(args.out))
        return 1
    if args.ref_panel is not None and not os.path.exists(args.ref_panel):
        common.WARNING("Error: %s does not exist"%args.ref_panel)
        return 1

    outtypes = []
    for outtype in args.outtype:
        try:
            ot = OutputFileTypes[outtype]
            outtypes.append(ot)
        except KeyError:
            common.WARNING("Invalid output type")
            return 1

    if args.vcftype != 'auto':
        try:
            checktype = trh.VcfTypes[args.vcftype]
        except:
            raise ValueError("Error: invalid vcftype specified")

    ###### Load reference panel info (optional) #######
    refpanel_metadata = None
    refreader = None
    if args.ref_panel is not None:
        refreader = utils.LoadSingleReader(args.ref_panel)
        if refreader is None:
            return 1
        if args.vcftype != 'auto':
            refpanel_vcftype = trh.VcfTypes[args.vcftype]
        else:
            refpanel_vcftype = trh.InferVCFType(refreader)
        if refpanel_vcftype == trh.VcfTypes.popstr:
            common.WARNING("Error: reference panel annotation not "
                           "currently supported for popSTR")
            return 1
        refpanel_metadata = LoadMetadataFromRefPanel(refreader, refpanel_vcftype)
        if len(refpanel_metadata.keys()) == 0:
            common.WARNING("Error: No TRs detected in reference panel. Quitting")
            return 1
        common.MSG("Loaded " + str(len(refpanel_metadata.keys())) + " TR loci from ref panel",
            debug=True)

    ###### Load reader #######
    reader = utils.LoadSingleReader(args.vcf)
    if reader is None:
        return 1
    if args.ref_panel is not None:
        vcftype = refpanel_vcftype # should be same as refpanel
    elif args.vcftype != 'auto':
        vcftype = trh.VcfTypes[args.vcftype]
    else:
        vcftype = trh.InferVCFType(reader)

    ##### Additional checks on input #####
    dosage_type = None
    if args.dosages is not None:
        try:
            dosage_type = trh.TRDosageTypes[args.dosages]
        except KeyError:
            common.WARNING("Error: invalid dosages argument")
            return 1
    if dosage_type == trh.TRDosageTypes.beagleap and not trh.IsBeagleVCF(reader):
        common.WARNING("Error: can only compute beagleap dosages on Beagle VCFs")
        return 1
    if dosage_type is None and np.all([ot in [OutputFileTypes.pgen] for ot in outtypes]):
        common.WARNING("Error: Output type pgen only supported "
                       "if using option --dosages")
        return 1
    if dosage_type not in [trh.TRDosageTypes.beagleap_norm, trh.TRDosageTypes.bestguess_norm] and \
        OutputFileTypes.pgen in outtypes:
        common.WARNING("Only normalized dosages are supported for PGEN output.")
        return 1

    ###### Set up writers #######
    if OutputFileTypes.vcf in outtypes:
        vcf_writer = GetVCFWriter(reader, args.out+".vcf", " ".join(sys.argv), vcftype,
                                    dosage_type=dosage_type, refreader=refreader)
        if vcf_writer is None:
            common.WARNING("Error: problem initializing vcf writer.")
            return 1
    # variant_ct needed for pgen
    # If using a ref panel, assume we have same number
    # of TRs as the panel
    # Otherwise, assume our file is all TRs and use record count
    if refpanel_metadata is not None:
        variant_ct = len(refpanel_metadata.keys())
    else:
        variant_ct = reader.num_records
    if OutputFileTypes.pgen in outtypes:
        pgen_writer, pvar_writer = GetPGenPvarWriter(reader, args.out, variant_ct)

    ###### Process each record #######
    num_variants_processed_batch = 0
    num_variants_processed = 0
    dosages_batch = np.empty((DEFAULT_PGEN_BATCHSIZE, len(reader.samples)), dtype=np.float32)
    for record in reader:
        # If using refpanel, first add required fields
        # In that case, only process records in the refpanel
        # Otherwise, process all records in the input VCF
        if refpanel_metadata is not None:
            locuskey = "{chrom}:{pos}:{ref}".format(
                chrom=record.CHROM,
                pos=record.POS,
                ref=record.REF
            )
            if locuskey not in refpanel_metadata.keys():
                continue
            for infofield in INFOFIELDS[vcftype]:
                record.INFO[infofield] = refpanel_metadata[locuskey][infofield]
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
        minlen = trrecord.min_allele_length
        maxlen = trrecord.max_allele_length
        if dosage_type is not None:
            dosages = trrecord.GetDosages(dosage_type)
            # Update record
            record.INFO["DSLEN"] = "{minlen},{maxlen}".format(minlen=minlen, maxlen=maxlen)
            record.set_format("TRDS", np.array(dosages))
            # Update batch
            dosages_batch[num_variants_processed_batch] = dosages

        # Write to VCF if using vcf output
        if OutputFileTypes.vcf in outtypes:
            vcf_writer.write_record(record)

        # Write pvar if using pgen output
        if OutputFileTypes.pgen in outtypes:
            WritePvarVariant(pvar_writer, record, minlen, maxlen)

        num_variants_processed += 1
        num_variants_processed_batch += 1

        # Reset batch, and write to pgen if using that
        if ((num_variants_processed_batch == DEFAULT_PGEN_BATCHSIZE) \
            or (num_variants_processed==variant_ct)):
            # Write batch
            common.MSG("Processed {numvars} variants".format(numvars=num_variants_processed))
            if OutputFileTypes.pgen in outtypes:
                pgen_writer.append_dosages_batch(dosages_batch[:num_variants_processed_batch])
            # Reset
            dosages_batch = np.empty((DEFAULT_PGEN_BATCHSIZE, len(reader.samples)), dtype=np.float32)
            num_variants_processed_batch = 0

    ###### Cleanup #######
    if OutputFileTypes.pgen in outtypes:
        pgen_writer.close()
        pvar_writer.close()
    if OutputFileTypes.vcf in outtypes:
        vcf_writer.close()

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()