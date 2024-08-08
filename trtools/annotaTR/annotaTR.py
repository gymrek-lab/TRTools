"""
Tool for annotating TR VCF files

TODO:

* add min/max len to pgen output (ask Tara about header)
* Is beagle annotation relevant for pgen output? add to the pvar file?
* Ask Tara why we had to skip loci with no period? probably bc not STRs?
* Add documentation to functions
* Add tests
* Add README and link to other docs
* Add beagle annotation
* Does input contain both SNPs and TRs for Beagle output? if so skip SNPs. 
this will change variant_ct for pgen writer so might need to do beagle pass first
to count


Beagle annotation steps:
1. Modify header of "imputed""
- change source/command header lines from ref panel to "preimputation_source" and "preimputation_command"
- copy contig, ALT, INFO/END lines from ref panel

2. Add relevant info fields from refpanel
advntr: RU VID
eh: RU VARID RL
gangstr: RU
hipstr: START END PERIOD

For beagle annotation can we re-use mergeSTR code to walk through in sorted order?
or if not too slow can build a big dictionary of chr:pos:ref->info
this will help precompute variant count for pgen

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

DEFAULT_PGEN_BATCHSIZE = 1000
DUMMY_REF = "A"
DUMMY_ALT = "T"

class OutputFileTypes(enum.Enum):
    """Different supported output file types."""
    vcf = "vcf"
    pgen = "pgen"
    tensorqtl = "tensorqtl"
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)

def GetVCFWriter(reader, fname, command, dosage_type=None):
    reader.add_to_header("##command-AnnotaTR=" + command)
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
    writer = cyvcf2.Writer(fname, reader)
    return writer

def GetPGenPvarWriter(reader, outprefix):
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
        len(reader.samples), variant_ct=reader.num_records, dosage_present=True)
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
        type=str, nargs="+", default="vcf")
    annot_group = parser.add_argument_group("Annotations")
    annot_group.add_argument(
        "--dosages", 
        help="Compute genotype dosages. " 
             "Optionally specify how. Options=%s"%[str(item) for item in trh.TRDosageTypes.__members__],
        type=str)
    annot_group.add_argument(
        "--annotate-beagle", 
        help="Annotate Beagle-imputed VCF with TR metadata from the reference panel",
        type=str)
    ver_group = parser.add_argument_group("Version")
    ver_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    args = parser.parse_args()
    return args

def main(args):
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
    if args.annotate_beagle is not None and not os.path.exists(args.annotate_beagle):
        common.WARNING("Error: %s does not exist"%args.annotate_beagle)
        return 1

    outtypes = []
    for outtype in args.outtype:
        try:
            ot = OutputFileTypes[outtype]
            outtypes.append(ot)
        except KeyError:
            common.WARNING("Invalid output type")
            return 1

    ###### Load reader #######
    reader = utils.LoadSingleReader(args.vcf)
    if reader is None:
        return 1
    if args.vcftype != 'auto':
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
    if dosage_type is None and np.all([ot in [OutputFileTypes.pgen, OutputFileTypes.tensorqtl] for ot in outtypes]):
        common.WARNING("Output types pgen and tensorflow only supported "
                       "if using option --dosages")
        return 1

    ###### Set up writers #######
    if OutputFileTypes.vcf in outtypes:
        vcf_writer = GetVCFWriter(reader, args.out+".vcf", " ".join(sys.argv),
                                    dosage_type=dosage_type)
    if OutputFileTypes.pgen in outtypes:
        pgen_writer, pvar_writer = GetPGenPvarWriter(reader, args.out)
    if OutputFileTypes.tensorqtl in outtypes:
        common.WARNING("TensorQTL output not yet implemented")
        return 1

    ###### Process each record #######
    num_variants_processed_batch = 0
    num_variants_processed = 0
    dosages_batch = np.empty((DEFAULT_PGEN_BATCHSIZE, len(reader.samples)), dtype=np.float32)
    for record in reader:
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
        minlen = None
        maxlen = None
        if dosage_type is not None:
            dosages, minlen, maxlen = trrecord.GetDosages(dosage_type)
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
            or (num_variants_processed==reader.num_records)):
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