"""
Tool for annotating TR VCF files

TODO:
* make sure to store min/max len used for each in pvar
  or VCF INFO field. how to add header/fields?
* Ask Tara why we had to skip loci with no period?
* Add documentation to functions
* Add tests
* Add README and link to other docs
* Add appropriate headers to VCF (e.g. dosage field)
* Add beagle annotation
"""

import argparse
import cyvcf2
import enum
import math
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

def GetVCFWriter(reader, fname, command):
    reader.add_to_header("##command-AnnotaTR=" + command)
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
    pvar_writer.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT"])+"\n")
    # Get pgen writer we will continue writing to
    pgen_writer = PgenWriter(bytes(outprefix+".pgen", "utf8"),
        len(reader.samples), variant_ct=reader.num_records, dosage_present=True)
    return pgen_writer, pvar_writer

def WritePvarVariant(pvar_writer, record):
    out_items = [record.CHROM, str(record.POS), record.ID, DUMMY_REF, DUMMY_ALT]
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
        type=str, default="vcf")
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
    if dosage_type is None and args.outtype in ["pgen","tensorqtl"]:
        common.WARNING("Output types pgen and tensorflow only supported "
                       "if using option --dosages")
        return 1

    ###### Set up writers #######
    try:
        outtype = OutputFileTypes[args.outtype]
    except KeyError:
        common.WARNING("Invalid output type")
        return 1
    if outtype == OutputFileTypes.vcf:
        vcf_writer = GetVCFWriter(reader, args.out+".vcf", " ".join(sys.argv))
    elif outtype == OutputFileTypes.pgen:
        pgen_writer, pvar_writer = GetPGenPvarWriter(reader, args.out)
    elif outtype == OutputFileTypes.tensorqtl:
        common.WARNING("TensorQTL output not yet implemented")
        return 1
    else:
        common.WARNING("Invalid output type specified")
        return 1

    ###### Process each record #######
    num_variants_processed_batch = 0
    num_variants_processed = 0
    dosages_batch = np.empty((DEFAULT_PGEN_BATCHSIZE, len(reader.samples)), dtype=np.float32)
    for record in reader:
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
        if dosage_type is not None:
            dosages = trrecord.GetDosages(dosage_type) # TODO change to dosages!!!
            print(dosages)
            # Update record - TODO
            # Update batch
            dosages_batch[num_variants_processed_batch] = dosages

        # Write to VCF if using vcf output
        if outtype == OutputFileTypes.vcf:
            vcf_writer.write_record(record)

        # Write pvar if using pgen output
        if outtype == OutputFileTypes.pgen:
            WritePvarVariant(pvar_writer, record)

        num_variants_processed += 1
        num_variants_processed_batch += 1

        # Reset batch, and write to pgen if using that
        if ((num_variants_processed_batch == DEFAULT_PGEN_BATCHSIZE) \
            or (num_variants_processed==reader.num_records)):
            # Write batch
            common.MSG("Processed {numvars} variants".format(numvars=num_variants_processed))
            if outtype == OutputFileTypes.pgen:
                pgen_writer.append_dosages_batch(dosages_batch[:num_variants_processed_batch])
            # Reset
            dosages_batch = np.empty((DEFAULT_PGEN_BATCHSIZE, len(reader.samples)), dtype=np.float32)
            num_variants_processed_batch = 0

    ###### Cleanup #######
    if outtype == OutputFileTypes.pgen:
        pgen_writer.close()
        pvar_writer.close()
    if outtype == OutputFileTypes.vcf:
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