"""
Tool for annotating TR VCF files
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
DUMMY_QUAL = "."
DUMMY_FILTER = "."

# Info fields copied from reference panel for each tool
INFOFIELDS = {
    trh.VcfTypes.hipstr: ["START","END","PERIOD"],
    trh.VcfTypes.longtr: ["START","END","PERIOD"],
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

class RefMatchTypes(enum.Enum):
    """Different supported output file types."""
    locid = "locid"
    rawalleles = "rawalleles"
    trimmedalleles = "trimmedalleles"
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)

def CheckAlleleCompatibility(record_ref, record_alt, panel_ref, panel_alt):
    """
    Check if the REF and ALT alleles of the record and
    reference panel are compatible.

    In the case of running imputation with Beagle followed by
    bcftools merge, allele order is maintained but sequences
    themselves may be trimmed by bcftools. This causes problems
    when harmonizing HipSTR records, since the START/END coords
    are not updated accordingly. Using annotaTR option --update-ref-alt
    can restore the original allele sequences from the refpanel.
    This function provides basic checks to make sure the ref/alts
    of the panel and target VCF are compatible. In particular:
    - is the number of ALT alleles the same
    - are all alleles offset by the same number of bp
    - are all the ALTs in the target VCF substrings of the refpanel ALTS.
    If any of these fail then the alleles are deemed incompatible.

    Parameters
    ----------
    record_ref : str
       REF allele of the target VCF
    record_alt : list of str
       ALT alleles of the target VCF
    panel_ref : str
       REF allele of the ref panel
    panel_alt : list of str
       ALT alleles of the ref panel

    Returns
    -------
    is_compatible : bool
       True if all checks pass, otherwise False
    """
    if len(record_alt) != len(panel_alt):
        return False
    len_offset = len(panel_ref)-len(record_ref)
    for i in range(len(panel_alt)):
        if (len(panel_alt[i])-len(record_alt[i])) != len_offset:
            return False
        if record_alt[i].upper() not in panel_alt[i].upper():
            return False
    return True

def UpdateVCFHeader(reader, command, vcftype, dosage_type=None, refreader=None):
    """
    Update the VCF header of the reader to include:
    - The annotatTR command used
    - new INFO and FORMAT fields we will add if annotating dosages
    (INFO/DSLEN and FORMAT/TRDS)
    - new INFO fields we will add if using refpanel annotation. The fields added
    depend on the vcftype and are listed in INFOFIELDS

    Note this function gets called even if we are not producing VCF output
    since in some cases we still might need to add these fields to the record
    during processing, and cyvcf2 will throw an error if the proper headers
    are not there

    Parameters
    ----------
    reader : cyvcf2.VCF
        Reader for the input VCF file
    command : str
        The annotaTR command used
    vcftype : trh.VcfTypes
        Which type of TR VCF file the reader/refreader are
    dosage_type : trh.TRDosageType
        The type of dosages to be annotated. None if not computing dosages
    refreader : cyvcf2.VCF
        Reader for the reference panel

    Returns
    -------
    success : bool
       True if adding header fields was successful, otherwise False

    """
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
                # cyvcf2 might add a dummy header
                if headerinfo["Description"].replace('"','') == "Dummy":
                    common.WARNING("Could not find required header field {field} in refpanel".format(field=infofield))
                    return False
                reader.add_info_to_header({
                    'ID': headerinfo["ID"],
                    'Description': headerinfo["Description"].replace('"',''),
                    'Type': headerinfo["Type"],
                    "Number": headerinfo["Number"]
                    })
            else:
                common.WARNING("Could not find required header field {field} in refpanel".format(field=infofield))
                return False
    return True

def TrimAlleles(ref_allele, alt_alleles):
    """
    Trim ref and alt alleles to remove common prefixes
    and suffixes present in all alleles

    Parameters
    ----------
    ref_allele : str
       Reference allele
    alt_allele : list of str
       List of alternate alleles

    Returns
    -------
    new_ref_allele : str
       Trimmed reference allele
    new_alt_allele : list of str
       List of trimmed alternate alleles
    """
    # First check for longest common suffix
    alleles_left = [ref_allele] + alt_alleles
    longest_common_suffix = os.path.commonprefix([item[::-1] for item in alleles_left])[::-1]
    new_alt_alleles = alt_alleles.copy()
    new_ref_allele = ref_allele
    if len(longest_common_suffix) > 0:
        new_ref_allele = new_ref_allele[:-1*len(longest_common_suffix)]
        for i in range(len(new_alt_alleles)):
            new_alt_alleles[i] = new_alt_alleles[i][:-1*len(longest_common_suffix)]
    # Now remove longest common prefix
    longest_common_prefix = os.path.commonprefix([new_ref_allele]+new_alt_alleles)
    new_ref_allele = new_ref_allele[len(longest_common_prefix):]
    for i in range(len(new_alt_alleles)):
        new_alt_alleles[i] = new_alt_alleles[i][len(longest_common_prefix):]
    # Replace empty string with "."
    if new_ref_allele == "": new_ref_allele = "."
    for i in range(len(new_alt_alleles)):
        if new_alt_alleles[i] == "":
            new_alt_alleles[i] = "."
    return new_ref_allele, new_alt_alleles

def GetLocusKey(record, match_on=RefMatchTypes.locid):
    """
    Get the key used to match refpanel loci to the target VCF

    Options to match on:

    - RefMatchTypes.locid: use the ID from the VCF file
    - RefMatchTypes.rawalleles: use chrom:pos:ref:alt where ref/alt are
       exactly those in the reference VCF
    - RefMatchTypes.trimmedalleles: use chrom:pos:ref:alt where ref/alt are
       trimmed to discard extra sequence, as is done in bcftools merge :(
       see: https://github.com/samtools/bcftools/issues/726

    Parameters
    ----------
    record : cyvcf2.Variant
       Record to get the locus key for
    match_on : RefMatchTypes
       way to generate the key (Default: locid)

    Returns
    -------
    locuskey : str
       String of the key
    """
    if match_on == RefMatchTypes.locid:
        if record.ID is None or record.ID == ".":
            raise ValueError("Error: {chrom}:{pos} cannot match on loci ID if ID=.".format(
                chrom=record.CHROM, pos=record.POS))
        return record.ID
    elif match_on == RefMatchTypes.rawalleles:
        return "{chrom}:{pos}:{ref}:{alt}".format(
            chrom=record.CHROM,
            pos=record.POS,
            ref=record.REF,
            alt=",".join(record.ALT)
        )
    elif match_on == RefMatchTypes.trimmedalleles:
        ref, alt = TrimAlleles(record.REF, record.ALT)
        return "{chrom}:{pos}:{ref}:{alt}".format(
            chrom=record.CHROM,
            pos=record.POS,
            ref=ref,
            alt=",".join(alt)
        )
    else:
        raise ValueError("Invalid match_refpanel_on=%s"%match_on)

def LoadMetadataFromRefPanel(refreader, vcftype, match_on=RefMatchTypes.locid,
        ignore_duplicates=False):
    """
    Load required INFO fields from the ref panel we will use to
    annotate the target VCF. The specific INFO fields loaded
    depends on the vcftype and are specified in INFOFIELDS

    The value of match_on determines what to use as the key in the
    returned dictionary. Options:

    - RefMatchTypes.locid: use the ID from the VCF file
    - RefMatchTypes.rawalleles: use chrom:pos:ref:alt where ref/alt are
       exactly those in the reference VCF
    - RefMatchTypes.trimmedalleles: use chrom:pos:ref:alt where ref/alt are
       trimmed to discard extra sequence, as is done in bcftools merge :(
       see: https://github.com/samtools/bcftools/issues/726
    
    Parameters
    ----------
    refreader : cyvcf2.VCF
        Reader for the reference panel
    vcftype : trh.VcfTypes
        Based on the TR genotyper used to generate the reference panel
    match_on : RefMatchTypes (Optional)
        What to use as the locus key to match target to ref panel loci
        Default: RefMatchTypes.id
    ignore_duplicates : bool
        If True, just output a warning about duplicates rather than giving up

    Returns
    -------
    metadata : Dict[str, str]
        The key depends on the match_on parameter (see above)
        Values is a Dict[str, str] with key=infofield and
        value=value of that info field in the reference panel
        Also includes REF/ALT to check against alleles in imputed VCF
    variant_ct : int
        Total number of variants

    Raises
    ------
    ValueError
        If a duplicate locus is found in the reference panel
        and ignore_duplicates=False
    """
    metadata = {} # chr:pos:ref->info
    variant_ct = 0
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
        locuskey = GetLocusKey(record, match_on=match_on)
        # Quit if we found a duplicate TR locus
        if locuskey in metadata.keys():
            if ignore_duplicates:
                common.WARNING("Warning: duplicate locus detected in refpanel: {locus}".format(locus=locuskey))
            else:
                raise ValueError(
                    "Error: duplicate locus detected in refpanel: {locus}".format(locus=locuskey)
                    )
        locdata["REF"] = record.REF
        locdata["ALT"] = record.ALT
        metadata[locuskey] = locdata
        variant_ct += 1
    return metadata, variant_ct

def GetPGenPvarWriter(reader, outprefix, variant_ct):
    """
    Generate a PGEN and corresponding PVAR writer.
    For PGEN, we return a pgenlib.PgenWriter instance
    For PVAR, we create a file object to which we will write
    info for each variant as we go. When initialized here we 
    add the DSLEN INFO header and also the header columns:
    #CHROM", "POS", "ID", "REF", "ALT", "INFO"

    In addition to the PGEN/PVAR writers, this function writes
    $outprefix.psam wtih sample information.

    Parameters
    ----------
    reader : cyvcf2.VCF
        Reader for the input VCF file
    outprefix : str
        Prefix to name output files
        Will generate $outprefix.pgen and $outprefix.pvar
    variant_ct : int
        Number of variants to be written to the PGEN output

    Returns
    -------
    pgen_writer : pgenlib.PgenWriter
        PGEN writer object
    pvar_writer : file object
        Writer for the PVAR file
    """
    # Write .psam 
    with open(outprefix+".psam", "w") as f:
        f.write("#IID\tSEX\n")
        for sample in reader.samples:
            f.write("{sample}\t0\n".format(sample=sample))
    # Get pvar writer
    pvar_writer = open(outprefix+".pvar", "w")
    pvar_writer.write("##fileformat=VCFv4.2\n") # Required if we want the pvar to be valid vcf
    pvar_writer.write("##INFO=<ID=DSLEN,Number=2,Type=Float,Description=\"Minimum and maximum dosages, used if normalization was applied\">\n")
    pvar_writer.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])+"\n")
    # Get pgen writer we will continue writing to
    pgen_writer = PgenWriter(bytes(outprefix+".pgen", "utf8"),
        len(reader.samples), variant_ct=variant_ct, dosage_present=True)
    return pgen_writer, pvar_writer

def WritePvarVariant(pvar_writer, record, minlen, maxlen):
    """
    Write variant metadata to a PVAR file
    Outputs CHROM, POS, ID, REF, ALT, INFO
    REF and ALT are set to dummy values (DUMMY_REF and DUMMY_ALT)
    INFO contains the DSLEN field with the min/max allele lengths
    observed

    Parameters
    ----------
    pvar_writer : file object
        Writer for the PVAR file
    record : cyvcf2.Variant
        Record object for the variant
    minlen : float
        Minimum TR allele length at this locus (in rpt. units)
    maxlen : float
        Maximum TR allele length at this locus (in rpt. units)
    """
    record_id = record.ID
    if record_id is None:
        record_id = "."
    out_items = [record.CHROM, str(record.POS), str(record_id), DUMMY_REF, DUMMY_ALT,
        DUMMY_QUAL, DUMMY_FILTER,
        "DSLEN=%.2f,%.2f"%(minlen, maxlen)]
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
    inout_group.add_argument("--vcf-outtype", help="Type of VCF output to produce. "
                             "z=compressed VCF, v=uncompressed VCF, "
                             "b=compressed BCF, u=uncompressed BCF, s=stdout", type=str, default="v")
    inout_group.add_argument("--region", help="Restrict analysis to this region. Syntax: chr:start-end", type=str)
    annot_group = parser.add_argument_group("Annotations")
    annot_group.add_argument(
        "--dosages", 
        help="Compute genotype dosages. " 
             "Optionally specify how. Options=%s"%[str(item) for item in trh.TRDosageTypes.__members__],
        type=str)
    annot_group.add_argument(
        "--warn-on-AP-error",
        help="Output a warning but don't crash on error computing on AP field",
        action="store_true"
    )
    annot_group.add_argument(
        "--ref-panel", 
        help="Annotate Beagle-imputed VCF with TR metadata from the reference panel. "
             "The reference must be the same VCF used for imputation. ", 
        type=str)
    annot_group.add_argument(
        "--match-refpanel-on",
        help="What to match loci on between refpanel and target VCF. "
             "Options=%s"%[str(item) for item in RefMatchTypes.__members__],
        type=str,
        default="locid"
        )
    annot_group.add_argument(
        "--ignore-duplicates",
        help="Output a warning but do not crash if duplicate loci in refpanel",
        action="store_true"
        )
    annot_group.add_argument("--update-ref-alt", help="Update the REF/ALT allele sequences from the "
                                                      "reference panel. Fixes issue with alleles being "
                                                      "chopped after bcftools merge. Use with caution "
                                                      "as this assumes allele order is exactly the same "
                                                      "between the refpanel and target VCF. Only works when "
                                                      "matching on locus id", action="store_true")
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument(
        "--chunk-size",
        help="If writing a PGEN file, load dosages "
             "in chunks of X variants; reduces memory. ",
        type=int,
        default=DEFAULT_PGEN_BATCHSIZE)
    other_group.add_argument("--debug", help="Run in debug mode", action="store_true")
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
    if args.ref_panel is not None and not os.path.exists(args.ref_panel):
        common.WARNING("Error: %s does not exist"%args.ref_panel)
        return 1
    if args.match_refpanel_on != "locid" and args.update_ref_alt:
        common.WARNING("Error: you cannot use --update-ref-alt unless "
                       " --match-refpanel-on is set to locid")
        return 1
    if args.update_ref_alt and args.ref_panel is None:
        common.WARNING("Error: --update-ref-alt only works with "
                       " --ref-panel.")
        return 1

    outtypes = set()
    for outtype in args.outtype:
        try:
            ot = OutputFileTypes[outtype]
            outtypes.add(ot)
        except KeyError:
            common.WARNING("Invalid output type")
            return 1
    if args.vcf_outtype not in ["z","v","u","b","s"]:
        common.WARNING("Invalid VCF output type specified: {vcf_outtype}".format(vcf_outtype=args.vcf_outtype))
        return 1
    if args.vcftype != 'auto':
        if args.vcftype not in trh.VcfTypes.__members__:
            common.WARNING("Invalid vcftype")
            return 1

    dosage_type = None
    if args.dosages is not None:
        try:
            dosage_type = trh.TRDosageTypes[args.dosages]
        except KeyError:
            common.WARNING("Error: invalid dosages argument")
            return 1
    if dosage_type is None and np.all([ot in [OutputFileTypes.pgen] for ot in outtypes]):
        common.WARNING("Error: Output type pgen only supported "
                       "if using option --dosages")
        return 1
    if dosage_type not in [trh.TRDosageTypes.beagleap_norm, trh.TRDosageTypes.bestguess_norm] and \
        OutputFileTypes.pgen in outtypes:
        common.WARNING("Only normalized dosages are supported for PGEN output.")
        return 1
    if args.dosages is None and args.ref_panel is None:
        common.WARNING("No operation specified")
        return 1


    ###### Load reference panel info (optional) #######
    refpanel_metadata = None
    refreader = None
    if args.ref_panel is not None:
        common.MSG("Loading reference panel", debug=True)
        refreader = utils.LoadSingleReader(args.ref_panel, lazy=True, samples=set())
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
        if args.region is not None:
            refreader = refreader(args.region)
        try:
            match_on = RefMatchTypes[args.match_refpanel_on]
        except KeyError:
            common.WARNING("Invalid argument to --match-refpanel-on")
            return 1
        refpanel_metadata, ref_variant_ct = LoadMetadataFromRefPanel(refreader, refpanel_vcftype, 
            match_on=match_on, ignore_duplicates=args.ignore_duplicates)
        if len(refpanel_metadata.keys()) == 0:
            common.WARNING("Error: No TRs detected in reference panel. Check: "
                           "Was the right vcftype specified? "
                           "Was an invalid region specified? Quitting")
            return 1
        common.MSG("Loaded " + str(ref_variant_ct) + " TR loci from ref panel",
            debug=True)

    ###### Load reader #######
    reader = utils.LoadSingleReader(args.vcf, checkgz=True)
    if reader is None:
        return 1
    if args.ref_panel is not None:
        vcftype = refpanel_vcftype # should be same as refpanel
    elif args.vcftype != 'auto':
        vcftype = trh.VcfTypes[args.vcftype]
    else:
        vcftype = trh.InferVCFType(reader)

    ##### More checks on input #####
    # This check waits until here sicne need to have loaded reader
    # and need to wait til refpanel loading to confirm vcf type
    if dosage_type in [trh.TRDosageTypes.beagleap, trh.TRDosageTypes.beagleap_norm] \
        and not trh.IsBeagleVCF(reader):
        common.WARNING("Error: can only compute beagleap dosages on Beagle VCFs")
        return 1

    ###### Set up writers #######
    # Update reader header, even if not writing VCF output
    # This is because we might add VCF fields for parsing
    # with TRHarmonizer along the way
    # Note need to set up new refreader in case we set the region 
    # above in which case refreader is an iterator
    tmp_refreader = None
    if args.ref_panel is not None:
        tmp_refreader = utils.LoadSingleReader(args.ref_panel, lazy=True, samples=set())
    if not UpdateVCFHeader(reader, " ".join(sys.argv), vcftype,
                        dosage_type=dosage_type, 
                        refreader=tmp_refreader):
        common.WARNING("Error: problem initializing vcf header.")
        return 1
    if OutputFileTypes.vcf in outtypes:
        if args.vcf_outtype == "v":
            vcf_writer = cyvcf2.Writer(args.out+".vcf", reader)
        elif args.vcf_outtype == "z":
            vcf_writer = cyvcf2.Writer(args.out+".vcf.gz", reader, mode="wz")
        elif args.vcf_outtype == "b":
            vcf_writer = cyvcf2.Writer(args.out+".bcf", reader, mode="wb")
        elif args.vcf_outtype == "u":
            vcf_writer = cyvcf2.Writer(args.out+".bcf", reader, mode="wbu")
        elif args.vcf_outtype == "s":
            vcf_writer = cyvcf2.Writer("-", reader)
        else:
            raise ValueError("Encountered invalid VCF output type")
    # variant_ct needed for pgen
    # If using a ref panel, assume we have same number
    # of TRs as the panel
    # Otherwise, assume our file is all TRs and use record count
    if refpanel_metadata is not None:
        variant_ct = ref_variant_ct
    else:
        variant_ct = reader.num_records
    if OutputFileTypes.pgen in outtypes:
        pgen_writer, pvar_writer = GetPGenPvarWriter(reader, args.out, variant_ct)

    ###### Process each record #######
    num_variants_processed_batch = 0
    num_variants_processed = 0
    num_samples = len(reader.samples)
    dosages_batch = np.empty((args.chunk_size, num_samples), dtype=np.float32)
    if args.region:
        reader = reader(args.region)
    for record in reader:
        # If using refpanel, first add required fields
        # In that case, only process records in the refpanel
        # Otherwise, process all records in the input VCF
        if refpanel_metadata is not None:
            locuskey = GetLocusKey(record, match_on=match_on)
            if locuskey not in refpanel_metadata.keys():
                # If this looks like a TR, but not in our panel
                # give up since that is suspicous
                # Otherwise it is probably a SNP and we will skip
                try:
                    checkrec = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
                    common.WARNING("Error: Detected a TR {chrom}:{pos} not in refpanel".format(chrom=record.CHROM, pos=record.POS))
                    return 1
                except:
                    pass
                if args.debug:
                    common.WARNING("Detected locus not in refpanel: %s"%locuskey)
                continue
            for infofield in INFOFIELDS[vcftype]:
                record.INFO[infofield] = refpanel_metadata[locuskey][infofield]
            if args.update_ref_alt:
                # Update allele sequences to be exactly as in the
                # reference panel. 
                if not CheckAlleleCompatibility(record.REF, record.ALT,
                    refpanel_metadata[locuskey]["REF"], refpanel_metadata[locuskey]["ALT"]):
                    raise ValueError("--update-ref-alt set but the REF/ALT fields"
                                     " at {chrom}:{pos} are incompatible between the"
                                     " refpanel and target VCF".format(chrom=record.CHROM, pos=record.POS))
                record.REF = refpanel_metadata[locuskey]["REF"]
                record.ALT = refpanel_metadata[locuskey]["ALT"]
        try:
            trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
        except:
            common.WARNING("Error converting {chrom}:{pos} to a TR record. "
                "If your file is a mix of SNPs/TRs (e.g. from Beagle) you "
                "must provide a reference panel.".format(chrom=record.CHROM, pos=record.POS))
            return 1
        minlen = trrecord.min_allele_length
        maxlen = trrecord.max_allele_length
        # Add this check to warn us when bad things happen when parsing alleles
        if minlen == maxlen and len(trrecord.ref_allele) < 5:
            common.WARNING("Warning: Suspicious allele lengths found at "
                "{chrom}:{pos}. If you imputed then used bcftools merge "
                "and alleles were trimmed, consider using option "
                "--update-ref-alt. Otherwise dosage values may be invalid. "
                "Parsed alleles: ref={ref}, alt={alt}".format(chrom=record.CHROM, pos=record.POS, \
                    ref=trrecord.ref_allele, alt=",".join(trrecord.alt_alleles)))
        if dosage_type is not None:
            dosages = trrecord.GetDosages(dosage_type, strict=(not args.warn_on_AP_error))
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
        if ((num_variants_processed_batch == args.chunk_size) \
            or (num_variants_processed==variant_ct)):
            # Write batch
            common.MSG("Processed {numvars} variants".format(numvars=num_variants_processed), debug=True)
            if OutputFileTypes.pgen in outtypes:
                pgen_writer.append_dosages_batch(dosages_batch[:num_variants_processed_batch])
            # Reset
            dosages_batch = np.empty((args.chunk_size, num_samples), dtype=np.float32)
            num_variants_processed_batch = 0

    ###### Cleanup #######
    if OutputFileTypes.pgen in outtypes:
        try:
            pgen_writer.close()
        except RuntimeError:
            common.WARNING("Error writing PGEN! The output file is likely invalid. "
                "Did you run on files merged with bcftools merge? If so try rerunning "
                "with option --match-refpanel-on trimmedalleles or --match-refpanel-on locid.")
            return 1
        pvar_writer.close()
    if OutputFileTypes.vcf in outtypes:
        vcf_writer.close()
    return 0

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()