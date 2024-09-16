"""
Use SISTR to estimate selection coefficients
"""

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

from . import sistr_utils as sutils
from . import abc as abc

def GetAlleleFreqsList(tr_allele_freqs, opt_allele, numbins):
    # Get list of alleles
    allele_list = [round(allele-opt_allele) for allele in tr_allele_freqs.keys()]

    # Get highest absolute value in the list
    max_allele = abs(max(allele_list, key=abs))

    # Get list of allele freqs
    # opt allele is in the middle
    allele_freqs_list = [0]*(2*max_allele+1)
    for a in tr_allele_freqs.keys():
        allele_freqs_list[max_allele + round(a-opt_allele)] += tr_allele_freqs[a]

    # Add 0s to list of num alleles < num bins
    if len(allele_freqs_list) < numbins:
        num_zeros_to_add = int((numbins - len(allele_freqs_list))/2)
        allele_freqs_list = [0]*num_zeros_to_add + allele_freqs_list + \
            [0]*num_zeros_to_add
    return allele_freqs_list

def main(args):
    ###### Check input options #######
    if args.vcf is None:
        common.WARNING("Error: No VCF file specified")
        return 1
    if args.vcftype != 'auto':
        if args.vcftype not in trh.VcfTypes.__members__:
            common.WARNING("Invalid vcftype")
            return 1
    outtypes = set()
    for outtype in args.outtype:
        try:
            ot = sutils.SISTROutputFileTypes[outtype]
            outtypes.add(ot)
        except KeyError:
            common.WARNING("Invalid output type")
            return 1
    samples = []
    if args.samples_file is not None:
        if not os.path.exists(args.samples_file):
            common.WARNING("{sfile} does not exist".format(sfile=args.samples_file))
            return 1
        else:
            samples = [item.strip() for item in open(args.samples_file, "r").readlines()]
    if args.samples is not None:
        samples.extend(args.samples.split(","))
    if args.sistr_index is None:
        if (args.abc_lookup_folder is None or args.lrt_lookup_folder is None):
            common.WARNING("Must specify either --sistr-index or both --abc-lookup-folder"
                           " and --lrt-lookup-folder")
        return 1

    ###### Set up VCF reader #######
    reader = utils.LoadSingleReader(args.vcf, checkgz=True)
    if reader is None:
        return 1
    if args.vcftype != 'auto':
        vcftype = trh.VcfTypes[args.vcftype]
    else:
        vcftype = trh.InferVCFType(reader)
    if args.region is not None:
        reader = reader(args.region)
    if len(samples) > 0:
        all_samples = np.array(reader.samples)
        sample_index = np.isin(all_samples, samples)
    else: sample_index = None # process all samples

    ###### Setup SISTR ABC object #######
    if args.sistr_index:
        abc_index_folder = args.sistr_index
    else: abc_index_folder = args.abc_lookup_folder
    sistrABC = abc.SistrABC(
        sistr_index=abc_index_folder,
        minfreq=args.minfreq,
        numbins=args.numbins
    )

    ###### Set up writers #######
    # TODO
    
    ###### Process one locus at a time from VCF #######
    for record in reader:
        ##### Step 1: Obtain period, optimal allele and frequency info #####
        if args.verbose:
            common.WARNING("Processing {chrom}:{pos}".format(chrom=record.CHROM, pos=record.POS))
        trrecord = trh.HarmonizeRecord(vcfrecord=record, vcftype=vcftype)
        period = len(trrecord.motif)
        tr_allele_freqs = trrecord.GetAlleleFreqs(sample_index=sample_index) #numrpts->freq
        # Opt allele is most common allele (in num. copies of rpt.)
        opt_allele = round(max(tr_allele_freqs, key=tr_allele_freqs.get))
        # Format allele_freqs as needed for ABC
        # Opt allele is in the middle
        allele_freqs_list = GetAlleleFreqsList(tr_allele_freqs, opt_allele, args.numbins)
        if not sistrABC.CheckForModel(period, opt_allele):
            common.WARNING("Skipping {chrom}:{pos}. No model found for "
                           "period={period} optallele={opt}".format(
                                chrom=record.CHROM,
                                pos=record.POS,
                                period=period,
                                opt=opt_allele
                            ))
            continue

        ##### Step 2: Compute summary stats #####
        obs_summ_stats = abc.GetSummStats(allele_freqs_list, args.minfreq, args.numbins)

        ##### Step 3: Perform ABC - TODO #####

        # Step 4: Perform LRT - TODO

        # Step 5: Output to file - TODO

    ###### Clean up writers #######
    # TODO
    return 0
