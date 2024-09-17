"""
Use SISTR to estimate selection coefficients
"""

import trtools.utils.common as common
import trtools.utils.tr_harmonizer as trh
import trtools.utils.utils as utils

from . import sistr_utils as sutils
from . import abc as abc

NASTRING = "N/A"

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

    ###### Setup SISTR index #######
    if args.sistr_index:
        abc_index_folder = args.sistr_index
        lrt_index_folder = args.sistr_index
    else:
        abc_index_folder = args.abc_lookup_folder
        lrt_index_folder = args.lrg_lookup_folder
    sistrABC = abc.SistrABC(
        abc_index=abc_index_folder,
        lrt_index=lrt_index_folder,
        minfreq=args.minfreq,
        numbins=args.numbins,
        eps_het_numerator=args.eps_het_numerator,
        eps_het_denominator=args.eps_het_denominator,
        eps_bins=args.eps_bins,
        min_abc_acceptance=args.min_abc_acceptance
    )

    ###### Set up writers #######
    if sutils.SISTROutputFileTypes.tab in outtypes:
        tab_writer = open(args.out + ".tab", "w")
        tab_writer.write("\t".join(["chrom", "start", "end", "total", "period", \
                "optimal_ru", "motif", "ABC_s_median", "ABC_s_95%_CI", \
                "Percent_s_accepted", "Likelihood_0", "Likelihood_s", \
                "LR", "LogLR", "LRT_p_value"])+"\n")
    # TODO - VCF writer

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
        # Opt allele is in the middle, other indices relative to opt
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

        ##### Step 3: Perform ABC #####
        abc_results = sistrABC.RunABC(obs_summ_stats, period, opt_allele)

        ##### Step 4: Perform LRT #####
        if abc_results["passed"]:
            ABC_conf_int = "(%.5f,%.5f)"%(abc_results["lower_bound"], abc_results["upper_bound"])
            lrt_results = sistrABC.LikelihoodRatioTest(abc_results["median_s"], obs_summ_stats, \
                period, opt_allele)
        else:
            lrt_results = {}
            ABC_conf_int = "N/A"

        ##### Step 5: Output to file #####
        if sutils.SISTROutputFileTypes.tab in outtypes:
            items = [
                record.CHROM,
                str(record.POS),
                str(record.INFO["END"]),
                str(record.num_called),
                str(period),
                str(opt_allele),
                trrecord.motif,
                str(abc_results.get("median_s", NASTRING)),
                ABC_conf_int,
                str(abc_results.get("num_accepted", NASTRING)),
                str(lrt_results.get("likelihood_0", NASTRING)),
                str(lrt_results.get("likelihood_s_ABC", NASTRING)),
                str(lrt_results.get("LR", NASTRING)),
                str(lrt_results.get("LogLR", NASTRING)),
                str(lrt_results.get("pval", NASTRING))
            ]
            tab_writer.write("\t".join(items)+"\n")

        # TODO - VCF
        # locus-level stats like tab output, but also
        # score individual alleles

    ###### Clean up writers #######
    if sutils.SISTROutputFileTypes.tab in outtypes:
        tab_writer.close()

    # TODO - VCF
    return 0
