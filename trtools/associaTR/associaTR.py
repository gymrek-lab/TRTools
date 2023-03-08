#!/usr/bin/env python3

import argparse
import datetime
import shutil
import time
from typing import Optional
import sys

import cyvcf2
import numpy as np
import scipy.stats
import statsmodels.api as sm
from statsmodels.regression.linear_model import OLS
import statsmodels.stats.weightstats

from . import load_and_filter_genotypes
import trtools
import trtools.utils.utils as utils

pval_precision = 2

def _merge_arrays(a, b):
    '''
    Return a left outer join b.
    Join is performed on the first column.

    Assume first column of each array is id, but not necessarily same order

    Parameters
    ----------
    a, b: np.ndarray
        2D arrays
    '''
    
    # TODO these assertions should be moved elsewhere and made into understandable error messages
    assert len(a.shape) == 2 and len(b.shape) == 2
    assert len(set(a[:, 0]).intersection(b[:,0])) > 0
    assert len(set(a[:, 0])) == a.shape[0]
    assert len(set(b[:, 0])) == b.shape[0]

    b = b[np.isin(b[:, 0], a[:, 0])] # drop all elements of b that aren't in a
    matches = np.isin(a[:, 0], b[:, 0]) # a indicies that are in b

    a_sort = np.argsort(a[matches, 0])
    b_match_sorted = np.searchsorted(a[matches, 0], b[:, 0], sorter=a_sort)

    new_data = np.full((a.shape[0], b.shape[1] - 1), np.nan) # to hold the data from b
    new_data[matches, :] = b[np.argsort(b_match_sorted), 1:][np.argsort(a_sort), :]

    return np.concatenate((
        a,
        new_data
    ), axis=1)

def _weighted_binom_conf(weights, successes, confidence):
    r'''
    return a weighted binomial confidence interval
    using the Wilson method described here
    https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Wilson_score_interval
    except rederiving for a weighted sum \sum_i w_i X_i
    where the X_i are iid binomial variables

    Derivation starts with
    \sigma = \sqrt(p(1-p)}
    and
    \frac{ \sum_i w_i (X_i - p) } { \sigma \sqrt{ \sum_i w_i^2 } } \sim z
    where z \sim N(0, 1)
    substituting in \sigam, solving for p using the quadratic formula
    and getting
    p = \frac 1 {1 + c^2/t^2}(\hat p + \frac c^2 {2 t^2}) +/-
      \frac{c/t}{1 + c^2/t^2} \sqrt{\hat p (1 - \hat p) + c^2/(4t^2)}
    where t = \sum_i w_i (e.g. n for noneighted problem)
    \hat p = \sum_i w_i X_i / t
    c = z \sqrt{ \sum_i w_i^2 }
    z = F_{N(0,1)}^{-1}(1 - \alpha)
    where
    \alpha is the desired confidence interval

    returns (mean prob, lower ci bound, upper ci bound)

    spot testing on 2021/08/04 shows this code is the same as
    statsmodels.stats.proportion.proportion_confint(method='wilson')
    when all the weights are 1,
    that setting the weigths of a bunch of elements to near zero
    gives about the same effect as them not being there (but slightly
    smaller intervals)
    and setting the weights of a bunch of elements to 0.5 gives about
    the same effect as having half as many of those observations, but
    slightly smaller intervals (in the case where that class of observations
    is already overrepresented)
    '''
    assert weights.shape == successes.shape
    assert len(weights.shape) == 1
    t = np.sum(weights)
    phat = np.dot(weights, successes)/t
    z = scipy.stats.norm.ppf(1 - confidence/2)
    c = z * np.sqrt(np.dot(weights, weights))

    divisor = 2 + 2*c**2/t**2
    center = (2*phat + c**2/t**2)/divisor
    interval_size = c/t*np.sqrt(4*phat*(1-phat) + c**2/t**2)/divisor

    return (phat, center - interval_size, center + interval_size)

#    outfile,
#    traits_fnames,
#    untransformed_phenotypes_fname,
#    get_genotype_iter,
#    phenotype,
#    samples,
#    binary,
#    region,
#    runtype,
#):
def perform_gwas_helper(
    outfile,
    all_samples,
    get_genotype_iter,
    phenotype_name,
    trait_fnames,
    same_samples,
    sample_fname,
    beagle_dosages,
    plotting_phenotype_fname,
    paired_genotype_plot,
    plot_phenotype_residuals,
    plotting_ci_alphas
):
    outfile.write(
        "chrom\tpos\talleles\tn_samples_tested\tlocus_filtered\tp_{}\tcoeff_{}\t".format(phenotype_name, phenotype_name)
    )
    #if binary != 'logistic':
    outfile.write('se_{}\tregression_R^2\t'.format(phenotype_name))
    outfile.flush()
   
    print('{} samples in the VCF'.format(len(all_samples)), flush=True)

    if not same_samples:
        covars = np.load(trait_fnames[0])
        if np.sum(np.isin(np.array(all_samples, dtype=float), covars[:, 0])) < 3:
            print(all_samples, covars[:, 0])
            print(
                'Less than 3 samples matched between the covars array and the VCF. '
                'Prehaps you meant to run with --same-samples? '
                'Erroring out.'
            )
            exit(1)
        for trait_fname in trait_fnames[1:]:
            new_covars = np.load(trait_fname)
            covars = _merge_arrays(covars, new_covars)
        covars = _merge_arrays(np.array(all_samples, dtype=float).reshape(-1, 1), covars)
    else:
        covars_array_list = []
        for trait_fname in trait_fnames:
            covars_array_list.append(np.load(trait_fname))
            if not covars_array_list[-1].shape[0] == len(all_samples):
                print("different number of samples in covariates file {trait_fname} than VCF, "
                      "and --same-samples was specified. Erroring out."
                )
                sys.exit(1)
        # need an empty first column that will be used for the genotypes
        covars = np.hstack([np.full((covars_array_list[0].shape[0], 1), -1), *covars_array_list])

    if sample_fname:
        with open(sample_fname) as sample_file:
            sample_subset = [line.strip() for line in sample_file.readlines()]
            sample_filter = np.isin(all_samples, sample_subset)
            print((
                '{} samples remain after subsetting to samples '
                'from the file {}.\n'
                '{} samples from the sample file '
                'were not present in the VCF and were discarded.'
            ).format(np.sum(sample_filter), sample_fname, len(sample_subset) - np.sum(sample_filter)))
    else:
        sample_filter = np.array([True]*len(all_samples))

    prev_n_samples = sum(sample_filter)
    sample_filter = sample_filter & ~np.any(np.isnan(covars), axis=1)
    current_n_samples = sum(sample_filter)
    print((
        'Removing {} samples which had missing '
        'phenotypes or covariates.\n'
        'Using {} for the regression.\n'
        'The number of samples used in each variant\'s regression will only be lower '
        'if that variant has missing calls.\n'
    ).format(prev_n_samples - current_n_samples, current_n_samples))

    covars = covars[sample_filter, :]
    pheno_std = np.std(covars[:, 1])
    covars = (covars - np.mean(covars, axis=0))/np.std(covars, axis=0)
    outcome = covars[:, 1].copy()
    covars[:, 1] = 1 # reuse the column that was the outcome as the intercept
   
    if plotting_phenotype_fname:
        plotting_phenotype = np.load(plotting_phenotype_fname)
        if not same_samples:
            plotting_phenotype = _merge_arrays(
                np.array(all_samples, dtype=float).reshape(-1, 1), plotting_phenotype
            )[sample_filter, 1]
        else:
            plotting_phenotype = plotting_phenotype[sample_filter, 0]

    genotype_iter = get_genotype_iter(sample_filter.copy())

    # first yield is special
    extra_detail_fields = next(genotype_iter)
    outfile.write('\t'.join(extra_detail_fields) + '\n')

    #if not binary:
    #    stat = 'mean'
    #else:
    #    stat = 'fraction'
    stat = 'mean'

    if plotting_phenotype_fname:
        residual = 'residual_' if plot_phenotype_residuals else ''

        if not beagle_dosages:
            outfile.write('\tsample_count_per_summed_length')
        else:
            outfile.write('\ttotal_dosage_per_summed_length')
        outfile.write('\t{}_{}{}_per_summed_length'.format(stat, residual, phenotype_name))
        for alpha in plotting_ci_alphas:
            outfile.write('\tsummed_length_{:.2g}_alpha_CI'.format(alpha))

        if paired_genotype_plot:
            if not beagle_dosages:
                outfile.write('\tsample_count_per_length_pair')
            else:
                outfile.write('\ttotal_dosage_per_length_pair')
            outfile.write('\t{}_{}{}_per_length_pair'.format(stat, residual, phenotype_name))
            for alpha in plotting_ci_alphas:
                outfile.write('\tlength_pair_{:.2g}_alpha_CI'.format(alpha))

        outfile.write('\n')
        outfile.flush()

    n_loci = 0
    batch_time = 0
    batch_size = 50
    total_time = 0
 
    start_time = time.time()
    for gts, unique_alleles, chrom, pos, called_samples_filter, locus_filtered, locus_details in genotype_iter:
        assert len(locus_details) == len(extra_detail_fields)

        covars[:, 0] = np.nan # reuse the column that was the ids as the genotypes

        n_loci += 1
        allele_names = ','.join(list(unique_alleles.astype(str)))
        outfile.write(
            "{}\t{}\t{}\t{}\t".format(
            chrom, pos, allele_names, np.sum(called_samples_filter)
        ))
        if not locus_filtered and covars.shape[1] >= np.sum(called_samples_filter):
            locus_filtered = 'n covars >= n samples'
        if locus_filtered:
            outfile.write('{}\tnan\tnan\tnan\tnan\t'.format(locus_filtered))
            outfile.write('\t'.join(locus_details))
            n_nans = (2 + len(plotting_ci_alphas)) * (int(bool(plotting_phenotype_fname)) + int(bool(paired_genotype_plot)))
            outfile.write('\tnan'*n_nans + '\n')
            outfile.flush()
            continue
        else:
            outfile.write('False\t')

        if not beagle_dosages:
            summed_gts = np.sum(gts, axis=1)
        else:
            summed_gts = np.sum([
                len_*np.sum(dosages, axis=1) for len_, dosages in gts.items()
            ], axis=0)
        std = np.std(summed_gts)
        summed_gts = (summed_gts - np.mean(summed_gts))/np.std(summed_gts)
        covars[called_samples_filter, 0] = summed_gts

        #if not binary or binary == 'linear':
        #do da regression
        model = OLS(
            outcome[called_samples_filter],
            covars[called_samples_filter, :],
            missing='drop',
        )
        reg_result = model.fit()
        pval = reg_result.pvalues[0]
        coef = reg_result.params[0]
        se = reg_result.bse[0]
        rsquared = reg_result.rsquared
        outfile.write(("{:." + str(pval_precision) + "e}\t{}\t{}\t{}\t").format(pval, coef/std*pheno_std, se/std*pheno_std, rsquared))
#        else:
#            model = sm.GLM(
#                outcome,
#                covars,
#                missing='drop',
#                family=sm.families.Binomial()
#            )
#            reg_result = model.fit()
#            pval = reg_result.pvalues[0]
#            coef = reg_result.params[0]
#            outfile.write(f'{pval:.2e}\t{coef/std}\tnan\tnan\t')

        outfile.write('\t'.join(locus_details))

        # ----- plot phenotype statistics -----

        if plotting_phenotype_fname:
            if not plot_phenotype_residuals:
                phenotypes = plotting_phenotype
            else:
                #if not binary or binary == 'linear':
                #do da regression
                untrans_model = OLS(
                    plotting_phenotype,
                    covars[:, 1:],
                    missing='drop',
                )
        #        else:
        #            untrans_model = sm.GLM(
        #                ori_phenotypes,
        #                covars[:, 1:],
        #                missing='drop',
        #                family=sm.families.Binomial()
        #            )
                untrans_reg_result = untrans_model.fit()
                phenotypes = plotting_phenotype - untrans_reg_result.fittedvalues

            #if not binary:
            summed_lengths = {}
            genod_dicts = [summed_lengths]
            if paired_genotype_plot:
                paired_gts = {}
                geno_dicts.append(paired_gts)
            if not beagle_dosages:
                for summed_len in np.unique(summed_gts):
                    summed_lengths[summed_len] = summed_gts == summed_len
                if paired_genotype_plot:
                    sorted_gts = np.sort(gts, axis=1)
                    for pair in np.unique(sorted_gts, axis=0):
                        paired_gts[tuple(pair)] = \
                            (sorted_gts[:, 0] == pair[0]) & (sorted_gts[:, 1] == pair[1])
            else:
                for len1 in unique_alleles:
                    for len2 in unique_alleles:
                        if len1 > len2:
                            continue
                        if len1 != len2:
                            dosages = (dosage_gts[len1][:, 0]*dosage_gts[len2][:, 1] +
                                       dosage_gts[len1][:, 1]*dosage_gts[len2][:, 0])
                        else:
                            dosages = dosage_gts[len1][:, 0]*dosage_gts[len1][:, 1]
                        if np.sum(dosages) <= 0:
                            continue
                        summedlen_ = len1 + len2
                        if summedlen_ not in summed_lengths:
                            summed_lengths[summedlen_] = dosages
                        else:
                            summed_lengths[summedlen_] += dosages

                        minlen = min(len1, len2)
                        maxlen = max(len1, len2)
                        dosages_per_length_pair[(minlen, maxlen)] = dosages

            for geno_dict in geno_dicts:
                outfile.write('\t' + load_and_filter_genotypes.dict_str({key: np.sum(arr) for key, arr in geno_dict.items()}))
                stats = {}
                CIs = {alpha: {} for _ in plotting_ci_alphas}
                for len_, weights in geno_dict.items():
                    if len(np.unique(phenotypes[weights != 0])) <= 1:
                        stats[len_] = phenotype
                        for alpha in plotting_ci_alphas:
                            CIs[alpha][len_] = (np.nan, np.nan)
                        continue
                    mean_stats = statsmodels.stats.weightstats.DescrStatsW(
                        phenotypes,
                        weights = weights
                    )
                    stats[len_] = mean_stats.mean
                    for alpha in plotting_ci_alphas:
                        CIs[alpha][len_] = mean_stats.tconfint_mean(alpha)
                    # binary
        #            else:
        #                for len_, dosages in dosages_per_summed_length.items():
        #                    if not np.any(dosages != 0):
        #                        continue
        #                    p, lower, upper = _weighted_binom_conf.weighted_binom_conf(
        #                        dosages, phenotypes, 0.05
        #                    )
        #                    summed_length_stat[len_] = p
        #                    summed_length_95_CI[len_] = (lower, upper)
        #                    _,  lower_gwas, upper_gwas = weighted_binom_conf.weighted_binom_conf(
        #                        dosages, phenotypes, 5e-8
        #                    )
        #                    summed_length_GWAS_CI[len_] = (lower_gwas, upper_gwas)
                outfile.write('\t' + load_and_filter_genotypes.dict_str(summed_length_stat))
                for CI in CIs:
                    outfile.write('\t' + load_and_filter_genotypes.dict_str(summed_length_GWAS_CI))

        outfile.write('\n')
        outfile.flush()

        duration = time.time() - start_time
        total_time += duration
        batch_time += duration
        if n_loci % batch_size == 0:
            print((
                "time/locus (last {}): {}s\ntime/locus ({} total loci): {}s\n"
            ).format(batch_size, batch_time/batch_size, n_loci, total_time/n_loci),
                flush = True
            )
            batch_time = 0
        start_time = time.time()
    if n_loci > 0:
        print(
            "Done.\nTotal loci: {}\nTotal time: {}s\ntime/locus: {}s\n".format(
                n_loci, total_time, total_time/n_loci
            ),
            flush=True
        )
    else:
        print("No variants found in the region being looked at\n", flush=True)

def perform_gwas(
        outfname,
        tr_vcf,
        phenotype_name,
        traits_fnames,
        vcftype,
        same_samples,
        sample_fname,
        region,
        non_major_cutoff,
        beagle_dosages,
        plotting_phenotype_fname,
        paired_genotype_plot,
        plot_phenotype_residuals,
        plotting_ci_alphas,
        imputed_ukb_strs_paper_period_check
):

    all_samples = cyvcf2.VCF(tr_vcf).samples
    get_genotype_iter = lambda samples: load_and_filter_genotypes.load_trs(
        tr_vcf, samples, region, non_major_cutoff, beagle_dosages, vcftype,
        imputed_ukb_strs_paper_period_check
    )

    print("Writing output to {}.temp".format(outfname), flush=True)
    with open(outfname + '.temp', 'w') as outfile:
        perform_gwas_helper(
            outfile,
            all_samples,
            get_genotype_iter,
            phenotype_name,
            traits_fnames,
            same_samples,
            sample_fname,
            beagle_dosages,
            plotting_phenotype_fname,
            paired_genotype_plot,
            plot_phenotype_residuals,
            plotting_ci_alphas
        )

    print("Moving {}.temp to {}".format(outfname, outfname), flush=True)
    shutil.move(
        outfname + '.temp',
        outfname
    )
    print("Done.", flush=True)

def run():
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('outfile')
    parser.add_argument('tr_vcf')
    parser.add_argument('phenotype_name', help='name of the phenotype being regressed against')
    parser.add_argument(
        'traits', nargs='+',
        help='At least one (possibly more) .npy 2d float array files, containing trait values for samples. '
        'The first trait from the first file is the phenotype to be regressed against, all other traits '
        'from that file are used as covariates. Additional files can be listed to add additional covariates. '
        'If --same-samples is not specified, the first column of each file must be the numeric sample ID. '
        'So the phenotype will correspond to the second column from the first file. If there are multiple '
        'files, they will be joined on sample ID. '
        'If --same-samples is specified, there must be the same number of rows in each array as the number '
        'of samples in the vcf. In that case, the first column of the first array is the phenotype. If there '
        'are multiple files, then they will be concatenated horizontally. Since IDs do not need to be stored '
        'in the npy arrays, --same-samples allows for non-numeric sample IDs. '
        'Traits and the phenotype will be standardized to mean 0 and std 1 prior to regression, but '
        'coefficients/standard errors are transformed back to the original scale before being written out.'
    )
    parser.add_argument('--vcftype', choices=['eh', 'hipstr', 'gangstr', 'popstr', 'advntr'],
                       help="Specify which caller produced the TR VCF, useful when the VCF is ambiguous "
                            "and the caller cannot be automatically inferred.")
    parser.add_argument('--same-samples', default=False, action='store_true',
                        help='see the traits help string')
    #parser.add_argument('--binary', default=None, choices={'linear', 'logistic'}, type=Optional[str])
    parser.add_argument(
        '--sample-list',
        help="File containing list of samples to use, one sample ID per line. "
             "If not specified, all samples are used."
    )
    # TODO add variants list
    parser.add_argument('--region', help="Restrict to \"chr:start-end\"")
    parser.add_argument(
        '--non-major-cutoff', type=float, default=20,
        help='If not --beagle-dosages, then this is just the non-major-allele-count cutoff. '
             'I.e. filter all loci with non-major-allele-count < cutoff.'
             'If working with dosages, this cutoff is applied to the dosage sums. '
             'As with the regression itself, for this purpose alleles are coallesced by length. '
             'Default of 20 per plink\'s best practices: '
             'https://www.cog-genomics.org/plink/2.0/assoc#glm '
             'Set to 0 to disable this filter. '
    )
    parser.add_argument(
        '--beagle-dosages', action='store_true', default=False,
        help="regress against Beagle dosages from the AP{1,2} fields instead of from the GT field. "
             "(The GP field is not supported)"
    )
    parser.add_argument(
        '--plotting-phenotype',
        help=argparse.SUPPRESS
        #help="An npy array with the same format as traits. If specified, statistics "
        #     "will be output to allow for plotting TR loci against the first trait in this file "
        #     "(the plotting phenotype). All other triats in the file will be ignored. "
        #     "This is useful because often you will wish to regress against rank-inverse-normalized "
        #     "phenotypes, but the axis when plotting those is mostly meaningless, so this allows "
        #     "for plotting against the untransformed phenotypes. If unspecified then "
        #     "plotting phenotype statistics will not be computed or written out. "
        #     "If not using dosages, then samples are binned according to the summed value of their "
        #     "two alleles' lengths. If using dosages, then each sample contributes to each sum-bin "
        #     'proportionally to the estimated probability of their having that summed allele length.'
    )
    parser.add_argument(
        '--paired-genotype-plot',
        action='store_true',
        default=False,
        help=argparse.SUPPRESS
        #help='Only used for the --plotting-phenotype option. If True, then in addition to '
        #     'generating statistics for summed genotype values, also generate statistics for '
        #     'paired genotype values. This could be important if you wish to look for examples of '
        #     'nonadditive/nonlinear effects. '
        #     "If not using dosages, then samples are binned according to the unordered pair of their "
        #     "two alleles' lengths. If using dosages, then each sample contributes to each unordered-pair-bin "
        #     'proportionally to the estimated probability of their having that unoredered-pair genotype.'
    )
    parser.add_argument(
        '--plot-phenotype-residuals',
        action='store_true',
        default=False,
        help=argparse.SUPPRESS
        #help='Only used for the --plotting-phenotype option. If True, then linearlly regress out '
        #     'all covariates from the plotting phenotype before calculating the statistics. '
        #     'I would recommend against using this feature, as (a) in my small experience it has not '
        #     'affected the end results and (b) it is confusing, as this linear regression is '
        #     'performed on the plotting phenotype (not the phenotype used for regression) '
        #     'and if the plotting phenotype is untransformed (as is the intent) then an linear '
        #     'regression may not be appropriate, thus rendering this moot anyway. '
    )
    parser.add_argument(
        '--plotting-ci-alphas',
        type=float,
        nargs='*',
        default=[],
        help=argparse.SUPPRESS
        #help='Only used for the --plotting-phenotype option. Will generate these confidence intervals '
        #     'of size 1-alpha for each alpha in this list for each of the statistics generated for '
        #     'plotting.'
    )
    parser.add_argument(
        #This is only for one specific use case, don't use this
        '--imputed-ukb-strs-paper-period-check',
        default=False,
        action='store_true',
        help=argparse.SUPPRESS
    )
    parser.add_argument("--version", action="version", version = '{}'.format(trtools.__version__))

    args = parser.parse_args()
    main(args)

def main(args):
    today = datetime.datetime.now().strftime("%Y_%m_%d")
    print('-------Running AssociaTR (trtools v{}) ----------'.format(trtools.__version__))
    print("Run date: {}".format(today))
    print(args, flush=True)

    """
    # TODO properly implement binary
    else:
        readme.write('Doing logistic regressions against the binary phenotype. No longer '
                     'using firth penalized logistic regression when MAC <= 400, should but '
                     "this doesn't apply to any trs in this dataset. Instead, always using "
                     'standard logistic regression.\n')
    """

    perform_gwas(
        args.outfile,
        args.tr_vcf,
        args.phenotype_name,
        args.traits,
        args.vcftype,
        args.same_samples,
        args.sample_list,
        args.region,
        args.non_major_cutoff,
        args.beagle_dosages,
        args.plotting_phenotype,
        args.paired_genotype_plot,
        args.plot_phenotype_residuals,
        args.plotting_ci_alphas,
        args.imputed_ukb_strs_paper_period_check
    )


if __name__ == '__main__':
    run()

