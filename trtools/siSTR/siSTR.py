"""
Tool for performing selection inference at STRs

TODO
* allow separate input for demographic model file (sistr_sims.GetEffectivePopSize). maybe end_samp_n should go there
* organize options based on which of the commands they are needed for
* Index: where did the TMRCA.txt file come from for variable num generations? and why max 5920?
  I see that file in the google drive
* Simplify functions in sistr_sims.py
* score allow to process either freqs file or VCF file

Example:
sistr index --out testindex/test --verbose
"""

import argparse
import enum
import sys

import trtools.utils.common as common
import trtools.utils.utils as utils
import trtools.utils.tr_harmonizer as trh
from trtools import __version__

from . import index as index
from . import score as score
from . import sistr_utils as sutils

class SISTRCommands(enum.Enum):
    """Possible SISTR commands to run."""
    index = "index"
    score = "score"
    def __repr__(self):
        return '<{}.{}>'.format(self.__class__.__name__, self.name)

def getargs(): # pragma: no cover
    parser = argparse.ArgumentParser(
        __doc__,
        formatter_class=utils.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("command", help="Options=%s"%[str(item) for item in SISTRCommands.__members__])
    inout_group = parser.add_argument_group("Input/output")
    inout_group.add_argument("--out", help="Prefix for output files", type=str, required=True)
    score_group = parser.add_argument_group("Scoring options")
    score_group.add_argument("--vcf", help="Input VCF file to score", type=str)
    score_group.add_argument("--vcftype", help="Options=%s"%[str(item) for item in trh.VcfTypes.__members__],
        type=str, default="auto")
    score_group.add_argument("--outtype", help="Type of output file. Options=%s"%[str(item) for item in sutils.SISTROutputFileTypes.__members__], \
        type=str, nargs="+", default=["tab"])
    score_group.add_argument("--region", help="Restrict analysis to this region. Syntax: chr:start-end", type=str)
    score_group.add_argument("--samples", help="Restrict to comma separated list of samples", type=str)
    score_group.add_argument("--samples-file", help="Restrict to list of samples in this file", type=str)
    score_group.add_argument("--vcf-outtype", help="Type of VCF output. z=compressed VCF, v=uncompressed VCF, b=compressed BCF, u=uncompressed BCF, s=stdout", type=str, default="v")
    score_group.add_argument("--sistr-index", help="Folder of SISTR index. Alternatively specify" 
                             " both --abc-lookup-folder and --lrt-lookup-folder", type=str)
    score_group.add_argument("--abc-lookup-folder", help="Folder for ABC lookup tables", type=str)
    score_group.add_argument("--lrt-lookup-folder", help="Folder for LRT lookup tables", type=str)
    score_group.add_argument("--minfreq", help="Minimum allele frequency to count an allele as common", type=float, default=0.05)
    score_group.add_argument("--numbins", help="Number of bins for summarizing allele frequencies", type=int, default=5)
    score_group.add_argument("--eps-het-numerator", help="Numerator used to tune epsilon for heterozygosity comparisons", type=float, default=0.005)
    score_group.add_argument("--eps-het-denominator", help="Denominator used to tune epsilon for heterozygosity comparisons", type=float, default=3.0)
    score_group.add_argument("--eps-bins", help="Used to tune epsilon for bins comparisons", type=float, default=0.3)
    score_group.add_argument("--min-abc-acceptance", help="Do not report results for loci with fewer than this many ABC acceptances", type=int, default=10)
    index_group = parser.add_argument_group("Indexing options")
    index_group.add_argument(
        "--config",
        help="JSON file with indexing options. Index command line arguments "
             "override individual parameter values. See example: trtools/siSTR/config.json",
        type=str
    )
    index_group.add_argument(
        "--periods",
        help="Comma-separated list of rpt. unit lengths (in bp) to include. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["periods"]])),
        type=str,
    )
    index_group.add_argument(
        "--opt-allele-ranges",
        help="Comma-separated list of optimal allele ranges to include for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item[0])+"-"+str(item[1]) for item in sutils.DEFAULTS["opt_allele_ranges"]])),
        type=str
    )
    index_group.add_argument(
        "--log10-mut-slopes",
        help="Slopes of log10 mutation rate vs. allele length for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["log10_mut_slopes"]])),
        type=str
    )
    index_group.add_argument(
        "--betas",
        help="Beta values for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["betas"]])),
        type=str
    )
    index_group.add_argument(
        "--rhos",
        help="Rho values for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["rhos"]])),
        type=str
    )
    index_group.add_argument(
        "--baseline-mus",
        help="Baseline mut.rate values for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["baseline_mus"]])),
        type=str
    )
    index_group.add_argument(
        "--baseline-mu-alleles",
        help="Alleles corresponding to baseline mut.rate values for each rpt. unit length. "
             "If not set defaults to {default}".format(default=",".join([str(item) for item in sutils.DEFAULTS["baseline_mu_alleles"]])),
        type=str
    )
    index_group.add_argument(
        "--n-effective",
        help="Effective population size. "
             "If not set defaults to {default}".format(default=sutils.DEFAULTS["n_effective"]),
        type=int
    )
    index_group.add_argument(
        "--num-gens",
        help="Maximum number of generations to simulate. "
             "If not set defaults to {default}".format(default=sutils.DEFAULTS["num_gens"]),
        type=int
    )
    index_group.add_argument(
        "--num-alleles",
        help="Number of possible alleles in simulations. "
             "If not set defaults to {default}".format(default=sutils.DEFAULTS["num_alleles"]),
        type=int
    )
    index_group.add_argument(
        "--s-prior-gamma-params",
        help="Alpha and beta for prior distribution on s. "
             "If not set defaults to {a},{b}".format(a=sutils.DEFAULTS["gamma_alpha"],
                b=sutils.DEFAULTS["gamma_beta"]),
        type=str
    )
    index_group.add_argument(
        "--abc-num-sims",
        help="Number of simulations for ABC lookup tables. "
             "If not set defaults to {numsim}".format(numsim=sutils.DEFAULTS["abc_num_sims"]),
        type=int
    )
    index_group.add_argument(
        "--min-mu",
        help="Don't allow mutation rates to go below this value. "
             "If not set defaults to {minmu}".format(minmu=sutils.DEFAULTS["min_mu"]),
        type=float
    )
    index_group.add_argument(
        "--max-mu",
        help="Don't allow mutation rates to go above this value. "
             "If not set defaults to {maxmu}".format(maxmu=sutils.DEFAULTS["max_mu"]),
        type=float
    )
    index_group.add_argument(
        "--dont-use-drift",
        help="Don't apply drift in sampling",
        action="store_true"
    )
    index_group.add_argument(
        "--end-samp-n",
        help="Sample size to use for end sampling step. "
             "If not set defaults to {endsampn}".format(endsampn=sutils.DEFAULTS["end_samp_n"]),
        type=int
    )
    index_group.add_argument(
        "--lrt-num-sims",
        help="Number of simulations to use for LRT lookup tables. "
             "If not set defaults to {lrtnumsim}".format(lrtnumsim=sutils.DEFAULTS["lrt_num_sims"]),
        type=int
    )
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument("--seed", help="Set seed for random number generation", type=int)
    other_group.add_argument("--version", action="version", version = '{version}'.format(version=__version__))
    other_group.add_argument("--verbose", help="Output helpful status messages", action="store_true")
    args = parser.parse_args()
    return args

def main(args):
    try:
        command = SISTRCommands[args.command]
    except KeyError:
        common.WARNING("Command {cmd} invalid".format(cmd=args.command))
        sys.exit(1)

    if command == SISTRCommands.index:
        return index.main(args)
    if command == SISTRCommands.score:
        return score.main(args)

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()

