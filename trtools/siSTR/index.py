"""
Generate the SISTR index
"""

import numpy as np
import sys

from . import sistr_utils as sutils
from . import sistr_sims as ssims

def GenerateIndex(period, opt_allele, config, outprefix, verbose=False):
	# Extract params
	ind = config["periods"].index(period)
	beta = config["betas"][ind]
	rho = config["rhos"][ind]
	L = config["log10_mut_slopes"][ind]
	baseline_mu = config["baseline_mus"][ind]
	baseline_mu_allele = config["baseline_mu_alleles"][ind]
	mu_prime = ssims.GetMuPrime(baseline_mu, baseline_mu_allele, L, opt_allele,
		config["min_mu"], config["max_mu"])

	# Set up transition matrix (constant)
	transition_matrix_transpose = ssims.GetTransitionMatrix(
		config["num_alleles"], mu_prime, beta, rho, L,
		config["min_mu"], config["max_mu"]).transpose()

	# Set up output file
	outf = open(outprefix + "_" + str(period) + "_" + str(opt_allele) + "_abc.txt", "w")
	outf.write("\t".join(["s","freqs"])+"\n")

	# Draw s values from prior
	s_values = np.random.gamma(config["gamma_alpha"], config["gamma_beta"],
		size=config["abc_num_sims"])

	# Run each simulation
	for i in range(config["abc_num_sims"]):
		if verbose and i%100==0:
			sys.stderr.write("- Simulation %s/%s\n"%(i, config["abc_num_sims"]))
		simres = ssims.RunSimulation(
			sval = s_values[i],
			transition_matrix_transpose = transition_matrix_transpose,
			max_iter=config["num_gens"],
			n_effective=config["n_effective"],
			use_drift=config["use_drift"],
			end_samp_n=config["end_samp_n"]
		)
		outf.write("\t".join([str(s_values[i]), simres["afreqs_string"]])+"\n")
		outf.flush()
	# Done
	outf.close()

def main(args):
	if args.seed is not None:
		np.random.seed(args.seed)
	#### Load configuration ####
	config = sutils.LoadSISTRConfig(args)
	if config is None:
		return 1
	if args.verbose:
		sutils.PrintConfigInfo(config)
	sutils.WriteConfig(config, args.out + ".json")

	#### Generate index for each period/opt allele ####
	for i in range(len(config["periods"])):
		period = config["periods"][i]
		min_opt_allele = config["opt_allele_ranges"][i][0]
		max_opt_allele = config["opt_allele_ranges"][i][1]
		for opt_allele in range(min_opt_allele, max_opt_allele+1):
			if args.verbose:
				sys.stderr.write("Generating index for period=%s opt=%s\n"%(period, opt_allele))
			GenerateIndex(period, opt_allele, config, args.out, verbose=args.verbose)

	return 0