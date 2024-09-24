"""
Generate the SISTR index
"""

import numpy as np
import stdpopsim
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

	# Set up demographic model
	if config["demog_model"] is not None:
		demog_model = stdpopsim.species.registered_species[config["species"]].get_demographic_model(config["demog_model"])
		popid = config["popid"]
		n_eff = None
	else:
		demog_model = None
		popid = None
		n_eff = config["n_effective"]
	# Function to run sims for a certain s value
	# since all other params are the same
	def RunSim(s):
		return ssims.RunSimulation(
			sval=s,
			transition_matrix_transpose=transition_matrix_transpose,
			max_iter=config["num_gens"],
			n_effective=n_eff,
			demog_model=demog_model, popid=popid,
			use_drift=config["use_drift"],
			end_samp_n=config["end_samp_n"]
		)

	########### LRT Lookup tables ################
	outf = open(outprefix + "_" + str(period) + "_" + str(opt_allele) + "_freqs.txt", "w")
	outf.write("\t".join(["s","freqs"])+"\n")
	for s in config["lrt_svals"]:
		allele_freq_results = []
		for i in range(config["lrt_num_sims"]):
			if verbose and i%100==0:
				sys.stderr.write("- LRT Simulation s=%s %s/%s\n"%(s, i, config["lrt_num_sims"]))
			simres = RunSim(s)
			if simres is None:
				sys.stderr.write("Error running simulation. Quitting")
				return False
			allele_freq_results.append(simres["afreqs_string"])
		outf.write("\t".join([str(s), ";".join(allele_freq_results)])+"\n")
	outf.close()

	########### LRT Lookup tables - zero ################
	outf = open(outprefix + "_" + str(period) + "_" + str(opt_allele) + "_zero_freqs.txt", "w")
	outf.write("\t".join(["s","freqs"])+"\n")
	allele_freq_results = []
	for i in range(config["lrt_num_sims"]):
		if verbose and i%100==0:
			sys.stderr.write("- LRT Simulation ZERO %s/%s\n"%(i, config["lrt_num_sims"]))
		simres = RunSim(0)
		allele_freq_results.append(simres["afreqs_string"])
	outf.write("\t".join(["0", ";".join(allele_freq_results)])+"\n")
	outf.close()
	
	########### ABC Lookup table ################
	# Set up output file
	outf = open(outprefix + "_" + str(period) + "_" + str(opt_allele) + ".txt", "w")
	outf.write("\t".join(["s","freqs"])+"\n")

	# Draw s values from prior
	s_values = np.random.gamma(config["gamma_alpha"], config["gamma_beta"],
		size=config["abc_num_sims"])

	# Run each simulation
	for i in range(config["abc_num_sims"]):
		if verbose and i%100==0:
			sys.stderr.write("- ABC Simulation %s/%s\n"%(i, config["abc_num_sims"]))
		simres = RunSim(s_values[i])
		outf.write("\t".join([str(s_values[i]), simres["afreqs_string"]])+"\n")
		outf.flush()
	outf.close()
	return True

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
			check = GenerateIndex(period, opt_allele, config, args.out, verbose=args.verbose)
			if not check:
				sys.stderr.write("Error generating tables")
				return 1

	return 0