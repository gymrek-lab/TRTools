"""
Generate the SISTR index
"""

import numpy as np
import sys

from . import sistr_utils as sutils
from . import sistr_sims as ssims

def GenerateIndex(period, opt_allele, config, outprefix):
	# Set up output file
	outf = open(outprefix + "_" + str(period) + "_" + str(opt_allele) + "_abc.txt", "w")

	# Draw s values from prior
	s_values = np.random.gamma(config["gamma_alpha"], config["gamma_beta"],
		size=config["abc_num_sims"])

	# Run each simulation
	for i in range(config["abc_num_sims"]):
		simres = ssims.RunSimulation(
			sval = s_values[i]
		)
		outf.write("\t".join([str(s_values[i]), simres["afreqs_string"]])+"\n")

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
			GenerateIndex(period, opt_allele, config, args.out)

	return 0