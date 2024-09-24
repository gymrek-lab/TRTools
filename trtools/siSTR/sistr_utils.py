"""
Common utilities for SISTR commands
"""

import enum
import json
import stdpopsim
import sys

import trtools.utils.common as common

# Default SISTR params
DEFAULTS = {
	"periods": [2, 3, 4],
	"opt_allele_ranges": [(11, 20), (5, 13), (7, 10)],
	"species": None,
	"demog_model": None,
	"popid": None,
	"log10_mut_slopes": [0.15, 0.33, 0.45],
	"betas": [0.3, 0.3, 0.3],
	"rhos": [0.6, 0.9, 0.9],
	"baseline_mus": [10**-5, 10**-5.5, 10**-6],
	"baseline_mu_alleles": [6, 5, 3],
	"n_effective": 7310,
	"num_gens": 55920,
	"num_alleles": 25,
	"gamma_alpha": 0.0881,
	"gamma_beta": 0.2541,
	"abc_num_sims": 10000,
	"min_mu": 10**-8,
	"max_mu": 10**-3,
	"use_drift": True,
	"end_samp_n": 6500,
	"lrt_svals": [0.0,1e-06,1e-05,2e-05,3e-05,4e-05,5e-05,6e-05,7e-05,8e-05,9e-05,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.4,0.6,0.8,1.0],
	"lrt_num_sims": 2000
}

def WriteConfig(config, fname):
	with open(fname, "w") as f:
		json.dump(config, f, indent=4)

def PrintConfigInfo(config):
	sys.stderr.write("**** Loaded SISTR config info: *****\n")
	sys.stderr.write("Mutation models:\n")
	for i in range(len(config["periods"])):
		sys.stderr.write("  Period={period}\n".format(period=config["periods"][i]))
		sys.stderr.write("    Optimal allele range={minval}-{maxval}\n".format(
				minval=config["opt_allele_ranges"][i][0],
				maxval=config["opt_allele_ranges"][i][1]
			))
		sys.stderr.write("    Log10 mut slope={slope}\n".format(slope=config["log10_mut_slopes"][i]))
		sys.stderr.write("    Beta={beta}\n".format(beta=config["betas"][i]))
		sys.stderr.write("    Rho={rho}\n".format(rho=config["rhos"][i]))
		sys.stderr.write("    Baseline mu={bmu}\n".format(bmu=config["baseline_mus"][i]))
		sys.stderr.write("    Baseline mu allele={bmua}\n".format(bmua=config["baseline_mu_alleles"][i]))
	sys.stderr.write("Demographic model:\n")
	sys.stderr.write("    species={species}".format(species=config["species"]))
	sys.stderr.write("    demog_model={demogmodel}".format(demogmodel=config["demog_model"]))
	sys.stderr.write("    popid={popdid}".format(popid=config["popid"]))
	sys.stderr.write("    n_effective={neff}\n".format(neff=config["n_effective"]))
	sys.stderr.write("    num_gens={numgens}\n".format(numgens=config["num_gens"]))
	sys.stderr.write("Simulation params:\n")
	sys.stderr.write("    num_alleles={numalleles}\n".format(numalleles=config["num_alleles"]))
	sys.stderr.write("    gamma_params={a},{b}\n".format(a=config["gamma_alpha"], b=config["gamma_beta"]))
	sys.stderr.write("    abc_num_sims={numsim}\n".format(numsim=config["abc_num_sims"]))
	sys.stderr.write("    lrt_num_sims={numsim}\n".format(numsim=config["lrt_num_sims"]))
	sys.stderr.write("    min_mu={min_mu}\n".format(min_mu=config["min_mu"]))
	sys.stderr.write("    max_mu={max_mu}\n".format(max_mu=config["max_mu"]))
	sys.stderr.write("    use_drift={drift}\n".format(drift=config["use_drift"]))
	sys.stderr.write("    end_samp_n={endsampn}\n".format(endsampn=config["end_samp_n"]))
	sys.stderr.write("    lrt_svals={svals}\n".format(svals=config["lrt_svals"]))
	sys.stderr.write("************************************\n")


def LoadSISTRConfig(args):
	"""
	Load SISTR configuration parameters.
	If a config file is specified, load from JSON.
	Override JSON provided values or defaults with
	command-line options

	Return None if we encountered a problem parsing config
	"""

	# Set up default config
	config = DEFAULTS

	# Load from JSON
	if args.config is not None:
		json_config = json.load(open(args.config, "r"))
		for key in DEFAULTS:
			if key in json_config.keys():
				config[key] = json_config[key]
		for key in json_config.keys():
			if key not in DEFAULTS.keys():
				common.WARNING("Warning: config key {key} ignored.".format(key=key))

	# Override any command line defaults
	if args.periods is not None:
		config["periods"] = [int(item) for item in args.periods.strip().split(",")]
	if args.opt_allele_ranges is not None:
		try:
			config["opt_allele_ranges"] = [(int(item.split("-")[0]), int(item.split("-")[1])) \
				for item in args.opt_allele_ranges.strip().split(",")]
		except (IndexError, ValueError) as e:
			common.WARNING("Error parsing opt-allele-ranges")
			return None
	if args.species is not None:
		config["species"] = args.species
	if args.demographic_model is not None:
		config["demog_model"] = args.demographic_model
	if args.popid is not None:
		config["popid"] = args.popid
	if args.log10_mut_slopes is not None:
		config["log10_mut_slopes"] = [float(item) for item in args.log10_mut_slopes.strip().split(",")]
	if args.betas is not None:
		config["betas"] = [float(item) for item in args.betas.strip().split(",")]
	if args.rhos is not None:
		config["rhos"] = [float(item) for item in args.rhos.strip().split(",")]
	if args.baseline_mus is not None:
		config["baseline_mus"] = [float(item) for item in args.rhos.strip().split(",")]
	if args.baseline_mu_alleles is not None:
		config["baseline_mu_alleles"] = [int(item) for item in args.periods.strip().split(",")]
	if args.n_effective is not None:
		config["n_effective"] = args.n_effective
	if args.num_gens is not None:
		config["num_gens"] = args.num_gens
	if args.num_alleles is not None:
		config["num_alleles"] = args.num_alleles
	if args.s_prior_gamma_params is not None:
		items = args.s_prior_gamma_params.strip().split(",")
		try:
			alpha = float(items[0])
			beta = float(items[1])
		except (IndexError, TypeError) as e:
			common.WARNING("Error parsing alpha,gamma from " + args.s_prior_gamma_params)
			return None
		config["gamma_alpha"] = alpha
		config["gamma_beta"] = beta
	if args.abc_num_sims is not None:
		config["abc_num_sims"] = args.abc_num_sims
	if args.min_mu is not None:
		config["min_mu"] = args.min_mu
	if args.max_mu is not None:
		config["max_mu"] = args.max_mu
	if args.dont_use_drift:
		config["use_drift"] = False
	if args.end_samp_n is not None:
		config["end_samp_n"] = args.end_samp_n
	if args.lrt_num_sims is not None:
		config["lrt_num_sims"] = args.lrt_num_sims

	# Checks on values
	if len(config["periods"]) != len(config["opt_allele_ranges"]):
		common.WARNING("Error: a different number of period and optimal allele "
					   "ranges specified.")
		return None
	if len(config["periods"]) != len(config["log10_mut_slopes"]):
		common.WARNING("Error: a different number of period and log10_mut_slopes "
					   "specified.")
		return None
	if len(config["periods"]) != len(config["betas"]):
		common.WARNING("Error: a different number of period and betas "
					   "specified.")
		return None
	if len(config["periods"]) != len(config["rhos"]):
		common.WARNING("Error: a different number of period and rhos "
					   "specified.")
		return None
	if len(config["periods"]) != len(config["baseline_mus"]):
		common.WARNING("Error: a different number of period and baseline_mus "
					   "specified.")
		return None
	if len(config["periods"]) != len(config["baseline_mu_alleles"]):
		common.WARNING("Error: a different number of period and baseline_mu_alleles "
					   "specified.")
		return None
	# Check of all or none for species/demog_model/popid specified
	num_stdpopsim_args = sum([config["species"] is None, config["demog_model"] is None,
		config["popid"] is None])
	if num_stdpopsim_args > 0 and num_stdpopsim_args != 3:
		common.WARNING("Error: if using a stdpopsim model must specify "
			           "all of --species, --demographic-model, --popid")
		return None
	if config["species"] is not None and \
			config["species"] not in stdpopsim.species.registered_species.keys():
		common.WARNING("Error: could not find stdpopsim model for "+config["species"])
		return None
	if config["demog_model"] is not None and \
			config["demog_model"] not in \
			[item.id for item in stdpopsim.species.registered_species[config["species"]].demographic_models]:
		common.WARNING("Error: could not find stdpopsim model with ID " + config["demog_model"])
		return None
	if config["popid"] is not None and \
			config["popid"] not in \
			[item.name for item in stdpopsim.species.registered_species[config["species"]].get_demographic_model(config["demog_model"]).populations]:
		common.WARNING("Error: could not find pop " + config["popid"] + " in the specified model")
		return None

	# Additional checks on command line arguments
	if args.periods is not None:
		if args.opt_allele_ranges is None:
			common.WARNING("Error: if you change --periods you must also set --opt-allele-ranges")
			return None
		if args.log10_mut_slopes is None:
			common.WARNING("Error: if you change --periods you must also set --log10-mut-slopes")
			return None
		if args.betas is None:
			common.WARNING("Error: if you change --periods you must also set --betas")
			return None
		if args.rhos is None:
			common.WARNING("Error: if you change --periods you must also set --rhos")
			return None
		if args.baseline_mus is None:
			common.WARNING("Error: if you change --periods you must also set --baseline-mus")
			return None
		if args.baseline_mus is None:
			common.WARNING("Error: if you change --periods you must also set --baseline-mu-alleles")
			return None
	for period in config["periods"]:
		if period < 0:
			common.WARNING("Error: cannot have a period < 0")
			return None
		if period > 6:
			common.WARNING("Error: only rpt. units of <=6 are supported")
			return None
	for ar in config["opt_allele_ranges"]:
		if len(ar) != 2:
			common.WARNING("Improperly formatted allele range {ar}.".format(ar=ar))
			return None
		if type(ar[0]) != int or type(ar[1]) != int:
			common.WARNING("Invalid allele range {ar}. Values must be integers".format(ar=ar))
			return None		
		minval, maxval = ar
		if maxval < minval:
			common.WARNING("Invalid allele range {ar}. Maxval < minval".format(ar=ar))
			return None
		if maxval <= 0 or minval <= 0:
			common.WARNING("Invalid allele range {ar}. min/max must be > 0".format(ar=ar))
			return None
	for beta in config["betas"]:
		if beta < 0 or beta > 1:
			common.WARNING("Error: beta values must be between 0 and 1")
			return None
	for rho in config["rhos"]:
		if rho < 0 or rho > 1:
			common.WARNING("Error: rho values must be between 0 and 1")
			return None
	for bmu in config["baseline_mus"]:
		if bmu < 0 or bmu > 1:
			common.WARNING("Error: mu values must be between 0 and 1")
			return None	
	for bmua in config["baseline_mu_alleles"]:
		if bmua <= 0 :
			common.WARNING("Error: baseline mu allele must be >0")
			return None
	if config["n_effective"] <= 0:
		common.WARNING("Error: --n-effective must be > 0")
		return None
	if config["num_gens"] <= 0:
		common.WARNING("Error: --num-gens must be > 0")
		return None
	if config["num_alleles"] <= 0:
		common.WARNING("Error: --num-alleles must be > 0")
		return None
	if config["abc_num_sims"] <= 0:
		common.WARNING("Error: --abc-num-sims must be > 0")
		return None
	if config["min_mu"] < 0 or config["min_mu"] > 1:
		common.WARNING("Error: --min-mu must be between 0 and 1")
		return None
	if config["max_mu"] < 0 or config["max_mu"] > 1:
		common.WARNING("Error: --max-mu must be between 0 and 1")
		return None
	if config["max_mu"] < config["min_mu"]:
		common.WARNING("Error: --min-mu cannot be larger than --max-mu")
		return None
	if config["end_samp_n"] <= 0:
		common.WARNING("Error: --end-samp-n must be > 0")
		return None
	if config["lrt_num_sims"] <= 0:
		common.WARNING("Error: --lrt-num-sims must be > 0")
		return None
	return config
