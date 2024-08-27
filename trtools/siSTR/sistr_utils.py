"""
Common utilities for SISTR commands
"""

import json

import trtools.utils.common as common

# Default SISTR params
DEFAULTS = {
	"periods": [2, 3, 4],
	"opt_allele_ranges": [(11, 20), (5, 13), (7, 10)]
}

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
		config["opt_allele_ranges"] = [(int(item.split("-")[0]), int(item.split("-")[1])) \
			for item in args.opt_allele_ranges.strip().split(",")]

	# Checks on inputs
	if len(config["periods"]) != len(config["opt_allele_ranges"]):
		common.WARNING("Error: a different number of period and optimal allele "
					   "ranges specified.")
		return None
	if args.periods is not None:
		if args.opt_allele_ranges is None:
			common.WARNING("Error: if you change --periods you must also set --opt-allele-ranges")
			return None
	for ar in config["opt_allele_ranges"]:
		minval, maxval = ar
		if maxval < minval:
			common.WARNING("Invalid allele range {ar}. Maxval < minval".format(ar=ar))
			return None
		if maxval < 0 or minval < 0:
			common.WARNING("Invalid allele range {ar}. min/max must be > 0".format(ar=ar))
			return None

	return config
