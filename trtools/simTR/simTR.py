#!/usr/bin/env python3
"""
Simulate NGS reads using ART, while capturing
STR-specific stutter errors
"""

import argparse
import os
import pyfaidx
import shutil
import sys
import trtools.utils.utils as utils
from trtools import __version__

def ParseCoordinates(coords):
	r"""
	Extract chrom, start, end from
	coordinate string

	Parameters
	----------
	coords : str
	   Coordinate string in the form
	   chrom:start-end

	Returns
	-------
	chrom : str
	   Chromosome name
	start : int
	   start coordinate
	end : int
	   end coordinate
	"""
	chrom = coords.split(":")[0]
	start = int(coords.split(":")[1].split("-")[0])
	end = int(coords.split(":")[1].split("-")[1])
	return chrom, start, end

def GetMaxDelta(sprob, rho, pthresh):
	r"""
	Compute the max delta for which the frequency would
	be great than pthresh

	based on freq = sprob*rho*(1-rho)**(delta-1)

	Parameters
	----------
	sprob : float
	   Stutter probability
	rho : float
	   Stutter step size parameters
	pthresh : float
	   Minimum frequency threshold

	Returns
	-------
	delta : int
	   Highest delta for which freq>prob
	"""
	delta = np.ceil(np.log(pthresh/(sprob*rho))/np.log(1-rho)+1)
	return delta

def CheckRepeatUnit(seq, repeat_unit):
	r"""
	Count the number of consecutive times
	the specified repeat unit is in the sequence
	"""
	# TODO - check if we have this function somewhere else
	return 0

def GetTempDir():
	r"""
	Create a temporary directory to store
	intermediate fastas and fastqs
	"""
	# TODO print out name of dir
	# TODO if debug mode, don't delete it (add debug mode)
	return "TODO" # TODO

def GetStutterProb(delta, u, d, rho):
	r"""
	Return the probability of seeing the given delta
	"""
	# TODO - don't we have this in mosaicSTR? can copy
	return 0

def CreateAlleleFasta(seq_preflank, seq_postflank, \
				seq_repeat, repeat_unit, delta, tmpdir):
	r"""
	Create fasta file for this allele
	Return the path to the fasta
	"""
	pass # TODO

def SimulateReads(newfasta, coverage, read_length,
		insert, sd, outprefix):
	r"""
	TODO
	"""
	pass # TODO

def WriteCombinedFastqs(fqfiles, fname):
	r"""
	concatenate fastq files to output
	"""
	pass # TODO
	
def main(args):
	if not os.path.exists(args.ref):
		common.WARNING("Error: {} does not exist".format(args.ref))
		return 1
	if args.art is not None:
		if not os.path.exists(args.art):
			common.WARNING("Error: ART path {} does not exist".format(args.art))
			return 1
	else:
		if shutil.which("art_illumina") is None:
			common.WARNING("Error: Could not find art_illumina executable")
			return 1
	if args.u < 0 or args.u > 1:
		common.WARNING("Error: --u u ({}) is not between 0 and 1".format(args.u))
		return 1
	if args.d < 0 or args.d > 1:
		common.WARNING("Error: --d ({}) is not between 0 and 1".format(args.d))
		return 1
	if args.rho < 0 or args.rho > 1:
		common.WARNING("Error: --rho ({}}) is not between 0 and 1".format(args.rho))
		return 1
	if args.p_thresh < 0 or args.p_thresh > 1:
		common.WARNING("Error: --p_thresh ({}) is not between 0 and 1".format(args.p_thresh))
		return 1
	if args.coverage < 0:
		common.WARNING("Error: --coverage ({}}) cannot be less than 0".format(args.coverage))
		return 1
	if args.read_length < 0:
		common.WARNING("Error: --read_length ({}) cannot be less than 0".format(args.read_length))
		return 1
	if args.insert < 0:
		common.WARNING("Error: --insert ({}) cannot be less than 0".format(args.insert))
		return 1
	if args.sd < 0:
		common.WARNING("Error: --sd ({}) cannot be less than 0".format(args.sd))
		return 1
	if args.window < 0:
		common.WARNING("Error: --window ({}) cannot be less than 0".format(args.window))
		return 1
	if not os.path.exists(os.path.dirname(os.path.abspath(args.outprefix))):
		common.WARNING("Error: The directory which contains the output location {} does"
			" not exist".format(args.outprefix))
		return 1

	# Parse coordinates
	chrom, start, end = ParseCoordinates(args.coords)

	# Determine range of deltas to consider
	highdelta = GetMaxDelta(args.u, args.rho, args.p_thresh)
	lowdelta = GetMaxDelta(args.d, args.rho, args.p_thresh)

	# Extract ref sequences
	refgenome = pyfaidx.Fasta(args.ref)
	seq_repeat = refgenome[chrom][start-1:end]
	seq_preflank = refgenome[chrom][start-args.window-1:start-1]
	seq_postflank = refgenome[chrom][end:end+window]
	check_rpt = CheckRepeatUnit(seq_repeat, args.repeat_unit)
	if check_rpt == 0:
		common.WARNING("Did not find the unit {} in the repeat region {}".format(repeat_unit, seq_repeat))
		return 1
	else:
		common.MSG("Found the repeat unit {} times in the repeat region".format(check_rpt))

	# Create folder structure
	tmpdir = GetTempDir()

	# Simulate reads from each potential allele
	fq1files = []
	fq2files = []
	for delta in range(-1*lowdelta, highdelta+1):
		sprob = GetStutterProb(delta, args.u, args.d, args.rho)
		newfasta = CreateAlleleFasta(seq_preflank, seq_postflank, \
				seq_repeat, args.repeat_unit, delta, tmpdir)
		fq1, fq2 = SimulateReads(newfasta, coverage=int*sprob*args.coverage,
			read_length=args.read_length, insert=args.insert, sd=args.sd, args.outprefix)
		fq1files.append(fq1)
		fq2files.append(fq2)

	# Combine all fastqs to single output
	WriteCombinedFastqs(fq1files, args.outprefix+"_1.fq")
	WriteCombinedFastqs(fq2files, args.outprefix+"_1.fq")

def getargs():
	parser = argparse.ArgumentParser(
		__doc__,
		formatter_class=utils.ArgumentDefaultsHelpFormatter
	)
	inout_group = parser.add_argument_group("Input/output")
	inout_group.add_argument("--ref", help="Path to reference genome", type=str, required=True)
	inout_group.add_argument("--coords", help="Coordinates for the target TR (chrom:start-end)", \
		type=str, required=True)
	inout_group.add_argument("--repeat-unit", help="Repeat unit of the target TR", \
		type=str, required=True)
	inout_group.add_argument("--outprefix", help="Prefix to name output files", type=str, required=True)
	stutter_group = parser.add_argument_group("Stutter simulation parameters")
	stutter_group.add_argument("--u", help="Probability of adding additional copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--d", help="Probability of deleting copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--rho", help="Size of stutter-induced changes", type=float, default=0.9)
	stutter_group.add_argument("--p_thresh", help="Ignore stutter alleles expected to have lower than this frequency", \
		type=float, default=0.01)
	seq_group = parser.add_argument_group("Sequencing parameters")
	seq_group.add_argument("--coverage", help="Target coverage level", type=int, default=1000)
	seq_group.add_argument("--read_length", help="Length of each read (bp)", type=int, default=100)
	seq_group.add_argument("--insert", help="Mean fragment length", type=int, default=350)
	seq_group.add_argument("--sd", help="Std. deviation of fragmen tlength", type=int, default=50)
	seq_group.add_argument("--window", help="Size of window around target TR to sequence (bp)", type=int, default=1000)
	other_group = parser.add_argument_group("Other options")
	other_group.add_argument("--art", help="Path to ART simulator package (Default: art_illumina)", \
		type=str, required=True)
	ver_group = parser.add_argument_group("Version")
	ver_group.add_argument("--version", action="version", \
		version='{version}'.format(version=__version__))
	args = parser.parse_args()
	return args

def run():  # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__":
    run()