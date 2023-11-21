#!/usr/bin/env python3
"""
Simulate NGS reads using ART, while capturing
STR-specific stutter errors
"""

import argparse
import numpy as np
import os
import pyfaidx
import random
import re
import shutil
import subprocess
import sys
import tempfile
import trtools.utils.common as common
import trtools.prancSTR as ms
import trtools.utils.utils as utils
from trtools import __version__

_MAXWINDOW = 1000000

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

	If we encounter an error parsing, then
	chrom, start, end are None
	"""
	if type(coords) != str:
		return None, None, None
	if re.match(r"\w+:\d+-\d+", coords) is None:
		return None, None, None
	chrom = coords.split(":")[0]
	start = int(coords.split(":")[1].split("-")[0])
	end = int(coords.split(":")[1].split("-")[1])
	if start >= end:
		common.WARNING("Problem parsing coordinates {}. start>=end".format(coords))
		return None, None, None
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
	   Return 0 if no such delta exists, which
	   can happen e.g. with low rho
	"""
	delta = np.ceil(np.log(pthresh/(sprob*rho))/np.log(1-rho)+1)
	if delta < 1: return 0
	return int(delta)

def GetTempDir(debug=False, dir=None):
	r"""
	Create a temporary directory to store
	intermediate fastas and fastqs

	Parameters
	----------
	debug : bool
	   Ignored for now
	dir : str
	   Directory in which to create the temporary directory

	Returns
	-------
	dirname : str
	   Path to the temporary directory
	   Return None if there was a problem creating the directory
	"""
	if not os.path.isdir(dir):
		common.WARNING("Error: The specified tmpdir {} does"
			" not exist".format(dir))
		return None
	dirname = tempfile.mkdtemp(dir=dir)
	return dirname

def GetAlleleSeq(seq_preflank, seq_postflank, \
				seq_repeat, repeat_unit, delta):
	r"""
	Generate a new allele with a change of 
	delta repeat units

	Parameters
	----------
	seq_preflank : str
	   Sequence upstream of the STR
	seq_postflank : str
	   Sequence downstream of the STR
	seq_repeat : str
	   Sequence of the STR region
	repeat_unit : str
	   Repeat unit sequence
	delta : int
	   Change in repeat units compared to ref
	tmpdir : str
	   Path to create the fasta in

	Returns
	-------
	newseq : str
	   New repeat allele sequence
	   Return None if there was a problem
	"""
	newseq = seq_preflank
	if delta == 0:
		newseq += seq_repeat
	elif delta > 0:
		newseq += seq_repeat + repeat_unit*delta
	else:
		subtract_size = -1*delta*len(repeat_unit)
		if subtract_size > len(seq_repeat):
			common.WARNING("Error: tried to delete {} bp but the "
				"total repeat is {} bp long".format(subtract_size, len(seq_repeat)))
			return None
		newseq += seq_repeat[:-1*subtract_size]
	newseq += seq_postflank
	return newseq

def CreateAlleleFasta(newseq, delta, tmpdir):
	r"""
	Create fasta file for this allele
	Return the path to the fasta

	Parameters
	----------
	newseq : str
	   New repeat allele sequence
	delta : int
	   Change in repeat units compared to ref
	tmpdir : str
	   Path to create the fasta in

	Returns
	-------
	fname : str
	   Path to created fasta file
	"""
	fname = os.path.join(tmpdir, "simTR_{}.fa".format(delta))
	with open(fname, "w") as f:
		f.write(">seq_{}\n".format(delta))
		f.write(newseq+"\n")
	return fname

def SimulateReads(newfasta, coverage, read_length,
		single, insert, sd, tmpdir, delta, art_cmd):
	r"""
	Run ART on our dummy fasta file
	with specified parameters

	Parameters
	----------
	newfasta : str
	   Path to dummy fasta file
	coverage : int
	   Desired coverage level (ART -f)
	read_length : int
	   Read length (ART -l)
	single : bool
	   Use single-end read mode
	insert : float
	   Mean fragment length (ART -m)
	sd : float
	   Std dev of fragment length distribution (ART -s)
	tmpdir : str
	   Path to create the fasta in
	delta : int
	   Difference in repeat units from reference
	   Used for naming files
	art_cmd : str
	   Command to run ART

	Returns
	-------
	fq1file, fq2file : str, str
	   Paths to fastq file output
	   for the two read pairs.
	   Return None, None if failed
	   If single end mode, fq2file is None
	"""
	outprefix = os.path.join(tmpdir, "artsim_{}_".format(delta))
	cmd = [art_cmd, \
			"-i", newfasta, \
			"-l", str(read_length), \
			"-f", str(coverage), \
			"-m", str(insert), \
			"-s", str(sd), \
			"-o", outprefix
		]
	if not single: cmd.append("-p")
	process = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
	if process.returncode != 0:
		common.WARNING(process.stdout)
		return None, None
	fq1file = outprefix+"1.fq"
	if single:
		fq2file = None
	else: fq2file = outprefix+"2.fq"
	return fq1file, fq2file

def WriteCombinedFastqs(fqfiles, fname):
	r"""
	Concatenate fastq files to output

	Parameters
	----------
	fqfiles : list of str
	   List of paths to fastqfiles to concatenate
	fname : str
	   Name of final output file
	"""
	with open(fname, "w") as outfile:
		for fqn in fqfiles:
			with open(fqn) as infile:
				for line in infile:
					outfile.write(line)
	return

def main(args):
	if not os.path.exists(args.ref):
		common.WARNING("Error: {} does not exist".format(args.ref))
		return 1
	if args.u < 0 or args.u > 1:
		common.WARNING("Error: --u u ({}) is not between 0 and 1".format(args.u))
		return 1
	if args.d < 0 or args.d > 1:
		common.WARNING("Error: --d ({}) is not between 0 and 1".format(args.d))
		return 1
	if (args.d + args.u) > 1:
		common.WARNING("Error: --d ({}) and --u ({}) can't add to more than 1".format(args.d, args.u))
		return 1
	if args.rho < 0 or args.rho > 1:
		common.WARNING("Error: --rho ({}) is not between 0 and 1".format(args.rho))
		return 1
	if args.p_thresh < 0 or args.p_thresh > 1:
		common.WARNING("Error: --p_thresh ({}) is not between 0 and 1".format(args.p_thresh))
		return 1
	if args.coverage < 0:
		common.WARNING("Error: --coverage ({}) cannot be less than 0".format(args.coverage))
		return 1
	if args.read_length < 0:
		common.WARNING("Error: --read_length ({}) cannot be less than 0".format(args.read_length))
		return 1
	if args.read_length > args.insert:
		common.WARNING("Error: --read_length ({}) must be shorter than"
			" --insert ({})".format(args.read_length, args.insert))
	if args.insert < 0:
		common.WARNING("Error: --insert ({}) cannot be less than 0".format(args.insert))
		return 1
	if args.sd < 0:
		common.WARNING("Error: --sd ({}) cannot be less than 0".format(args.sd))
		return 1
	if args.window < 0:
		common.WARNING("Error: --window ({}) cannot be less than 0".format(args.window))
		return 1
	if args.window > _MAXWINDOW:
		common.WARNING("Error: --window ({}) must be <= {}".format(args.window, _MAXWINDOW))
		return 1
	if args.window < args.insert:
		common.WARNING("Error: --window ({}) must be greater than the fragment length".format(args.window))
		return 1
	if not os.path.exists(os.path.dirname(os.path.abspath(args.outprefix))):
		common.WARNING("Error: The directory which contains the output location {} does"
			" not exist".format(args.outprefix))
		return 1
	if args.seed is not None:
		random.seed(args.seed)
	art_path = None
	if args.art is not None:
		if not os.path.exists(args.art) and not shutil.which(args.art):
			common.WARNING("Error: ART path {} does not exist".format(args.art))
			return 1
		else: art_path = args.art
	else:
		if shutil.which("art_illumina") is None:
			common.WARNING("Error: Could not find art_illumina executable")
			return 1
		else: art_path = "art_illumina"
	common.MSG("Using this command for ART: {}".format(art_path), debug=args.debug)
	# Parse coordinates
	chrom, start, end = ParseCoordinates(args.coords)
	if chrom is None:
		common.WARNING("Error: could not extract coordinates")
		return 1

	# Determine range of deltas to consider
	highdelta = GetMaxDelta(args.u, args.rho, args.p_thresh)
	lowdelta = GetMaxDelta(args.d, args.rho, args.p_thresh)

	# Extract ref sequences
	refgenome = pyfaidx.Fasta(args.ref)
	if chrom not in refgenome.records:
		common.WARNING("Could not find {} in {}".format(chrom, args.ref))
		return 1
	seq_repeat = str(refgenome[chrom][start-1:end]).upper()
	seq_preflank = str(refgenome[chrom][start-args.window-1:start-1]).upper()
	seq_postflank = str(refgenome[chrom][end:end+args.window]).upper()

	# Require the whole region to be at least as long as window
	seq_len = len(seq_preflank + seq_repeat + seq_postflank)
	if seq_len <= args.window:
		common.WARNING("Extracted sequence length shorter {} than window {}".format(seq_len, args.window))
		return 1

	check_rpt = utils.LongestPerfectRepeat(seq_repeat, args.repeat_unit, check_reverse=False)
	if check_rpt <= len(args.repeat_unit)*2:
		common.WARNING("Did not find the unit {} a sufficient "
			"number of times in the repeat region {}".format(args.repeat_unit, seq_repeat))
		return 1
	else:
		common.MSG("Found a {} bp stretch with a perfect match to the repeat unit".format(check_rpt), \
			debug=args.debug)

	# Create folder structure
	tmpdir = GetTempDir(debug=args.debug, dir=args.tmpdir)
	if tmpdir is None:
		common.WARNING("ERROR: could not create temoporary directory")
		return 1
	common.MSG("Created temporary directory at {}".format(tmpdir), debug=args.debug)

	# Simulate reads from each potential allele
	fq1files = []
	fq2files = []

	for delta in range(-1*lowdelta, highdelta+1):
		sprob = ms.StutterProb(delta, args.u, args.d, args.rho)
		cov = np.random.binomial(args.coverage, sprob)
		newseq = GetAlleleSeq(seq_preflank, seq_postflank, seq_repeat, \
				args.repeat_unit, delta)
		if newseq is None:
			common.WARNING("Problem getting allele sequence for delta={}".format(delta))
			return 1
		newfasta = CreateAlleleFasta(newseq, delta, tmpdir)
		fq1, fq2 = SimulateReads(newfasta, cov,
			args.read_length, args.single, args.insert, args.sd, 
			tmpdir, delta, art_path)
		if fq1 is None:
			return 1
		if args.single:
			common.MSG("Created {}".format(fq1), debug=args.debug)
		else: common.MSG("Created {} and {}".format(fq1, fq2), debug=args.debug)
		fq1files.append(fq1)
		fq2files.append(fq2)

	# Combine all fastqs to single output
	WriteCombinedFastqs(fq1files, args.outprefix+"_1.fq")
	common.MSG("Output fastq file {}".format(args.outprefix+"_1.fq", debug=args.debug))

	if not args.single:
		WriteCombinedFastqs(fq2files, args.outprefix+"_2.fq")
		common.MSG("Output fastq file {}".format(args.outprefix+"_2.fq", debug=args.debug))

	return 0

def getargs(): # pragma: no cover
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
	inout_group.add_argument("--tmpdir", help="Temporary directory to store intermediate "
		"results. Default: {}".format(os.environ.get("TMPDIR","/tmp")), type=str, \
			default=os.environ.get("TMPDIR","/tmp"))
	stutter_group = parser.add_argument_group("Stutter simulation parameters")
	stutter_group.add_argument("--u", help="Probability of adding additional copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--d", help="Probability of deleting copy of repeat", type=float, default=0.05)
	stutter_group.add_argument("--rho", help="Size of stutter-induced changes", type=float, default=0.9)
	stutter_group.add_argument("--p-thresh", help="Ignore stutter alleles expected to have lower than this frequency", \
		type=float, default=0.001)
	stutter_group.add_argument("--seed", help="Set the seed to make runs reproducible", type=int)
	seq_group = parser.add_argument_group("Sequencing parameters")
	seq_group.add_argument("--coverage", help="Target coverage level", type=int, default=1000)
	seq_group.add_argument("--read-length", help="Length of each read (bp)", type=int, default=100)
	seq_group.add_argument("--insert", help="Mean fragment length", type=int, default=350)
	seq_group.add_argument("--sd", help="Std. deviation of fragment length", type=int, default=50)
	seq_group.add_argument("--window", help="Size of window around target TR to sequence (bp)", type=int, default=1000)
	seq_group.add_argument("--single", help="Generate single-end reads (default is paired)", action="store_true")
	other_group = parser.add_argument_group("Other options")
	other_group.add_argument("--art", help="Path to ART simulator package (Default: art_illumina)", \
		type=str, required=False)
	other_group.add_argument("--debug", help="Run in debug mode", action="store_true")
	ver_group = parser.add_argument_group("Version")
	ver_group.add_argument("--version", action="version", \
		version='{version}'.format(version=__version__))
	args = parser.parse_args()
	return args

def run(): # pragma: no cover
    args = getargs()
    if args == None:
        sys.exit(1)
    else:
        retcode = main(args)
        sys.exit(retcode)

if __name__ == "__main__": # pragma: no cover
    run()
