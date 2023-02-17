#!/usr/bin/env python

# To test: ./clumpSTR.py --summstats-snps tests/eur_gwas_pvalue_chr19.LDL.glm.linear --clump-snp-field ID --clump-field p-value --clump-chrom-field CHROM --clump-pos-field position --clump-p1 0.2 --out test.clump
import click
import scipy.stats
import sys

class Variant:

    """
    Initialize a clumpSTR object

    Parameters
    ----------
    varid: Variants ID
        Variant's identification (Ej. rs1234)
    chrom: Chromosome
        Variant's chrosmosome identification (Ej. 1:22)
    pos: Position
        Variant's Position (Ej. 12345634)
    pval: P value from the variant
        Variant's P-value (Ej. 0.05)
    vartype: Variant Type
        Type of Variant (SNP or STR)
    """

    def __init__(self, varid, chrom, pos, pval, vartype):
        self.varid = varid
        self.chrom = chrom
        self.pos = pos
        self.pval = pval
        self.vartype = vartype

    def __str__(self):
        return "%s %s %s %s %s"%(self.varid, self.chrom, self.pos, self.pval, self.vartype)

class SummaryStats:
    """
    Initialize summary statistics 

    Parameters
    ----------

    summstats: Initialize Summary Statistics with constructor 
        keep track of summary statistics

    """

    def __init__(self):
        self.summstats = []

    def Load(self, statsfile, vartype="SNP", pthresh=1.0,
            snp_field="SNP", p_field="P",
            chrom_field="CHR", pos_field="POS"):
        """
        Load summary statistics
        Ignore variants with pval < pthresh

        Parameters
        ----------
        self: constructor
            load function construtor
        statsfile: summary statistical file
            Summary statistic file to be input by the user 
        vartype: variant Type
            Type of Variant (SNP or STR)
        pthresh: p-value treshold
            P value input as treshold for loading data (Ej. 0.05)
        snp_field: SNP column
            Location of the SNP in the summary statistics
        p_field: P-value column
            Location of the P-value of the snps in the summary statistics
        chrom_field: chromosome column
            Location of the chromosome of the snps in the summary statistics
        pos_field: position column
            Location of the position of the snps in the summary statistics
        """
        
        summstats = [] # List of Variants

        # First, parse header line to get col. numbers
        f = open(statsfile, "r")
        header_items = [item.strip() for item in f.readline().split()]
        try:
            snp_col = header_items.index(snp_field)
        except ValueError:
            print("Could not find %s in header"%snp_field)
            sys.exit(1)
        try:
            p_col = header_items.index(p_field)
        except ValueError:
            print("Could not find %s in header"%p_field)
            sys.exit(1)
        try:
            chrom_col = header_items.index(chrom_field)
        except ValueError:
            print("Could not find %s in header"%chrom_field)
            sys.exit(1)
        try:
            pos_col = header_items.index(pos_field)
        except ValueError:
            print("Could not find %s in header"%pos_field)
            sys.exit(1)

        # Now, load in stats. Skip things with pval>pthresh
        line = f.readline()
        while line.strip() != "":
            items = [item.strip() for item in line.strip().split()]
            if float(items[p_col]) > pthresh:
                line = f.readline()
                continue
            summstats.append(Variant(items[snp_col], items[chrom_col],
                int(items[pos_col]), float(items[p_col]), vartype))
            line = f.readline()
        f.close()
        self.summstats.extend(summstats)

    def GetNextIndexVariant(self, index_pval_thresh):
        """
        Get the next index variant, which is the 
        variant with the best p-value

        If no more variants below the clump-p1 threshold,
        return None

        Not yet implemented
        """
        best_var = None
        best_var_p = 1.0
        for variant in self.summstats:
            if variant.pval < best_var_p and variant.pval<index_pval_thresh:
                best_var = variant
                best_var_p = variant.pval
        return best_var

    def QueryWindow(self, indexvar, window_kb):
        """
        Find all candidate variants in the specified
        window around the index variant

        Not yet implemented
        """
        # First get stats on the indexvariant
        chrom = indexvar.chrom
        pos = indexvar.pos

        # Find candidates in the window
        candidates = []
        for variant in self.summstats:
            if variant.chrom == chrom and abs(variant.pos-pos)/1000 < window_kb:
                candidates.append(variant)

        return candidates

    def RemoveClump(self, clumpvars):
        """
        Remove the variants from a clump 
        from further consideration
        """
        keepvars = []
        for variant in self.summstats:
            if variant not in clumpvars:
                keepvars.append(variant)
        self.summstats = keepvars

def GetOverlappingSamples(snpgts, strgts):
    """
    Get overlapping set of samples
    """
    return [] # TODO

def LoadVariant(var, snpgts, strgts, samples):
    """
    Extract vector of genotypes for this variant
    """
    return [] # TODO

def ComputeLD(var1, var2, snpgts, strgts, samples):
    """
    Compute the LD between two variants
    """
    # Load genotypes
    gts1 = LoadVariant(var1, snpgts, strgts, samples)
    gts2 = LoadVariant(var2, snpgts, strgts, samples)
    # TODO - possibly check for NAs and remove them
    # Compute and return Pearson r2
    return scipy.stats.pearsonr(gts1, gts2)[0]**2

def WriteClump(indexvar, clumped_vars, outf):
    """
    Write a clump to the output file
    Not yet implemented
    """
    outf.write("\t".join([indexvar.varid, indexvar.chrom, str(indexvar.pos),
        str(indexvar.pval), indexvar.vartype, 
        ",".join([str(item) for item in clumped_vars])])+"\n")

@click.command()
@click.option('--summstats-snps', type=str, help='File to load snps summary statistics')
@click.option('--summstats-strs', type=str, help='File to load strs summary statistics')
@click.option('--gts-snps', type=str, help='SNP genotypes (VCF or PGEN)')
@click.option('--gts-strs', type=str, help='STR genotypes (VCF)')
@click.option('--clump-p1', type=float, default=0.0001, help='Max pval to start a new clump')
@click.option('--clump-p2', type=float, default=0.01, help='Filter for pvalue less than')
@click.option('--clump-snp-field', type=str, default='SNP', help='Column header of the variant ID')
@click.option('--clump-field', type=str, default='P', help='Column header of the p-values')
@click.option('--clump-chrom-field', type=str, default='CHR', help='Column header of the chromosome')
@click.option('--clump-pos-field', type=str, default='POS', help='Column header of the position')
@click.option('--clump-kb', type=float, default=250, help='clump kb radius')
@click.option('--clump-r2', type=float, default=0.5, help='r^2 threshold')
@click.option('--out', type=str, required=True, help='Output filename')
def clumpstr(summstats_snps, summstats_strs, gts_snps, gts_strs, clump_p1, clump_p2,
    clump_snp_field, clump_field, clump_chrom_field, clump_pos_field,
    clump_kb, clump_r2, out):
    ###### User checks ##########
    # TODO - need one of summstats_snps or summstats_strs
    # TODO - if summstats_snps, also need gts_snps
    # TODO - if summstats_strs, also need gts_strs

    ###### Load summary stats ##########
    summstats = SummaryStats()
    if summstats_snps is not None:
        summstats.Load(summstats_snps, vartype="SNP", pthresh=clump_p2,
            snp_field=clump_snp_field, p_field=clump_field,
            chrom_field=clump_chrom_field, pos_field=clump_pos_field)
    if summstats_strs is not None:
        summstats.Load(summstats_strs, vartype="STR", pthresh=clump_p2,
            snp_field=clump_snp_field, p_field=clump_field,
            chrom_field=clump_chrom_field, pos_field=clump_pos_field)

    ###### Set up genotypes ##########
    snpgts = None
    strgts = None
    if gts_snps is not None:
        pass # TODO - load with haptools data
    if gts_strs is not None:
        pass # TODO - load with haptools? TRTools?
    samples = GetOverlappingSamples(snpgts, strgts)

    ###### Setup output file ##########
    outf = open(out, "w")
    outf.write("\t".join(["ID","CHROM","POS","P","VARTYPE","CLUMPVARS"])+"\n")

    ###### Perform clumping ##########
    indexvar = summstats.GetNextIndexVariant(clump_p1)
    while indexvar is not None:
        candidates = summstats.QueryWindow(indexvar, clump_kb)
        clumpvars = []
        for c in candidates:
            r2 = ComputeLD(c, indexvar, snpgts, strgts, samples)
            if r2 > clump_r2:
                clumpvars.append(c)
        WriteClump(indexvar, clumpvars, outf)
        summstats.RemoveClump(clumpvars+[indexvar])
        indexvar = summstats.GetNextIndexVariant(clump_p1)
    outf.close()

if __name__ == '__main__':
    clumpstr()