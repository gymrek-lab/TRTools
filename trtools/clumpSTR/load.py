import click
import sys
import pandas as pd

class SummaryStats:
    """
    Keep track of summary statistics

    TODO: add detailed class and methods documentation
    """

    def __init__(self):
	self.snp_summstats = None
	self.str_summstats = None

    def Load(self, statsfile, vartype="SNP", pthresh=1.0):
        """
        Load summary statistics
        Ignore variants with pval < pthresh
        Not yet implemented
        """
   	df = pd.read_csv(statsfile, delimiter='\t') ##optimize this function to work better
   	df = df[df['p-value'] < pthresh]
	if vartype =="SNP":
		self.snp_summstats = df
	elif vartype == "STR"
		self.str_summstats = df
	else:
		print("Invalid vartype")
		sys.exit(1) ##Cleaning exit the program function
   	click.echo(df)

    def GetNextIndexVariant(self, index_pval_thresh):
        """
        Get the next index variant, which is the 
        variant with the best p-value

        If no more variants below the clump-p1 threshold,
        return None

        Not yet implemented
        """
        return None # TODO

    def QueryWindow(self, indexvar, window_kb):
        """
        Find all candidate variants in the specified
        window around the index variant

        Not yet implemented
        """
        return [] # TODO

    def RemoveClump(self, indexvar, clumpvars):
        """
        Remove the variants from a clump 
        from further consideration
        """
        pass # TODO

def ComputeLD(var1, var2):
    """
    Compute the LD between two variants
    """
    return 0 # TODO

def WriteClump(indexvar, clumped_vars):
    """
    Write a clump to the output file
    Not yet implemented
    """
    pass # TODO

@click.command()
@click.option('--summstats_snps', type=click.File('r'), help='File to load snps summary statistics')
@click.option('--summstats_strs', type=click.File('r'), help='File to load strs summary statistics')
@click.option('--clump_p2', type=float, default=0.01, help='Filter for pvalue less than')


if __name__ == '__main__':

    summstats = SummaryStats()
    if summstats_snps is not None:
    	summstats.Load(summstats_snps, vartype="SNP", pthresh=clump_p2)
    if summstats_strs is not None:
    	summstats.Load(summstats_strs, vartype="STR", pthresh=clump_p2)


	#load_file()
	#summstats = SummaryStats()