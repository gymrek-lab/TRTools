#!/usr/bin/env python
"""

Author: Shubham Saini
shubhamsaini@eng.ucsd.edu


Plink style LD-based clumping of GWAS results from plinkSTR
Example:
./clump.py --clump-p1 0.2 --clump-p2 0.2 --clump-r2 0.1 --clump-kb 3000 \
--assoc /storage/s1saini/gwas_simu/results/assoc_output.100.assoc \
--vcf /storage/s1saini/hipstr_allfilters/phased_feb18/hipstr.chr21.phased.vcf.gz

"""


import argparse
import sys
import numpy as np
import pandas as pd
from cyvcf2 import VCF

def PrintLine(text, f):
    f.write(text+"\n")
    f.flush()

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--clump-p1", help="Significance threshold for index SNPs", required=True, type=str)
    parser.add_argument("--clump-p2", help="Secondary significance threshold for clumped SNPs", required=True, type=str)
    parser.add_argument("--clump-r2", help="LD threshold for clumping", required=True, type=str)
    parser.add_argument("--clump-kb", help="Physical distance threshold for clumping (kb)", required=True, type=str)
    parser.add_argument("--assoc", help="Association Results from plinkSTR", required=True, type=str)
    parser.add_argument("--vcf", help="VCF File", required=True, type=str)
    parser.add_argument("--out", help="Output file name. Default: stdout", required=False, type=str)
    parser.add_argument("--region", help="Region", required=False, type=str)
    args = parser.parse_args()

    CLUMP_P1 = float(args.clump_p1)
    CLUMP_P2 = float(args.clump_p2)
    CLUMP_R2 = float(args.clump_r2)
    CLUMP_KB = int(args.clump_kb)*1000

    ASSOC_FILE = args.assoc
    VCF_FILE = args.vcf

    # Prepare output
    if args.out: outf = open(args.out,"w")
    else: outf = sys.stdout

    PrintLine("CHR SNP BP P TOTAL SP2", outf)

    assoc_results = pd.read_csv(ASSOC_FILE, delim_whitespace=True)
    assoc_results['keep'] = 0
    assoc_results['P'] = pd.to_numeric(assoc_results['P'],errors='coerce')
    assoc_results.loc[assoc_results['P'] < CLUMP_P2,'keep'] = 2
    assoc_results.loc[assoc_results['P'] < CLUMP_P1,'keep'] = 1
    assoc_results = assoc_results[assoc_results['keep']>0]
    assoc_results['keep'] = -1

    if args.region:
        region_chrom = int(args.region.split(":")[0])
        assoc_results = assoc_results[assoc_results['CHROM'] == region_chrom]

        ### implement position

    while(assoc_results[(assoc_results['keep']==-1)].shape[0]>0):
        strs_retain = []
        strs_discard = []

        index_str_pmin = np.min(assoc_results[assoc_results['keep']==-1]['P'])
        index_str_id = assoc_results[(assoc_results['P']==index_str_pmin) & (assoc_results['keep']==-1)]['SNP'].values[0]
        index_str_bp = int(assoc_results[assoc_results['SNP']==index_str_id]['BP'])
        chrom = int(assoc_results[assoc_results['SNP']==index_str_id]['CHROM'])
        assoc_results.loc[assoc_results['SNP']==index_str_id,'keep'] = 1

        start = int(index_str_bp - CLUMP_KB)
        end = int(index_str_bp + CLUMP_KB)
        curr_assoc_results = assoc_results[(assoc_results['BP']>=start) & (assoc_results['BP']<=end) & (assoc_results['keep']==-1)]
        if curr_assoc_results.shape[0] == 0:
            #PrintLine('No STRs to clump with the index STR %d:%d'%(chrom,index_str_bp), outf)
            PrintLine('%d %s %d %s %d NONE'%(chrom, index_str_id, index_str_bp, str(index_str_pmin), len(strs_discard)), outf)
            continue;

        assoc_ids = curr_assoc_results['SNP'].values.tolist()
        vcf = VCF(VCF_FILE)
        index_str_gt = []
        for v in vcf('%d:%d-%d'%(chrom,index_str_bp,index_str_bp)):
            if str(v.ID) != index_str_id:
                #PrintLine('Index STR %s not found in VCF file'%(index_str_id), outf)
                continue
            str_gt = []
            for gt in v.gt_bases:
                str_gt.append(np.sum([len(i)-len(v.REF) for i in gt.split("|")]))
            index_str_gt = str_gt

        vcf = VCF(VCF_FILE)
        vcf_strids = []
        for v in vcf('%d:%d-%d'%(chrom,start,end)):
            if str(v.ID) == index_str_id:
                continue
	    else:
		str_gt = []
                gt_array = []
                if str(v.ID) in assoc_ids:
                    	vcf_strids.append(str(v.ID))
                    	if CLUMP_R2 == 0:
				strs_discard.append(str(v.ID))
			else:
				for gt in v.gt_bases:
                        		str_gt.append(np.sum([len(i)-len(v.REF) for i in gt.split("|")]))
                    		corr_matrix = np.corrcoef(str_gt,index_str_gt)[0,1] ** 2
                    		if corr_matrix < CLUMP_R2:
                        		strs_retain.append(str(v.ID))
                    		else:
                        		strs_discard.append(str(v.ID))


        assoc_results.loc[assoc_results['SNP'].isin(strs_discard),'keep'] = 0
        PrintLine('%d %s %d %s %d %s'%(chrom, index_str_id, index_str_bp, str(index_str_pmin), len(strs_discard), ",".join(strs_discard)), outf)

if __name__ == "__main__":
    main()
