#!/usr/bin/env python3

import random
import pathlib

import cyvcf2
import numpy as np

import subprocess as sp

random.seed(11)

SCRIPT_DIR = pathlib.Path(__file__).parent.resolve()

# biallelic

vcf = cyvcf2.VCF(str(SCRIPT_DIR / 'many_samples_biallelic.vcf.gz'))
samples = vcf.samples

with open(str(SCRIPT_DIR / 'gp_dosages.tsv'), 'w') as gp_out, open(str(SCRIPT_DIR / 'ap1_dosages.tsv'), 'w') as ap1_out, open(str(SCRIPT_DIR / 'ap2_dosages.tsv'), 'w') as ap2_out:
    for var in vcf:
        gp_out.write('{}\t{}\t{}'.format(var.CHROM, var.POS, var.POS))
        ap1_out.write('{}\t{}\t{}'.format(var.CHROM, var.POS, var.POS))
        ap2_out.write('{}\t{}\t{}'.format(var.CHROM, var.POS, var.POS))
        gts = var.genotype.array()
        for i in range(gts.shape[0]):
            if -1 in gts[i, :]:
                gp_out.write('\t.')
                ap1_out.write('\t.')
                ap2_out.write('\t.')
                continue
            # plink doesn't store dosages with much precision,
            # so reduce the precision of the values we emit
            p1 = round(random.random(), 2)
            p2 = round(random.uniform(0, 1 - p1), 2)
            p3 = max(1 - p1 - p2, 0.0)
            ps = [p1, p2 ,p3]
            maxp = np.max(ps)
            ps.pop(np.argmax(ps))
            ps.insert(0, maxp)
            gt = np.sum(gts[i, :-1])
            ordered_ps = []
            for j in range(3):
                if j == 0:
                    gp_out.write('\t')
                if j == gt:
                    val = ps[0]
                else:
                    val = ps.pop()
                ordered_ps.append(val)
                gp_out.write('{:.2}'.format(val))
                if j != 2:
                    gp_out.write(',')
            total_dosage = ordered_ps[1] + 2*ordered_ps[2]
            if gts[i, 0] == 0:
                ap1 = random.uniform(max(0, total_dosage - 1), min(total_dosage, 0.5))
            else:
                ap1 = random.uniform(max(0.5, total_dosage - 1), min(total_dosage, 1))
            ap2 = total_dosage - ap1
            assert 0 <= ap2 <= 1
            ap1_out.write('\t{:.10}'.format(ap1))
            ap2_out.write('\t{:.10}'.format(ap2))
        gp_out.write('\n')
        ap1_out.write('\n')
        ap2_out.write('\n')

cmd = (
    'bash -c "'
    'bgzip -f gp_dosages.tsv && '
    'bgzip -f ap1_dosages.tsv && '
    'bgzip -f ap2_dosages.tsv && '
    'bcftools annotate '
        '-a gp_dosages.tsv.gz '
        '-h gp_dosage_header.hdr '
        '-S <(tail -n +2 samples.txt) '
        '-c CHROM,FROM,TO,FMT/GP '
        'many_samples_biallelic.vcf.gz | '
    'bcftools annotate '
        '-a ap1_dosages.tsv.gz '
        '-h ap1_dosage_header.hdr '
        '-S <(tail -n +2 samples.txt) '
        '-c CHROM,FROM,TO,FMT/AP1 '
        '- | '
    'bcftools annotate '
        '-a ap2_dosages.tsv.gz '
        '-h ap2_dosage_header.hdr '
        '-S <(tail -n +2 samples.txt) '
        '-c CHROM,FROM,TO,FMT/AP2 '
        '- > many_samples_biallelic_dosages.vcf && '
    'bgzip -f many_samples_biallelic_dosages.vcf && '
    'tabix -f many_samples_biallelic_dosages.vcf.gz '
    '"'
)
sp.run(cmd, shell = True, check=True, cwd=str(SCRIPT_DIR))

# multiallelic

vcf = cyvcf2.VCF(str(SCRIPT_DIR / 'many_samples_multiallelic.vcf.gz'))
samples = vcf.samples

with open(str(SCRIPT_DIR / 'ap1_multi_dosages.tsv'), 'w') as ap1_out, open(str(SCRIPT_DIR / 'ap2_multi_dosages.tsv'), 'w') as ap2_out:
    for var in vcf:
        ap1_out.write('{}\t{}\t{}'.format(var.CHROM, var.POS, var.POS))
        ap2_out.write('{}\t{}\t{}'.format(var.CHROM, var.POS, var.POS))
        gts = var.genotype.array()
        for i in range(gts.shape[0]):
            if -1 in gts[i, :]:
                ap1_out.write('\t.')
                ap2_out.write('\t.')
                continue
            for idx, out in enumerate([ap1_out, ap2_out]):
                max_ap = random.uniform(1/3, 1)
                small_ap_1 = random.uniform(max(0, 1 - 2*max_ap), min(1-max_ap, max_ap))
                small_ap_2 = 1 - max_ap - small_ap_1
                assert 0 <= small_ap_2
                if gts[i, idx] == 0:
                    aps = [small_ap_1, small_ap_2]
                else:
                    aps = [small_ap_1]
                    aps.insert(gts[i, idx]-1, max_ap)
                out.write('\t{:.2},{:.2}'.format(aps[0], aps[1]))
        ap1_out.write('\n')
        ap2_out.write('\n')

cmd = (
    'bash -c "'
    'bgzip -f ap1_multi_dosages.tsv && '
    'bgzip -f ap2_multi_dosages.tsv && '
    'bcftools annotate '
        '-a ap1_multi_dosages.tsv.gz '
        '-h ap1_dosage_header.hdr '
        '-S <(tail -n +2 samples.txt) '
        '-c CHROM,FROM,TO,FMT/AP1 '
        'many_samples_multiallelic.vcf.gz | '
    'bcftools annotate '
        '-a ap2_multi_dosages.tsv.gz '
        '-h ap2_dosage_header.hdr '
        '-S <(tail -n +2 samples.txt) '
        '-c CHROM,FROM,TO,FMT/AP2 '
        '- > many_samples_multiallelic_dosages.vcf && '
    'bgzip -f many_samples_multiallelic_dosages.vcf && '
    'tabix -f many_samples_multiallelic_dosages.vcf.gz '
    '"'
)
sp.run(cmd, shell = True, check=True, cwd=str(SCRIPT_DIR))

