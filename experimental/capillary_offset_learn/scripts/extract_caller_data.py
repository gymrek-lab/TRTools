import pandas as pd
import shlex, json
import subprocess
methods = ['advntr'] # eh
suffix = {'gangstr': '_merged.vcf.gz',
        'hipstr': '.vcf.gz',
        'advntr': '_merged.vcf.gz',
        }
dir_root = "/projects/ps-gymreklab" #"/gymreklab-tscc" # 
outs_root = '/%s/mousavi/results/1000genomes/'%dir_root


# Load product sizes to get list of samples and primer/locus ids
col_list = ["PrimerID", "SampleID"]
ps_df = pd.read_csv("../data/csv/1000g_product_sizes.csv", usecols=col_list)
samples_list = list(set(ps_df['SampleID']))
samples_list.remove('Reference')
samples_list.remove('reference')
# print(samples_list)
locus_id_list = list(set(ps_df['PrimerID']))
# print(locus_id_list)

# Load list of loci to get chrom, start, end
col_list = ["LocusID", "Chrom", "Start", "End", "Motif"]
lc_df = pd.read_csv("../data/csv/1000g_loci.csv", usecols=col_list)
def get_chrom_start_end_motif(lid):
    q_df = lc_df.query('LocusID == @lid')
    if len(q_df) == 0:
        return None, None, None, None
    return list(q_df['Chrom'])[0], list(q_df['Start'])[0], list(q_df['End'])[0], list(q_df['Motif'])[0]

# Load pop and spop for each sample
sm_df = pd.read_csv('/%s/mousavi/analysis/1000genomes/samples/all_samples.txt'%dir_root, sep='\t', names = ['SampleID', 'SubPop', 'SuperPop', 'Sex', 'CramPath'])
def get_spop_pop(sid):
    q_df = sm_df.query('SampleID == @sid')
    if len(q_df) == 0:
        return None, None
    return list(q_df['SuperPop'])[0], list(q_df['SubPop'])[0]

def convert_data_to_dict(proc_line, meth, motif):
    if meth == 'gangstr':
        ref = proc_line[0]
        ref_len = len(ref)

        gb = []
        for c in proc_line[1].split(','):
            gb.append(str(int(c) * len(motif) - ref_len))
        return {'cn': proc_line[1], 
            'ci': proc_line[2], 
            'gb': ','.join(gb),
            'q': proc_line[4]}

    elif meth == 'hipstr':
        if proc_line[2] == '.':
            return {}
        ref = proc_line[0]
        alts = proc_line[1].split(',')
        if '|' in proc_line[2]:
            gt = proc_line[2].split('|')
            gb = proc_line[3].split('|')
        elif '/' in proc_line[2]:
            gt = proc_line[2].split('/')
            gb = proc_line[3].split('/')
        else:
            raise ValueError('Invalid GT: ', proc_line[2])
        q = proc_line[4]
        alleles = []
        allele_lens = []
        cn = []
        for a in gt:
            if int(a) == 0:
                alleles.append(ref)
                allele_lens.append(str(len(ref)))
                cn.append(str(len(ref) / len(motif)))
            else:
                alleles.append(alts[int(a) - 1])
                allele_lens.append(str(len(alts[int(a) - 1])))
                cn.append(str(len(alts[int(a) - 1]) / len(motif)))
        return {'alleles': ','.join(alleles), 
            'allele_lens': ','.join(allele_lens),
            'gb': ','.join(gb),
            'cn': ','.join(cn),
            'q': q}
    elif meth == 'eh':
        ref_len = int(proc_line[0]) * len(motif)
        gb = []
        for c in proc_line[1].split('/'):
            gb.append(str(int(c) * len(motif) - ref_len))
        return {'cn': proc_line[1], 
            'ci': proc_line[2], 
            'gb': ','.join(gb)
            }
    elif meth == 'advntr':
        if proc_line[2] == '.':
            return {}
        ref = proc_line[0]
        alts = proc_line[1].split(',')
        if '|' in proc_line[2]:
            gt = proc_line[2].split('|')
        elif '/' in proc_line[2]:
            gt = proc_line[2].split('/')
        else:
            raise ValueError('Invalid GT: ', proc_line[2])

        alleles = []
        allele_lens = []
        cn = []
        for a in gt:
            if int(a) == 0:
                alleles.append(ref)
                allele_lens.append(str(len(ref)))
                cn.append(str(len(ref) / len(motif)))
            else:
                alleles.append(alts[int(a) - 1])
                allele_lens.append(str(len(alts[int(a) - 1])))
                cn.append(str(len(alts[int(a) - 1]) / len(motif)))
        gb = []
        for c in cn:
            gb.append(str(float(c) * len(motif) - len(ref)))

        return {'alleles': ','.join(alleles), 
            'allele_lens': ','.join(allele_lens),
            'gb': ','.join(gb),
            'cn': ','.join(cn)}
    return {}




# create a dictionary for output
out_data = {}
for meth in methods:
    out_data[meth] = {}



for meth in methods:
    meth_outs = outs_root + meth + '_outs/'
    for lid in locus_id_list:
        out_data[meth][lid] = {}
        chrom, start, end, motif = get_chrom_start_end_motif(lid)
        if chrom is None: 
            continue
        # if lid != 'CACNA1A':
        #     continue
        print(lid + '\t' + chrom + ':' + str(start) + '-' + str(end))

        # iterate over all samples
        for sample in samples_list:
            spop, pop = get_spop_pop(sample)
            if spop is None:
                continue
            if meth != 'eh':
                vcf_dir = meth_outs + '/'.join([spop, pop, 'merged', pop + suffix[meth]])
            else:
                vcf_dir = meth_outs + '/'.join(['eh-polymorphic--1kg-all-hg38','merged', pop + '_merged.vcf.gz'])

            proc_line = -1
            if meth == 'gangstr':
                p = subprocess.Popen(shlex.split("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%REF@%REPCN@%REPCI@%RC@%Q]\n'"\
                    .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample)),
                    stdout=subprocess.PIPE)
                p.wait()
                out_lines = p.stdout.readline()
                if out_lines is None:
                    continue
                if len(out_lines) == 0:
                    continue
                proc_line = out_lines.strip().decode("utf-8").split('@')
            elif meth == 'hipstr':
                p = subprocess.Popen(shlex.split("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%REF@%ALT@[%GT@%GB@%Q]]\n'"\
                    .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample)),
                    stdout=subprocess.PIPE)
                # print("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%REF|%ALT|[%Q]]\n'"\
                #     .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample))
                p.wait()
                out_lines = p.stdout.readline()
                if out_lines is None:
                    continue
                if len(out_lines) == 0:
                    continue
                proc_line = out_lines.strip().decode("utf-8").split('@')
            elif meth == 'eh':
                p = subprocess.Popen(shlex.split("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%INFO/REF@%REPCN@%REPCI]\n'"\
                    .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample)),
                    stdout=subprocess.PIPE)
                p.wait()
                out_lines = p.stdout.readline()
                if out_lines is None:
                    continue
                if len(out_lines) == 0:
                    continue
                proc_line = out_lines.strip().decode("utf-8").split('@')
            elif meth == 'advntr':
                p = subprocess.Popen(shlex.split("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%REF@%ALT@[%GT]]\n'"\
                    .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample)),
                    stdout=subprocess.PIPE)
                # print("bcftools query {vcf} -r {ch}:{st} -s {samp} -f '[%REF|%ALT|[%Q]]\n'"\
                #     .format(vcf=vcf_dir, ch=chrom, st=start, samp=sample))
                p.wait()
                out_lines = p.stdout.readline()
                if out_lines is None:
                    continue
                if len(out_lines) == 0:
                    continue
                proc_line = out_lines.strip().decode("utf-8").split('@')
            if sample not in out_data[meth][lid]:
                if proc_line != -1:
                    out_data[meth][lid][sample] = convert_data_to_dict(proc_line, meth, motif)
            else:
                raise ValueError('Duplicate sample:', sample)
    with open('../data/json/' + meth  + '_calls.json', 'w') as outjson:
        json.dump(out_data[meth], outjson, indent=4)
print(out_data)

