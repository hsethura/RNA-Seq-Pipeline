import os
import re
import collections
import pandas as pd


def get_trans_to_gene(trans_to_gene_file):
    trans_to_gene = dict() 
    reg_str = "\t"
    with open(trans_to_gene_file) as tgf:
        for line0 in tgf:
            line = line0.rstrip()
            parts = re.split(reg_str, line)        
            trans = parts[0]
            gene = parts[1]
            trans_to_gene[trans] = gene
            # print("transcript: " + trans + " gene: " + gene)    
    return trans_to_gene


def get_gene_to_frag(rpkm_file, trans_to_gene):
    print("rpkm_file: " + rpkm_file)
    skip_count = 5
    reg_str = "\s+"
    gene_to_frag = dict()
    with open(rpkm_file) as ff:
        for _ in xrange(skip_count):
            next(ff)
        for line0 in ff:
            line = line0.rstrip()
            parts = re.split(reg_str, line)
            trans = parts[0]
            frag_count_str = parts[6]
            frag_count = int(frag_count_str)
            gene = trans_to_gene[trans]
            if gene not in gene_to_frag:
                gene_to_frag[gene] = frag_count
            elif gene in gene_to_frag:
                gene_to_frag[gene] += frag_count
    #print(gene_to_frag)
    gene_to_frag_s = collections.OrderedDict(sorted(gene_to_frag.items()))
    gene_to_frag_df = pd.DataFrame(gene_to_frag_s.items())
    filename = os.path.basename(os.path.normpath(rpkm_file))
    file_pre = filename.replace("_rpkm.txt", "")
    gene_to_frag_df.columns = ['gene', file_pre]
    gene_series = gene_to_frag_df.iloc[:,0]
    fragc_series = gene_to_frag_df.iloc[:,1]
    return gene_series, fragc_series


def all_gene_to_frag(rpkm_dir, outfile, trans_to_gene):
    gene_series_ori = pd.Series()
    count_series_ori = pd.Series()
    df = pd.DataFrame()
    for lfile in os.listdir(rpkm_dir):
        if lfile.endswith("_rpkm.txt"):
            lfile_abs = rpkm_dir + "/" + lfile
            gene_series, count_series = get_gene_to_frag(lfile_abs, trans_to_gene)
            if gene_series_ori.empty:
                gene_series_ori = gene_series
                df[gene_series_ori.name] = gene_series_ori
                count_series_ori = count_series
            if not gene_series_ori.equals(gene_series):
                print(count_series_ori)
                print(count_series)
                raise ValueError('These two geneid_series are not equal: ' 
                  + count_series_ori.name + ' and ' + count_series.name)
            df[count_series.name] = count_series
    return df

            
trans_to_gene_file = 'Mouse_mm10_gene_transcript.txt' 
rpkm_dir = '/broad/hptmp/RNASeq_proj/nirmalya/Hannan1/tests/test4'
outfile = rpkm_dir + "/all_rpkm.tsv"
trans_to_gene = get_trans_to_gene(trans_to_gene_file)
df = all_gene_to_frag(rpkm_dir, outfile, trans_to_gene)
df.to_csv(outfile, sep = '\t', index = False)

