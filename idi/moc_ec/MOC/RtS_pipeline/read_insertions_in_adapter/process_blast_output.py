import pandas as pd
import os
import numpy as np
import linecache
import argparse
from tqdm import tqdm
import pdb

import plotly.graph_objects as go

# csv columns - 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'

def process_rts(df, args):
     # params
    wb_len = 9 # well barcode length (including the terminal A/T)
    index1_len  = 33 # TGTGAGAAAGGGATGTGCTGCGAGAAGGCTAGA
    p5_len = 8

    read_start_pos_list = []
    for iter in tqdm(df.groupby('query_seqid')):
        query_seqid, row_hits = iter
        if 'index1' in row_hits['sub_seqid'].values:
            # Filter only rows with Index1 adapter match and sort those rows by evalue
            row_hit = row_hits[row_hits['sub_seqid'] == 'index1'].sort_values(by='evalue')

            # We are finding the starting position of the read (query) wrt to the adapter. If start position is 0, it means that the read starts from the well-barcode.
            # If start position is wb_len, then the read starts from Index1 adapter
            # If the start position is negative, then the read has DNA molecule sequence. 
            # Therefore, if the start position is negative, then the absolute value of the start position indicates the length of the DNA molecule in the read
            read_start_pos = row_hit['query_start'].values[0] - row_hit['sub_start'].values[0] - wb_len
            read_start_pos = -read_start_pos
        else:
            # Should only contain rows with matches to TruSeq adapter
            # Sort these rows by evalue
            row_hit = row_hits.sort_values(by='evalue')

            # We are finding the starting position of the read (query) wrt to the adapter. If start position is 0, it means that the read starts from the well-barcode.
            # If start position is wb_len, then the read starts from Index1 adapter
            # If start position is wb_len + index1_len, then the read starts from P5 adapter
            # If start position is wb_len + index1_len + P5 len, then the read starts from TruSeq adapter
            # If the start position is negative, then the read has DNA molecule sequence. 
            # Therefore, if the start position is negative, then the absolute value of the start position indicates the length of the DNA molecule in the read
            read_start_pos = row_hit['query_start'].values[0] - row_hit['sub_start'].values[0] - p5_len - index1_len - wb_len
            read_start_pos = -read_start_pos

        read_start_pos_list.append(read_start_pos)

    read_start_pos_list = np.array(read_start_pos_list)

    dna = read_start_pos_list[read_start_pos_list < 0]
    wb = read_start_pos_list[(read_start_pos_list >= 0) & (read_start_pos_list < 9)]
    index1 = read_start_pos_list[(read_start_pos_list >= 9) & (read_start_pos_list < 42)]
    p5 = read_start_pos_list[(read_start_pos_list >= 42) & (read_start_pos_list < 50)]
    truseq = read_start_pos_list[read_start_pos_list >= 50]

    # total_reads = len(list(SeqIO.parse(args.fasta, "fasta")))
    total_reads = sum(1 for _ in linecache.getlines(args.fasta)) // 2
    print ('Total number of reads:', total_reads)

    reads_with_adapter = len(read_start_pos_list)
    print ('Reads containing adapter:', reads_with_adapter)

    reads_with_only_adapter = len(read_start_pos_list) - len(dna)
    print ('Reads that contain only adapter (no DNA):', reads_with_only_adapter)

    exp_name = os.path.basename(args.fasta)

    fig =go.Figure(go.Sunburst(
        labels=['Total reads<br>{}'.format(total_reads), 'No Adapter seq<br>{}'.format(total_reads-reads_with_adapter), 'Has adapter seq<br>{}'.format(reads_with_adapter), 'Has DNA<br>{}'.format(reads_with_adapter-reads_with_only_adapter), 'No DNA<br>{}'.format(reads_with_only_adapter)],
        parents=['', 'Total reads<br>{}'.format(total_reads), 'Total reads<br>{}'.format(total_reads), 'Has adapter seq<br>{}'.format(reads_with_adapter), 'Has adapter seq<br>{}'.format(reads_with_adapter)],
        values=[total_reads, total_reads-reads_with_adapter, reads_with_adapter, reads_with_adapter-reads_with_only_adapter, reads_with_only_adapter],
        branchvalues="total",
    ))
    fig.update_layout(margin = dict(t=0, l=0, r=0, b=0), title_text='Pie-chart distribution of Reads ({})'.format(exp_name))
    op_fpath = args.csv.replace('.csv','_pie_chart.png')
    fig.write_image(op_fpath)
    print ('Pie chart generated at', op_fpath)


    fig = go.Figure()
    fig.add_trace(go.Histogram(x=truseq, name='TruSeq'))
    fig.add_trace(go.Histogram(x=p5, name='P5 (8)'))
    fig.add_trace(go.Histogram(x=index1, name='Index1 (33)'))
    fig.add_trace(go.Histogram(x=wb, name='well-barcode (9)', xbins=dict(size=1)))
    fig.add_trace(go.Histogram(x=dna, name='DNA'))

    # The two histograms are drawn on top of another
    fig.update_layout(
        barmode='stack', 
        title_text='Reads with adapter seqs ({})'.format(exp_name),
        xaxis_title_text='Read Starting Positions',
        yaxis_title_text='Count')
    op_fpath = args.csv.replace('.csv','_histogram.png')
    fig.write_image(op_fpath)
    print ('Histogram generated at', op_fpath)

def process_scr(df, args):
    # params
    umi_len = 8
    wb_len = 13
    truseq_len = 47
    bb_len = 16

    read_start_pos_list = []
    for iter in tqdm(df.groupby('query_seqid')):
        query_seqid, row_hits = iter
        if 'TruSeq' in row_hits['sub_seqid'].values:
            # Filter only rows with TruSeq adapter match and sort those rows by evalue
            row_hit = row_hits[row_hits['sub_seqid'] == 'TruSeq'].sort_values(by='evalue')

            # We are finding the starting position of the read (query) wrt to the adapter. If start position is 0, it means that the read starts from the well-barcode.
            # If start position is umi_len+wb_len, then the read starts from TruSeq adapter
            # If the start position is negative, then the read has DNA molecule sequence. 
            # Therefore, if the start position is negative, then the absolute value of the start position indicates the length of the DNA molecule in the read
            read_start_pos = row_hit['query_start'].values[0] - row_hit['sub_start'].values[0] - umi_len - wb_len
            read_start_pos = -read_start_pos
        else:
            # Should only contain rows with matches to Index2 adapter
            # Sort these rows by evalue
            row_hit = row_hits.sort_values(by='evalue')

            # We are finding the starting position of the read (query) wrt to the adapter. If start position is 0, it means that the read starts from the well-barcode.
            # If start position is umi_len+wb_len, then the read starts from TruSeq adapter
            # If start position is umi_len+wb_len+truseq_len+bb_len, then the read starts from Index2 or P5 adapter
            # If the start position is negative, then the read has DNA molecule sequence. 
            # Therefore, if the start position is negative, then the absolute value of the start position indicates the length of the DNA molecule in the read
            read_start_pos = row_hit['query_start'].values[0] - row_hit['sub_start'].values[0] - bb_len - truseq_len - umi_len - wb_len
            read_start_pos = -read_start_pos

        read_start_pos_list.append(read_start_pos)

    read_start_pos_list = np.array(read_start_pos_list)

    dna = read_start_pos_list[read_start_pos_list < 0]
    wb = read_start_pos_list[(read_start_pos_list >= 0) & (read_start_pos_list < 13)]
    umi = read_start_pos_list[(read_start_pos_list >= 13) & (read_start_pos_list < 21)]
    truseq = read_start_pos_list[(read_start_pos_list >= 21) & (read_start_pos_list < 54)]
    boa = read_start_pos_list[(read_start_pos_list >= 54) & (read_start_pos_list < 68)]
    bb = read_start_pos_list[(read_start_pos_list >= 68) & (read_start_pos_list < 84)]
    index2 = read_start_pos_list[read_start_pos_list >= 84]

    # total_reads = len(list(SeqIO.parse(args.fasta, "fasta")))
    total_reads = sum(1 for _ in linecache.getlines(args.fasta)) // 2
    print ('Total number of reads:', total_reads)

    reads_with_adapter = len(read_start_pos_list)
    print ('Reads containing adapter:', reads_with_adapter)

    reads_with_only_adapter = len(read_start_pos_list) - len(dna)
    print ('Reads that contain only adapter (no DNA):', reads_with_only_adapter)

    exp_name = os.path.basename(args.fasta)

    fig =go.Figure(go.Sunburst(
        labels=['Total reads<br>{}'.format(total_reads), 'No Adapter seq<br>{}'.format(total_reads-reads_with_adapter), 'Has adapter seq<br>{}'.format(reads_with_adapter), 'Has DNA<br>{}'.format(reads_with_adapter-reads_with_only_adapter), 'No DNA<br>{}'.format(reads_with_only_adapter)],
        parents=['', 'Total reads<br>{}'.format(total_reads), 'Total reads<br>{}'.format(total_reads), 'Has adapter seq<br>{}'.format(reads_with_adapter), 'Has adapter seq<br>{}'.format(reads_with_adapter)],
        values=[total_reads, total_reads-reads_with_adapter, reads_with_adapter, reads_with_adapter-reads_with_only_adapter, reads_with_only_adapter],
        branchvalues="total",
    ))
    fig.update_layout(margin = dict(t=0, l=0, r=0, b=0), title_text='Pie-chart distribution of Reads ({})'.format(exp_name))
    op_fpath = args.csv.replace('.csv','_pie_chart.png')
    fig.write_image(op_fpath)
    print ('Pie chart generated at', op_fpath)


    fig = go.Figure()
    fig.add_trace(go.Histogram(x=index2, name='Index2'))
    fig.add_trace(go.Histogram(x=bb, name='bead-barcode (16)'))
    fig.add_trace(go.Histogram(x=boa, name='Boa (14)'))
    fig.add_trace(go.Histogram(x=truseq, name='TruSeq (33)'))
    fig.add_trace(go.Histogram(x=umi, name='UMI (8)'))
    fig.add_trace(go.Histogram(x=wb, name='well-barcode (13)', xbins=dict(size=1)))
    fig.add_trace(go.Histogram(x=dna, name='DNA'))

    # The two histograms are drawn on top of another
    fig.update_layout(
        barmode='stack', 
        title_text='Reads with adapter seqs ({})'.format(exp_name),
        xaxis_title_text='Read Starting Positions',
        yaxis_title_text='Count')
    op_fpath = args.csv.replace('.csv','_histogram.png')
    fig.write_image(op_fpath)
    print ('Histogram generated at', op_fpath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-csv', type=str, help='Blast output file (csv format)')
    parser.add_argument('-fasta', type=str, help='Reads in fasta format')
    parser.add_argument('-lc_method', type=str, help='Protocol (SCR, RtS, etc)')
    args = parser.parse_args()

    # In this program, query are the reads and subject are the adapters
    df = pd.read_csv(args.csv, header=None)
    df.columns = ['query_seqid', 'sub_seqid', 'percent_identity', 'alignment_length', 'mismatch', 'gapopen', 'query_start', 'query_end', 'sub_start', 'sub_end', 'evalue', 'bitscore']

    # Making the positions 0-indexed
    df['query_start'] = df['query_start']-1
    df['query_end'] = df['query_end']-1
    df['sub_start'] = df['sub_start']-1
    df['sub_end'] = df['sub_end']-1

    if 'scr' in args.lc_method.lower():
        process_scr(df, args)
    elif 'rts' in args.lc_method.lower():
        process_rts(df, args)
    else:
        raise ValueError("LC method neither contains scr or rts")

    






