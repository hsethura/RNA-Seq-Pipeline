import pandas as pd
import numpy as np
import argparse

def generate_pool_wise_metrics(options):
    print ('Generating Poolwise metrics by parsing: {}'.format(options.sample_metrics_file))
    # sample-wise data frame
    df_sample = pd.read_csv(options.sample_metrics_file, sep='\t')

    # dict with elements to add for each pool
    pool_wise_dict = {
        'MOCP_ID': [],
        'Project_ID': [],
        'Plate/Box_ID': [],
        'Pool_ID': [],
        'pool_total_reads_assigned_to_bcs': [],
        'pool_total_reads_including_unassigned': [],
        'pool_total_frags_counted': [],
        'pool_pcnt_mapped_to_bcs': [],
        'pool_pcnt_aligned': [],
        'pool_pcnt_properly_mapped_pairs': [],
        'pool_average_insert_len': [],
        'pool_pcnt_sense': [],
        'pool_CDS_total_counts_for_replicon': [],
        'pool_CDS_pcnt_of_counted': [],
        'pool_CDS_pcnt_sense': [],
        'pool_rRNA_pcnt_of_counted': [],
        'pool_rRNA_pcnt_sense': [],
        'pool_misc_RNA_pcnt_of_counted': [],
        'pool_misc_RNA_pcnt_sense': [],
        'pool_tRNA_pcnt_of_counted': [],
        'pool_tRNA_pcnt_sense': [],
        'pool_IGR_pcnt_of_counted': [],
        'pool_percent_trimmed_read1s': [],
        'pool_percent_trimmed_read2s': [],
        'pool_percent_empty_read1s': [],
        'pool_percent_empty_read2s': [],
        'pool_percent_trimmed_read_pairs': [],
        'pool_percent_empty_read_pairs': [],
    }

    # columns that are generic to all samples in the pool
    generic_cols = ['MOCP_ID', 'Project_ID', 'Plate/Box_ID', 'Pool_ID']

    # Iterate over each pool
    for pool_id, pool_rows in df_sample.groupby('Pool_ID'):
        for col in generic_cols:
            pool_wise_dict[col].append(pool_rows[col].values[0])
        
        # 'Pcnt_bc_in_pool' represents the percent of sequences in the pool that are mapped to a corresponding sample. 
        # This mapping is based on the inline/barcode sequence. 
        # Since there are sequences with barcodes that do not map to any recognized inline seq (ambiguous/no match), 'Pcnt_bc_in_pool' does not sum up to 1
        pcnt_bc_in_pool = pool_rows['Pcnt_bc_in_pool'].values

        # percent of reads whose sequences have a valid inline sequence 
        pool_pcnt_mapped_to_bcs = np.sum(pcnt_bc_in_pool)
        
        # Normalizing this column, so that it adds to 1
        if np.sum(pcnt_bc_in_pool) != 0:
            pcnt_bc_in_pool_normalized = pcnt_bc_in_pool*100 / np.sum(pcnt_bc_in_pool)
        else:
            # since there are no sequences assinged to the pool, every value is zero
            pcnt_bc_in_pool_normalized = pcnt_bc_in_pool
        
        # As explained above, since 'Pcnt_bc_in_pool' represents the percent of sequences in the pool that are mapped to a corresponding sample. 
        pcnt_sample_in_pool_normalized = pcnt_bc_in_pool_normalized
        
        pool_total_reads_assigned_to_bcs = np.sum(pool_rows['Total_reads'].values)
        pool_total_frags_counted = np.sum(pool_rows['Total_frags_counted'].values)

        pool_total_reads_including_unassigned = int (pool_total_reads_assigned_to_bcs * 100.0 / pool_pcnt_mapped_to_bcs)
        
        pool_pcnt_aligned = np.average(pool_rows['pcnt_aligned'].values, weights=pcnt_sample_in_pool_normalized)
        
        # The below two metrics are not calculated accurately. The weights given are %reads in the sample, we rather need %fragments in the sample
        pool_pcnt_properly_mapped_pairs = np.average(pool_rows['pcnt_properly_mapped_pairs'].values, weights=pcnt_sample_in_pool_normalized)
        pool_average_insert_len = np.average(pool_rows['average_insert_len'].values, weights=pcnt_sample_in_pool_normalized)
        
        # number of aligned reads per sample
        aligned_reads = (pool_rows['Total_reads'].values * pool_rows['pcnt_aligned'].values / 100).astype(np.int64)
        
        # number of aligned sense reads per sample
        aligned_sense_reads = (aligned_reads * pool_rows['pcnt_sense'].values / 100).astype(np.int64)
        
        # percent of aligned reads that are in sense direction
        if  np.sum(aligned_reads) != 0:
            pool_pcnt_sense = 100 * np.sum(aligned_sense_reads).astype(np.float) / np.sum(aligned_reads)
        else:
            pool_pcnt_sense = 0
        
        # CDS
        pool_CDS_total_counts_for_replicon = np.sum(pool_rows['CDS_total_counts_for_replicon'].values)
        
        if pool_total_frags_counted != 0:
            pool_CDS_pcnt_of_counted = 100 * float(pool_CDS_total_counts_for_replicon) / pool_total_frags_counted
        else:
            pool_CDS_pcnt_of_counted = 0
        
        num_CDS_sense = (pool_rows['CDS_pcnt_sense'].values * pool_rows['CDS_total_counts_for_replicon'].values / 100).astype(np.int64)
        if pool_CDS_total_counts_for_replicon != 0:
            pool_CDS_pcnt_sense = 100 * np.sum(num_CDS_sense).astype(np.float) / pool_CDS_total_counts_for_replicon
        else:
            pool_CDS_pcnt_sense = 0
        
        # rRNA
        num_rRNA_frags = (pool_rows['rRNA_pcnt_of_counted'].values * pool_rows['Total_frags_counted'].values / 100).astype(np.int64)
        if pool_total_frags_counted != 0:
            pool_rRNA_pcnt_of_counted = 100 * np.sum(num_rRNA_frags).astype(np.float) / pool_total_frags_counted
        else:
            pool_rRNA_pcnt_of_counted = 0
        
        num_rRNA_sense = (pool_rows['rRNA_pcnt_sense'].values * num_rRNA_frags / 100).astype(np.int64)
        if np.sum(num_rRNA_frags) != 0:
            pool_rRNA_pcnt_sense = 100 * np.sum(num_rRNA_sense).astype(np.float) / np.sum(num_rRNA_frags)
        else:
            pool_rRNA_pcnt_sense = 0
        
        # misc_RNA
        num_misc_RNA_frags = (pool_rows['misc_RNA_pcnt_of_counted'].values * pool_rows['Total_frags_counted'].values / 100).astype(np.int64)
        if pool_total_frags_counted != 0:
            pool_misc_RNA_pcnt_of_counted = 100 * np.sum(num_misc_RNA_frags).astype(np.float) / pool_total_frags_counted
        else:
            pool_misc_RNA_pcnt_of_counted = 0
        
        num_misc_RNA_sense = (pool_rows['misc_RNA_pcnt_sense'].values * num_misc_RNA_frags / 100).astype(np.int64)
        if np.sum(num_misc_RNA_frags) != 0:
            pool_misc_RNA_pcnt_sense = 100 * np.sum(num_misc_RNA_sense).astype(np.float) / np.sum(num_misc_RNA_frags)
        else:
            pool_misc_RNA_pcnt_sense = 0
        
        # tRNA
        num_tRNA_frags = (pool_rows['tRNA_pcnt_of_counted'].values * pool_rows['Total_frags_counted'].values / 100).astype(np.int64)
        if pool_total_frags_counted != 0:
            pool_tRNA_pcnt_of_counted = 100 * np.sum(num_tRNA_frags).astype(np.float) / pool_total_frags_counted
        else:
            pool_tRNA_pcnt_of_counted = 0
        
        num_tRNA_sense = (pool_rows['tRNA_pcnt_sense'].values * num_tRNA_frags / 100).astype(np.int64)
        if np.sum(num_tRNA_frags) != 0:
            pool_tRNA_pcnt_sense = 100 * np.sum(num_tRNA_sense).astype(np.float) / np.sum(num_tRNA_frags)
        else:
            pool_tRNA_pcnt_sense = 0
        
        # IGR
        num_IGR_frags = (pool_rows['IGR_pcnt_of_counted'].values * pool_rows['Total_frags_counted'].values / 100).astype(np.int64)
        if pool_total_frags_counted != 0:
            pool_IGR_pcnt_of_counted = 100 * np.sum(num_IGR_frags).astype(np.float) / pool_total_frags_counted
        else:
            pool_IGR_pcnt_of_counted = 0

        try:
            pool_read1s = np.sum(pool_rows['total_read1s'].values).astype(np.float)
            pool_read2s = np.sum(pool_rows['total_read2s'].values).astype(np.float)
            pool_read_pairs = np.sum(pool_rows['total_read_pairs'].values).astype(np.float)

            pool_trimmed_read1s = np.sum(pool_rows['trimmed_read1s'].values).astype(np.float)
            pool_trimmed_read2s = np.sum(pool_rows['trimmed_read2s'].values).astype(np.float)
            pool_trimmed_read_pairs = np.sum(pool_rows['trimmed_read_pairs'].values).astype(np.float)

            pool_empty_read1s = np.sum(pool_rows['empty_read1s'].values).astype(np.float)
            pool_empty_read2s = np.sum(pool_rows['empty_read2s'].values).astype(np.float)
            pool_empty_read_pairs = np.sum(pool_rows['empty_read_pairs'].values).astype(np.float)
        except KeyError:
            # if it enters this block, it means trimmomatic has not been run at the sample level
            pool_read1s = 0
            pool_read2s = 0
            pool_read_pairs = 0
            

        if pool_read1s != 0:
            pool_percent_trimmed_read1s = 100 * pool_trimmed_read1s / pool_read1s
            pool_percent_empty_read1s = 100 * pool_empty_read1s / pool_read1s
        else:
            pool_percent_trimmed_read1s = 0
            pool_percent_empty_read1s = 0

        if pool_read2s != 0:
            pool_percent_trimmed_read2s = 100 * pool_trimmed_read2s / pool_read2s
            pool_percent_empty_read2s = 100 * pool_empty_read2s / pool_read2s
        else:
            pool_percent_trimmed_read2s = 0
            pool_percent_empty_read2s = 0
        
        if pool_read_pairs != 0:
            pool_percent_trimmed_read_pairs = 100 * pool_trimmed_read_pairs / pool_read_pairs
            pool_percent_empty_read_pairs = 100 * pool_empty_read_pairs / pool_read_pairs
        else:
            pool_percent_trimmed_read_pairs = 0
            pool_percent_empty_read_pairs = 0

        # add pool information to dict
        pool_wise_dict['pool_total_reads_assigned_to_bcs'].append(pool_total_reads_assigned_to_bcs)
        pool_wise_dict['pool_total_reads_including_unassigned'].append(pool_total_reads_including_unassigned)
        pool_wise_dict['pool_total_frags_counted'].append(pool_total_frags_counted)
        pool_wise_dict['pool_pcnt_mapped_to_bcs'].append(pool_pcnt_mapped_to_bcs)
        pool_wise_dict['pool_pcnt_aligned'].append(pool_pcnt_aligned)
        pool_wise_dict['pool_pcnt_properly_mapped_pairs'].append(pool_pcnt_properly_mapped_pairs)
        pool_wise_dict['pool_average_insert_len'].append(pool_average_insert_len)
        pool_wise_dict['pool_pcnt_sense'].append(pool_pcnt_sense)
        pool_wise_dict['pool_CDS_total_counts_for_replicon'].append(pool_CDS_total_counts_for_replicon)
        pool_wise_dict['pool_CDS_pcnt_of_counted'].append(pool_CDS_pcnt_of_counted)
        pool_wise_dict['pool_CDS_pcnt_sense'].append(pool_CDS_pcnt_sense)
        pool_wise_dict['pool_rRNA_pcnt_of_counted'].append(pool_rRNA_pcnt_of_counted)
        pool_wise_dict['pool_rRNA_pcnt_sense'].append(pool_rRNA_pcnt_sense)
        pool_wise_dict['pool_misc_RNA_pcnt_of_counted'].append(pool_misc_RNA_pcnt_of_counted)
        pool_wise_dict['pool_misc_RNA_pcnt_sense'].append(pool_misc_RNA_pcnt_sense)
        pool_wise_dict['pool_tRNA_pcnt_of_counted'].append(pool_tRNA_pcnt_of_counted)
        pool_wise_dict['pool_tRNA_pcnt_sense'].append(pool_tRNA_pcnt_sense)
        pool_wise_dict['pool_IGR_pcnt_of_counted'].append(pool_IGR_pcnt_of_counted)
        pool_wise_dict['pool_percent_trimmed_read1s'].append(pool_percent_trimmed_read1s)
        pool_wise_dict['pool_percent_empty_read1s'].append(pool_percent_empty_read1s)
        pool_wise_dict['pool_percent_trimmed_read2s'].append(pool_percent_trimmed_read2s)
        pool_wise_dict['pool_percent_empty_read2s'].append(pool_percent_empty_read2s)
        pool_wise_dict['pool_percent_trimmed_read_pairs'].append(pool_percent_trimmed_read_pairs)
        pool_wise_dict['pool_percent_empty_read_pairs'].append(pool_percent_empty_read_pairs)

    # pool-wise data frame
    df_pool = pd.DataFrame(pool_wise_dict)

    cols = ['MOCP_ID','Project_ID','Plate/Box_ID','Pool_ID','pool_total_reads_assigned_to_bcs', 'pool_total_reads_including_unassigned', 'pool_total_frags_counted', 'pool_pcnt_mapped_to_bcs','pool_pcnt_aligned','pool_pcnt_properly_mapped_pairs','pool_average_insert_len','pool_pcnt_sense','pool_CDS_total_counts_for_replicon','pool_CDS_pcnt_of_counted','pool_CDS_pcnt_sense','pool_rRNA_pcnt_of_counted','pool_rRNA_pcnt_sense','pool_misc_RNA_pcnt_of_counted','pool_misc_RNA_pcnt_sense','pool_tRNA_pcnt_of_counted','pool_tRNA_pcnt_sense','pool_IGR_pcnt_of_counted','pool_percent_trimmed_read1s','pool_percent_empty_read1s','pool_percent_trimmed_read2s','pool_percent_empty_read2s','pool_percent_trimmed_read_pairs','pool_percent_empty_read_pairs']

    # Rearrange the columns
    df_pool = df_pool[cols]

    df_pool.to_csv(options.outfile, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the options.')
    parser.add_argument('--sample_metrics_file', type = str, required = True, help ='Sample-wise metrics file')
    parser.add_argument('--outfile', type = str, required = True, help ='Output file')
    options = parser.parse_args()

    generate_pool_wise_metrics(options)
    
