import pandas as pd

import argparse
import os

def join_trimmomatic_and_key_metrics(options):
    samples_counter = 0
    sample_ids = []
    fnames = os.listdir(options.sample_trimmomatic_stats_dir)
    for fname in fnames:
        if not fname.endswith('.csv'):
            continue

        fpath = os.path.join(options.sample_trimmomatic_stats_dir, fname)
        sample_id = fname.replace('__trimlog_stats.csv', '')

        if samples_counter == 0:
            df_samples_trim = pd.read_csv(fpath, header=0)
        else:
            df_sample_trim = pd.read_csv(fpath, header=0)
            df_samples_trim = pd.concat([df_samples_trim, df_sample_trim], ignore_index=True)

        
        sample_ids.append(sample_id)
        samples_counter += 1
    df_samples_trim['Sample_ID'] = sample_ids

    df_key_metrics = pd.read_csv(options.sample_metrics_file, sep='\t', header=0)
    
    df_merged = pd.merge(df_key_metrics, df_samples_trim, how='left', on='Sample_ID')

    # op_fpath = options.sample_metrics_file.replace('.txt', '_with_trimmomatic.txt')
    df_merged.to_csv(options.outfile, index=False)

    print ("Output file written to: {}".format(options.outfile))        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the options.')
    parser.add_argument('--sample_metrics_file', type = str, required = True, help ='Sample-wise metrics file')
    parser.add_argument('--sample_trimmomatic_stats_dir', type = str, required = True, help ='Sample-wise trimmomatic stats file')
    parser.add_argument('--outfile', type = str, required = True, help ='Output file')
    options = parser.parse_args()

    join_trimmomatic_and_key_metrics(options)