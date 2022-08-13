from unittest.mock import patch
import pandas as pd
import dill as pickle
from pathlib import Path
import multiprocessing as mp
import numpy as np
def get_cdf(read_len_dist,l):
    return read_len_dist.loc[read_len_dist.index<l,'PDF'].sum()/read_len_dist['PDF'].sum()
def get_all_reads_isoform_cond_prob_kde_score(args):
    worker_id,isoform_df,kde_path,output_path = args
    with open(kde_path,'rb') as f:
        kde_model = pickle.load(f)
    with open(f'{output_path}/temp/{worker_id}','rb') as f:
        reads_isoform_info,_ = pickle.load(f)
    all_reads_isoform_cond_prob = {}
    all_samples = []
    all_read_isoform = []
    for read in reads_isoform_info:
        for isoform,read_info in reads_isoform_info[read].items():
            all_samples.append(read_info)
            all_read_isoform.append((read,isoform))
    all_log_read_prob = kde_model.score_samples(all_samples)
    for (read,isoform),log_read_prob in zip(all_read_isoform,all_log_read_prob):
        isoform_len = isoform_df.loc[isoform,'isoform_len']
        log_isoform_len_prob = np.log(isoform_df.loc[isoform_df['isoform_len'] == isoform_len,'theta'].sum())
        cond_prob = log_read_prob - log_isoform_len_prob
        if cond_prob <= -300:
            cond_prob = 0 
        else:
            cond_prob = np.e**(cond_prob)
        if read not in all_reads_isoform_cond_prob:
            all_reads_isoform_cond_prob[read] = {}
        all_reads_isoform_cond_prob[read][isoform] = cond_prob
    Path(f'{output_path}/temp/cond_prob/').mkdir(exist_ok=True,parents=True)
    with open(f'{output_path}/temp/cond_prob/{worker_id}','wb') as f:
        pickle.dump(all_reads_isoform_cond_prob,f)
#     return all_reads_isoform_cond_prob
def get_cond_prob_kde_score(threads,output_path,isoform_df,kde_path):
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = worker_id,isoform_df,kde_path,output_path 
        futures.append(pool.apply_async(get_all_reads_isoform_cond_prob_kde_score,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()

