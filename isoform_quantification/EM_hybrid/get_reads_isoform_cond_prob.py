import pandas as pd
import pickle
from pathlib import Path
import multiprocessing as mp
import numpy as np
import glob
def get_cdf(read_len_dist,l):
    return read_len_dist.loc[read_len_dist.index<l,'PDF'].sum()/read_len_dist['PDF'].sum()
def get_Sm_dict(read_len_dist,isoform_df):
    descending_read_len_dist = read_len_dist.sort_index(ascending=False)
    Sm_dict = {}
    previous_Sm= 0
    previous_read_len = 0
    first_read_len = True
    for read_len in descending_read_len_dist.index:
        if first_read_len:
            sub_isoform_df = isoform_df[(isoform_df['isoform_len']>=read_len)].copy()
        else:
            sub_isoform_df = isoform_df[(isoform_df['isoform_len']>=read_len) & (isoform_df['isoform_len']<previous_read_len)].copy()
        sub_isoform_len_set = set(sub_isoform_df['isoform_len'].values)
        sub_isoform_len_cdf_dict = {}
        for isoform_len in sub_isoform_len_set:
            sub_isoform_len_cdf_dict[isoform_len] = get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1) 
        if sub_isoform_df.shape[0] != 0:
            sub_isoform_df['sum'] = sub_isoform_df.apply(lambda rows:rows['theta']/sub_isoform_len_cdf_dict[rows['isoform_len']],axis=1)
            Lm = sub_isoform_df['sum'].sum()
        else:
            Lm = 0
        if first_read_len:
            Sm_dict[read_len] = Lm
            first_read_len = False
        else:
            Sm_dict[read_len] = previous_Sm + Lm
        previous_Sm = Sm_dict[read_len]
        previous_read_len = read_len
    return Sm_dict
def get_all_reads_isoform_cond_prob_LIQA_modified(args):
    worker_id,output_path,isoform_df,read_len_dist,Sm_dict = args
    for fpath in glob.glob(f'{output_path}/temp/LR_alignments/{worker_id}_*'):
        batch_id = fpath.split('/')[-1].split('_')[-1]
        with open(fpath,'rb') as f:
            [reads_isoform_info,read_len_dict,_,_,_] = pickle.load(f)
        all_reads_isoform_cond_prob = {}
        for read in reads_isoform_info:
            all_reads_isoform_cond_prob[read] = {}
            read_len = read_len_dict[read]
            for isoform in reads_isoform_info[read]:
                isoform_len = isoform_df.loc[isoform,'isoform_len']
                pdf = read_len_dist.loc[read_len,'PDF']/read_len_dist['PDF'].sum()
                cdf = get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1)
                Sm = Sm_dict[read_len]
                if cdf == 0 or Sm == 0:
                    cond_prob = pdf
                else:
                    cond_prob = pdf/(get_cdf(read_len_dist,isoform_len+1) - get_cdf(read_len_dist,1)) / Sm_dict[read_len]
                if isoform_len - read_len + 1 != 0:
                    cond_prob /= isoform_len - read_len + 1
                else:
                    cond_prob = 0
                all_reads_isoform_cond_prob[read][isoform] = cond_prob
        rows = []
        MIN_PROB = 1e-100
        for read in all_reads_isoform_cond_prob:
            for isoform,prob in all_reads_isoform_cond_prob[read].items():
                rows.append([read,isoform,prob])
        cond_prob_df = pd.DataFrame(rows,columns=['read','isoform','cond_prob']).set_index(['read','isoform'])
        cond_prob_df.loc[(cond_prob_df['cond_prob'] == MIN_PROB) | (cond_prob_df['cond_prob'] == np.float('inf')),'cond_prob'] = 0
        Path(f'{output_path}/temp/cond_prob/').mkdir(exist_ok=True,parents=True)
        with open(f'{output_path}/temp/cond_prob/{worker_id}_{batch_id}','wb') as f:
            pickle.dump(cond_prob_df,f)
#     return all_reads_isoform_cond_prob
def get_cond_prob_MT_LIQA_modified(threads,output_path,isoform_df,read_len_dist,Sm_dict):
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = worker_id,output_path,isoform_df,read_len_dist,Sm_dict
        futures.append(pool.apply_async(get_all_reads_isoform_cond_prob_LIQA_modified,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
