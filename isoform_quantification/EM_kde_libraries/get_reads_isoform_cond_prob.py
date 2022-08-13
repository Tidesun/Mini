from unittest.mock import patch
import pandas as pd
import dill as pickle
from pathlib import Path
import multiprocessing as mp
import numpy as np
import bisect
import glob
def get_cdf(read_len_dist,l):
    return read_len_dist.loc[read_len_dist.index<l,'PDF'].sum()/read_len_dist['PDF'].sum()
def get_isoform_len_samples_dict(args):
    list_of_isoform_len,worker_id,kde_path,output_path = args
    with open(kde_path,'rb') as f:
        kde_model = pickle.load(f)
    isoform_len_samples_dict = {}
    name = '{}-{}_samples'.format(str(list_of_isoform_len[0]),str(list_of_isoform_len[-1]))
    for isoform_len in list_of_isoform_len:
        eps = isoform_len * 0.01
        min_isoform_num_reads = 10**5
        isoform_samples = np.array([[],[],[]]).T
        num_iters = 0
        while isoform_samples.shape[0] < min_isoform_num_reads:
            samples = kde_model.sample(min_isoform_num_reads)
            isoform_samples = np.concatenate([isoform_samples,samples[np.abs(samples[:,0] - isoform_len) < eps]],axis=0)
            num_iters += 1
            if num_iters > 10**2:
                break
        isoform_len_samples_dict[isoform_len] = isoform_samples
    Path(f'{output_path}/temp/random_samples/').mkdir(exist_ok=True,parents=True)
    with open(f'{output_path}/temp/random_samples/{name}.pkl','wb') as f:
        pickle.dump(isoform_len_samples_dict,f)
def prepare_kde_random_samples(isoform_len_set,threads):
    # min_isoform_len = 150
    isoform_len_list = sorted(list(isoform_len_set))
    list_of_isoform_len = []
    for isoform_len in isoform_len_list:
        # if isoform_len >= min_isoform_len:
        list_of_isoform_len.append(isoform_len)
    chunksize, extra = divmod(len(list_of_isoform_len), threads)
    if extra:
        chunksize += 1
    all_list_of_isoform_len = []
    for i in range(threads - 1):
        all_list_of_isoform_len.append(list_of_isoform_len[i*chunksize:(i+1)*chunksize])
    all_list_of_isoform_len.append(list_of_isoform_len[(threads-1)*chunksize:])
    return all_list_of_isoform_len
def get_samples_file_path(isoform_len,st_samples_index_dict,ed_samples_index_dict,output_path):
    start_index = bisect.bisect_right(sorted(list(st_samples_index_dict.keys())), isoform_len)
    end_index = bisect.bisect_left(sorted(list(ed_samples_index_dict.keys())), isoform_len)
    if start_index == end_index:
        return np.array([])
    else:
        index = sorted(list(st_samples_index_dict.keys()))[end_index:start_index][0]
        with open('{}/temp/random_samples/{}_samples.pkl'.format(output_path,st_samples_index_dict[index]),'rb') as f:
            samples = pickle.load(f)
        return samples[isoform_len]
def get_all_reads_isoform_cond_prob_kde(args):
    worker_id,st_samples_index_dict,ed_samples_index_dict,read_len_dist,output_path = args
    with open(f'{output_path}/temp/{worker_id}','rb') as f:
        reads_isoform_info,_ = pickle.load(f)
    all_reads_isoform_cond_prob = {}
    eps = 1
    for read in reads_isoform_info:
        all_reads_isoform_cond_prob[read] = {}
        for isoform,read_info in reads_isoform_info[read].items():
            [isoform_len,five_end,three_end] = read_info
            isoform_samples = get_samples_file_path(isoform_len,st_samples_index_dict,ed_samples_index_dict,output_path)
            if isoform_samples.shape[0] == 0:
                read_len = isoform_len - five_end - three_end
                cond_prob = read_len_dist.loc[read_len,'PDF'] / read_len_dist.loc[read_len_dist.index<isoform_len,'PDF'].sum()
            else:
                num_read_samples = isoform_samples[(np.abs(isoform_samples[:,1] - five_end) < eps) & (np.abs(isoform_samples[:,2] - three_end) < eps),:].shape[0]
                cond_prob = num_read_samples/isoform_samples.shape[0] / (2*eps)**2
            all_reads_isoform_cond_prob[read][isoform] = cond_prob
    Path(f'{output_path}/temp/cond_prob/').mkdir(exist_ok=True,parents=True)
    with open(f'{output_path}/temp/cond_prob/{worker_id}','wb') as f:
        pickle.dump(all_reads_isoform_cond_prob,f)
def get_kde_random_samples(threads,output_path,kde_path,isoform_len_set):
    all_list_of_isoform_len = prepare_kde_random_samples(isoform_len_set,threads)
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = all_list_of_isoform_len[worker_id],worker_id,kde_path,output_path
        futures.append(pool.apply_async(get_isoform_len_samples_dict,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
    st_samples_index_dict = {}
    ed_samples_index_dict = {}
    for fpath in glob.glob(f'{output_path}/temp/random_samples/*_samples.pkl'):
        name = fpath.split('/')[-1].split('_samples.pkl')[0]
        [st_isoform_len,ed_isoform_len] = name.split('-')
        st_isoform_len = int(st_isoform_len)
        ed_isoform_len = int(ed_isoform_len)
        st_samples_index_dict[st_isoform_len] = name
        ed_samples_index_dict[ed_isoform_len] = name
    return st_samples_index_dict,ed_samples_index_dict
def get_cond_prob_kde_random_samples(threads,output_path,read_len_dist,kde_path,isoform_len_set):
    st_samples_index_dict,ed_samples_index_dict = get_kde_random_samples(threads,output_path,kde_path,isoform_len_set)
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = worker_id,st_samples_index_dict,ed_samples_index_dict,read_len_dist,output_path
        futures.append(pool.apply_async(get_all_reads_isoform_cond_prob_kde,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()

