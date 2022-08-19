from pathlib import Path
import pandas as pd
import dill as pickle
import numpy as np
import multiprocessing as mp
import glob
from EM_SR.prepare_hits import prepare_hits
def E_step_MT(args):
#     MIN_PROB = 1e-100
    q = {}
    isoform_q = {}
    worker_id,isoform_df,output_path = args
    for fpath in glob.glob(f'{output_path}/temp/hits_dict/{worker_id}_*'):
        with open(fpath,'rb') as f:
            hits_dict = pickle.load(f)
        for read in hits_dict:
            sum_q = 0
            q[read] = {}
            for hit in hits_dict[read]:
                isoform = hit['isoform']
                q[read][isoform] = hit['ant'] * isoform_df.loc[isoform,'len_theta_product']
                sum_q += q[read][isoform]
            if sum_q != 0:
                for isoform,prob in q[read].items():
                    q[read][isoform] /=sum_q
        for read in q:
            for isoform,prob in q[read].items():
                if isoform not in isoform_q:
                    isoform_q[isoform] = 0
                isoform_q[isoform] += prob
    isoform_q_df = pd.Series(isoform_q)
    return isoform_q_df
def E_step(isoform_df,output_path,threads):
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        args = worker_id,isoform_df,output_path
        futures.append(pool.apply_async(E_step_MT,(args,)))
    list_of_isoform_q_df = []
    for future in futures:
        new_isoform_q_df = future.get()
        list_of_isoform_q_df.append(new_isoform_q_df)
    pool.close()
    pool.join()
    return list_of_isoform_q_df
def M_step(list_of_isoform_q_df,isoform_df):
    isoform_q_df = pd.Series(0,index=isoform_df.index)
    for single_thread_isoform_q_df in list_of_isoform_q_df:
        isoform_q_df = isoform_q_df.add(single_thread_isoform_q_df,fill_value=0)
    new_theta_df = isoform_q_df*isoform_df['eff_len'] /((isoform_df['eff_len']*isoform_df['theta']).sum() * isoform_q_df.sum())
    new_theta_df = new_theta_df/new_theta_df.sum()
    return new_theta_df
def get_eff_len_dict(eff_len_file):
    df = pd.read_csv(eff_len_file,sep='\t')
    eff_len_df = df.set_index('target_id')['eff_length']
    eff_len_dict =eff_len_df.to_dict()
    return eff_len_dict
def EM_algo_SR(SR_sam,eff_len_file,output_path,threads):
    eff_len_dict = get_eff_len_dict(eff_len_file)
    theta_df = prepare_hits(SR_sam,eff_len_dict,output_path,threads)
    num_iters = 10000
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    for i in range(num_iters):
        isoform_df = pd.DataFrame({'eff_len':eff_len_dict,'theta':theta_df})
        isoform_df = isoform_df.fillna(0)
        isoform_df['len_theta_product'] = isoform_df['eff_len'] * isoform_df['theta'] / ((isoform_df['eff_len'] * isoform_df['theta']).sum())
        min_diff = 1e-3
        list_of_isoform_q_df = E_step(isoform_df,output_path,threads)
        new_theta_df = M_step(list_of_isoform_q_df,isoform_df)
        list_of_isoform_q_df = []
        diff = np.abs(theta_df - new_theta_df)/new_theta_df
        diff[new_theta_df==0] = 0
        if i % 10 == 0:
            print(f'Iteration:{i}')
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max())
            theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')

        if diff[diff > min_diff].shape[0] == 0:
            print(f'Iteration:{i}')
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max())
            theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')
            break
        theta_df = new_theta_df.copy()
    TPM_df = (theta_df/theta_df.sum())*1e6
    TPM_df.name = 'TPM'
    TPM_df.index.name = 'Isoform'
    TPM_df.to_csv(f'{output_path}/EM_expression.out',sep='\t')