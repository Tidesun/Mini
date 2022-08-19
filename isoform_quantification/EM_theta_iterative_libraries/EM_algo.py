import multiprocessing as mp
from EM_theta_iterative_libraries.get_reads_isoform_cond_prob import get_all_reads_isoform_cond_prob_LIQA_modified,get_all_reads_isoform_cond_prob_LIQA,get_Sm_dict
from EM_libraries.get_reads_isoform_info import get_reads_isoform_info,get_read_len_dist
from EM_libraries.prepare_MT import prepare_MT
import dill as pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from pathlib import Path
import gc

def E_step_MT(args):
    MIN_PROB = 1e-100
    EM_choice = args[-2]
    if EM_choice == 'LIQA_modified':
        worker_id,theta_df,isoform_df,Sm_dict,read_len_dist,EM_choice,output_path = args
        all_reads_isoform_cond_prob = get_all_reads_isoform_cond_prob_LIQA_modified((worker_id,output_path,isoform_df,read_len_dist,Sm_dict))
    elif EM_choice == "LIQA":
        worker_id,theta_df,isoform_df,read_len_dist,EM_choice,output_path = args
        all_reads_isoform_cond_prob = get_all_reads_isoform_cond_prob_LIQA((worker_id,output_path,isoform_df,read_len_dist))
    # with open(f'{output_path}/temp/cond_prob/{worker_id}','rb') as f:
    #     all_reads_isoform_cond_prob = pickle.load(f)
    q = {}
    for read in all_reads_isoform_cond_prob:
        q[read] = {}
        sum_q = 0
        for isoform,prob in all_reads_isoform_cond_prob[read].items():
            if prob <= MIN_PROB or prob == float('inf'):
                prob = 0
            isoform_expression = theta_df[isoform]
            q[read][isoform] = prob*isoform_expression
            sum_q += q[read][isoform]
        try:
            if sum_q != 0:
                for isoform,prob in q[read].items():
                    q[read][isoform] /=sum_q
        except Exception as e:
            print(sum_q)
            print(all_reads_isoform_cond_prob[read])
            raise e
    isoform_q = {}
    for read in q:
        for isoform,prob in q[read].items():
            if isoform not in isoform_q:
                isoform_q[isoform] = 0
            isoform_q[isoform] += prob
    isoform_q_df = pd.Series(isoform_q)
    return isoform_q_df
def E_step(theta_df,threads,output_path,EM_choice,isoform_df,read_len_dist,Sm_dict=None):
    pool = mp.Pool(threads)
    futures = []
    for worker_id in range(threads):
        if EM_choice == 'LIQA':
            args = worker_id,theta_df,isoform_df,read_len_dist,EM_choice,output_path
        elif EM_choice == 'LIQA_modified':
            args = worker_id,theta_df,isoform_df,Sm_dict,read_len_dist,EM_choice,output_path
        futures.append(pool.apply_async(E_step_MT,(args,)))
    list_of_isoform_q_df = []
    for future in futures:
        new_isoform_q_df = future.get()
        list_of_isoform_q_df.append(new_isoform_q_df)
    pool.close()
    pool.join()
    return list_of_isoform_q_df
def M_step(list_of_isoform_q_df,theta_df):
    isoform_q_df = pd.Series(0,index=theta_df.index)
    for single_thread_isoform_q_df in list_of_isoform_q_df:
        isoform_q_df = isoform_q_df.add(single_thread_isoform_q_df,fill_value=0)
    new_theta_df = isoform_q_df/isoform_q_df.sum()
    return new_theta_df
def EM_algo(threads,theta_df,isoform_len_df,read_len_dist,output_path,EM_choice):
    num_iters = 2000
    min_diff = 1e-3
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    for i in range(num_iters):
        isoform_df = pd.DataFrame({'isoform_len':isoform_len_df,'theta':theta_df})
        isoform_df = isoform_df.fillna(0)
        Sm_dict = None
        if EM_choice == 'LIQA_modified':
            Sm_dict = get_Sm_dict(read_len_dist,isoform_df)
        list_of_isoform_q_df = E_step(theta_df,threads,output_path,EM_choice,isoform_df,read_len_dist,Sm_dict)
        new_theta_df = M_step(list_of_isoform_q_df,theta_df)
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
    return theta_df
def get_isoform_df(isoform_len_df,expression_dict):
    theta_df = pd.Series(expression_dict) 
    theta_df = theta_df/theta_df.sum()
    isoform_df = pd.DataFrame({'isoform_len':isoform_len_df,'theta':theta_df})
    isoform_df = isoform_df.fillna(0)
    return isoform_df,theta_df
def EM_algo_theta_iter_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_mapping,LR_gene_regions_dict,threads,output_path,EM_choice):
    _,reads_isoform_info,_,expression_dict,_ = get_reads_isoform_info(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_mapping,LR_gene_regions_dict)
    read_len_dict,read_len_dist,isoform_len_df = get_read_len_dist(reads_isoform_info,isoform_len_dict)
    prepare_MT(reads_isoform_info,read_len_dict,threads,output_path)
    del reads_isoform_info
    del read_len_dict
    del gene_regions_read_mapping
    del LR_gene_regions_dict
    del isoform_exon_dict
    del isoform_len_dict
    del strand_dict
    _,theta_df = get_isoform_df(isoform_len_df,expression_dict)
    del expression_dict
    gc.collect()
    final_theta_df = EM_algo(threads,theta_df,isoform_len_df,read_len_dist,output_path,EM_choice)


    # if EM_choice == 'LIQA_modified':
    #     Sm_dict = get_Sm_dict(read_len_dist,isoform_df)
    #     get_cond_prob_MT_LIQA_modified(threads,output_path,isoform_df,read_len_dist,Sm_dict)
    # elif EM_choice == "LIQA":
    #     get_cond_prob_MT_LIQA(threads,output_path,isoform_df,read_len_dist)
    # final_theta_df = EM_algo(threads,theta_df,output_path)
    TPM_df = (final_theta_df/final_theta_df.sum())*1e6
    TPM_df.name = 'TPM'
    TPM_df.index.name = 'Isoform'
    TPM_df.to_csv(f'{output_path}/EM_expression.out',sep='\t')