# hybrid
import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from pathlib import Path
import glob
import os
import gc
from EM_hybrid.EM_LR import prepare_LR
import scipy
import config
import datetime
def E_step_LR(cond_prob,theta_arr):
    q = cond_prob.multiply(theta_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def E_step_MT(args):
    os.nice(10)
#     MIN_PROB = 1e-100
    worker_id,pipe,queue,output_path,LR_num_batches,eff_len_arr = args
    cond_prob = scipy.sparse.load_npz(f'{output_path}/temp/cond_prob/{worker_id}_cond_prob.npz')
    p_out,p_in = pipe
#     p_out.close()
    while True:
        msg = p_out.recv()
        if msg == 'kill':
            break
        else:
            # new iteration
            theta_arr_fpath,dtype,shape = msg
            # with np.load(theta_arr_fpath) as f:
            #     theta_arr = f['theta']
            theta_arr = np.memmap(theta_arr_fpath, dtype=dtype, shape=shape, mode='c')
            theta_eff_len_product_arr = theta_arr * eff_len_arr /((theta_arr * eff_len_arr).sum())
            isoform_q_arr_LR = E_step_LR(cond_prob,theta_arr)
            
            
            out_fpath = f'{output_path}/temp/EM_isoform_q_LR/{worker_id}.npy'
            mm = np.memmap(out_fpath, dtype=isoform_q_arr_LR.dtype, shape=isoform_q_arr_LR.shape, mode='w+')
            mm[:] = isoform_q_arr_LR
            # np.savez_compressed(out_fpath,q=isoform_q_arr_LR)
            queue.put(('LR',out_fpath,isoform_q_arr_LR.dtype, isoform_q_arr_LR.shape))
            
            queue.put('done')
    return
def M_step(isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha,alpha_df):
    # if alpha_df is None:
    new_theta_arr = (isoform_q_arr_LR_all) / (isoform_q_arr_LR_all.sum())
    # else:
    #     new_theta_arr = ((1-alpha_df) * isoform_q_arr_SR_all + alpha_df * isoform_q_arr_LR_all) / ((1-alpha_df) * isoform_q_df_SR.sum() * ss + alpha_df * isoform_q_arr_LR_all.sum())
    new_theta_arr = new_theta_arr/new_theta_arr.sum()
    new_theta_arr[new_theta_arr<1e-100] = 0
    return new_theta_arr
def EM_listener(watcher_args):
    threads,eff_len_arr,theta_arr,output_path,queue,all_pipes,output_df = watcher_args
    min_diff = 1e-3
    num_iters = config.EM_SR_num_iters
    Path(f'{output_path}/temp/EM_isoform_q_LR/').mkdir(exist_ok=True,parents=True)
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    # with open(f'{output_path}/EM_log.txt','w') as logger:
    alpha_df_path = config.alpha_df_path
    if alpha_df_path is not None:
        alpha_df = pd.read_csv(alpha_df_path,sep='\t',skiprows=1,header=None)
        alpha_df.columns = ['Isoform','Alpha']
        alpha_df = alpha_df.set_index('Isoform')['Alpha']
        alpha = None
    else:
        alpha = float(config.alpha)
        alpha_df = None
    for i in range(num_iters):
        # define isoform_df
        # isoform_df = pd.DataFrame({'eff_len':eff_len_dict,'theta':theta_df})
        # isoform_df = isoform_df.fillna(0)
        # isoform_df['len_theta_product'] = isoform_df['eff_len'] * isoform_df['theta'] / ((isoform_df['eff_len'] * isoform_df['theta']).sum())
        theta_arr_fpath = f'{output_path}/temp/EM_theta_arr.npy'
        mm = np.memmap(theta_arr_fpath, dtype=theta_arr.dtype, shape=theta_arr.shape, mode='w+')
        mm[:] = theta_arr
        # with open(isoform_df_fpath,'wb') as f:
        #     pickle.dump(isoform_df,f)
        # np.savez_compressed(theta_arr_fpath,theta=theta_arr)
        for (p_out,p_in) in all_pipes:
            p_in.send((theta_arr_fpath,theta_arr.dtype,theta_arr.shape))
#             print('Send isoform df')
        isoform_q_arr_LR_all = None
        num_workers_done = 0
        while True:
            msg = queue.get()
            if msg == 'done':
                num_workers_done += 1
                if num_workers_done >= threads:
                    break
            else:
                seq,new_isoform_q_arr_path,dtype,shape = msg
                new_isoform_q_arr = np.memmap(new_isoform_q_arr_path, dtype=dtype, shape=shape, mode='c')

                if isoform_q_arr_LR_all is None:
                    isoform_q_arr_LR_all = new_isoform_q_arr
                else:
                    isoform_q_arr_LR_all = isoform_q_arr_LR_all + new_isoform_q_arr
        new_theta_arr = M_step(isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha,alpha_df)
        diff = np.abs(theta_arr[new_theta_arr>1e-7] - new_theta_arr[new_theta_arr>1e-7])/new_theta_arr[new_theta_arr>1e-7]
        # diff = diff[new_theta_arr>1e-7]
        # diff[(new_theta_arr==0).flatten()] = 0
        if i % config.EM_output_frequency == 0:
            print(f'Iteration:{i}')
            print(datetime.datetime.now())
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max(),flush=True)
            TPM_arr = (new_theta_arr.flatten()/new_theta_arr.flatten().sum())*1e6
            write_result_to_tsv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',output_df,TPM_arr)
            #     theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')

        if diff[diff > min_diff].shape[0] == 0:
            print(f'Iteration:{i}')
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max())
            print('Converge!',flush=True)
            TPM_arr = (new_theta_arr.flatten()/new_theta_arr.flatten().sum())*1e6
            write_result_to_tsv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',output_df,TPM_arr)
            # theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')
            break
        theta_arr = new_theta_arr.copy()
    # finished EM
    for (p_out,p_in) in all_pipes:
        p_in.send('kill')
    TPM_arr = (new_theta_arr.flatten()/new_theta_arr.flatten().sum())*1e6
    write_result_to_tsv(f'{output_path}/EM_expression.out',output_df,TPM_arr)
def callback_error(result):
    print('ERR:', result,flush=True)
def EM_manager(threads,eff_len_arr,theta_arr,output_path,LR_num_batches_dict,output_df):
    pool = mp.Pool(threads+1)
    manager = mp.Manager()
    queue = manager.Queue()
    futures = []
    all_pipes = []
    for worker_id in range(threads):
        pipe = mp.Pipe()
        args = worker_id,pipe,queue,output_path,LR_num_batches_dict[worker_id],eff_len_arr
        futures.append(pool.apply_async(E_step_MT,(args,),error_callback=callback_error))
        all_pipes.append(pipe)
    watcher_args = threads,eff_len_arr,theta_arr,output_path,queue,all_pipes,output_df
    watcher = pool.apply_async(EM_listener,(watcher_args,),error_callback=callback_error)
    for future in futures:
        future.get()
    watcher.get()
    for (p_out,p_in) in all_pipes:
        p_in.close()
        p_out.close()
    pool.close()
    pool.join()
def write_result_to_tsv(result_path,output_df,TPM_arr):
    output_df['TPM'] = TPM_arr
    output_df = output_df[['Isoform','TPM','Effective_length']]
    output_df.to_csv(result_path,sep='\t',index=False)
def EM_algo_LR_alone(isoform_len_dict,SR_sam,output_path,threads,EM_choice):
    # prepare arr
    isoform_len_df = pd.Series(isoform_len_dict)
    isoform_list = sorted(isoform_len_dict.keys())
    isoform_index_dict = {}
    isoform_len_arr = []
    for i,isoform in enumerate(isoform_list):
        isoform_index_dict[isoform] = i
        isoform_len_arr.append(isoform_len_dict[isoform])
    isoform_index_series = pd.Series(isoform_index_dict)
    isoform_index_series.name = 'Index'
    isoform_index_series.index.name = 'Isoform'
    output_df = isoform_index_series.to_frame().reset_index().sort_values(by="Index")
    isoform_len_arr = np.array(isoform_len_arr)
    eff_len_arr = isoform_len_arr.copy()
    output_df['Effective_length'] = eff_len_arr
    theta_LR_arr,_,LR_num_batches_dict = prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path)
    num_LRs = theta_LR_arr.sum()
    print(f'Number of LRs:{num_LRs}')
    print(f'Pseudo_count_LR:'+str(config.pseudo_count_LR),flush=True)
    if config.eps_strategy == 'add_eps_small':
        theta_LR_arr[theta_LR_arr<config.inital_theta_eps] += config.pseudo_count_LR
    elif config.eps_strategy == 'add_eps':
        theta_LR_arr += config.pseudo_count_LR
    theta_LR_arr = theta_LR_arr / theta_LR_arr.sum()
    if config.inital_theta in ['LR','LR_unique']:
        theta_arr = theta_LR_arr
    elif config.inital_theta == 'uniform':
        theta_arr = np.ones(shape=theta_LR_arr.shape)
        theta_arr = theta_arr/theta_arr.sum()
    elif config.inital_theta == 'random':
        theta_arr = np.random.dirichlet(np.ones(shape=theta_LR_arr.flatten().shape))
        if config.eps_strategy == 'add_eps_small':
            count_arr = theta_arr * num_LRs
            count_arr[theta_arr<config.inital_theta_eps] += config.pseudo_count_LR
            theta_arr = count_arr
        elif config.eps_strategy == 'add_eps':
            theta_arr = theta_arr * num_LRs
            theta_arr += config.inital_theta_eps
        theta_arr = np.expand_dims(theta_arr,0)
        theta_arr = theta_arr/theta_arr.sum()
    np.savez_compressed(f'{output_path}/initial_theta',theta=theta_arr)
    print('Using {} as initial theta'.format(config.inital_theta))
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    theta_arr = theta_arr/theta_arr.sum()
    EM_manager(threads,eff_len_arr,theta_arr,output_path,LR_num_batches_dict,output_df)
    
