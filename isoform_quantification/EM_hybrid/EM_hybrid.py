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
from EM_hybrid.EM_SR import prepare_hits
import scipy
import config
import datetime
def E_step_SR(ANT,theta_eff_len_product_arr):
    q = ANT.multiply(theta_eff_len_product_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def E_step_LR(cond_prob,theta_arr):
    q = cond_prob.multiply(theta_arr)
    q_sum = q.sum(axis=1)
    q_sum[q_sum == 0] = 1
    isoform_q_arr = q.multiply(1/q_sum).sum(axis=0)
    return isoform_q_arr
def E_step_MT(args):
    os.nice(10)
#     MIN_PROB = 1e-100
    worker_id,pipe,queue,output_path,SR_num_batches,LR_num_batches,eff_len_arr = args
    # all_ANT_paths = [f'{output_path}/temp/hits_dict/{worker_id}_{batch_id}_ANT.npz' for batch_id in range(SR_num_batches)]
    # all_cond_prob_paths = [f'{output_path}/temp/cond_prob/{worker_id}_{batch_id}_cond_prob.npz' for batch_id in range(LR_num_batches)]
    # ANT = scipy.sparse.vstack([scipy.sparse.load_npz(fpath) for fpath in all_ANT_paths])
    ANT = scipy.sparse.load_npz(f'{output_path}/temp/hits_dict/{worker_id}_ANT.npz')
    # cond_prob = scipy.sparse.vstack([scipy.sparse.load_npz(fpath) for fpath in all_cond_prob_paths])
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
            
            isoform_q_arr_SR = E_step_SR(ANT,theta_eff_len_product_arr)
            out_fpath = f'{output_path}/temp/EM_isoform_q_SR/{worker_id}.npy'
            # np.savez_compressed(out_fpath,q=isoform_q_arr_SR)
            mm = np.memmap(out_fpath, dtype=isoform_q_arr_SR.dtype, shape=isoform_q_arr_SR.shape, mode='w+')
            mm[:] = isoform_q_arr_SR
            queue.put(('SR',out_fpath,isoform_q_arr_SR.dtype, isoform_q_arr_SR.shape))
            queue.put('done')
    return
def M_step(isoform_q_arr_SR_all,isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha,alpha_df):
    ss = eff_len_arr / ((theta_arr * eff_len_arr).sum())
    # if alpha_df is None:
    new_theta_arr = ((1-alpha) * isoform_q_arr_SR_all + alpha * isoform_q_arr_LR_all) / ((1-alpha) * isoform_q_arr_SR_all.sum() * ss + alpha * isoform_q_arr_LR_all.sum())
    # else:
    #     new_theta_arr = ((1-alpha_df) * isoform_q_arr_SR_all + alpha_df * isoform_q_arr_LR_all) / ((1-alpha_df) * isoform_q_df_SR.sum() * ss + alpha_df * isoform_q_arr_LR_all.sum())
    new_theta_arr = new_theta_arr/new_theta_arr.sum()
    new_theta_arr[new_theta_arr<1e-100] = 0
    return new_theta_arr
def EM_listener(watcher_args):
    threads,eff_len_arr,theta_arr,output_path,queue,all_pipes,output_df = watcher_args
    min_diff = 1e-3
    num_iters = config.EM_SR_num_iters
    Path(f'{output_path}/temp/EM_isoform_q_SR/').mkdir(exist_ok=True,parents=True)
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
        isoform_q_arr_SR_all = None
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

                # with np.load(new_isoform_q_arr_path) as f:
                #     new_isoform_q_arr = f['q'] 
                if seq == 'LR':
                    if isoform_q_arr_LR_all is None:
                        isoform_q_arr_LR_all = new_isoform_q_arr
                    else:
                        isoform_q_arr_LR_all = isoform_q_arr_LR_all + new_isoform_q_arr
                elif seq == 'SR':
                    if isoform_q_arr_SR_all is None:
                        isoform_q_arr_SR_all = new_isoform_q_arr
                    else:
                        isoform_q_arr_SR_all = isoform_q_arr_SR_all + new_isoform_q_arr
        new_theta_arr = M_step(isoform_q_arr_SR_all,isoform_q_arr_LR_all,theta_arr,eff_len_arr,alpha,alpha_df)
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
            # if i % 100 == 0:
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
def EM_manager(threads,eff_len_arr,theta_arr,output_path,SR_num_batches_dict,LR_num_batches_dict,output_df):
    pool = mp.Pool(threads+1)
    manager = mp.Manager()
    queue = manager.Queue()
    futures = []
    all_pipes = []
    for worker_id in range(threads):
        pipe = mp.Pipe()
        args = worker_id,pipe,queue,output_path,SR_num_batches_dict[worker_id],LR_num_batches_dict[worker_id],eff_len_arr
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
def EM_algo_hybrid(isoform_len_dict,SR_sam,output_path,threads,EM_choice):
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
    # prepare SR
    theta_SR_arr,eff_len_arr,SR_num_batches_dict = prepare_hits(SR_sam,output_path,isoform_index_dict,threads)
    theta_LR_arr,_,LR_num_batches_dict = prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path)
    num_SRs = theta_SR_arr.sum()
    num_LRs = theta_LR_arr.sum()
    print(f'Number of SRs/eff_len:{num_SRs}')
    print(f'Number of LRs:{num_LRs}')
    print(f'Pseudo_count_SR:'+str(config.pseudo_count_SR))
    print(f'Pseudo_count_LR:'+str(config.pseudo_count_LR),flush=True)
    write_result_to_tsv(f'{output_path}/SR_count.tsv',output_df,theta_SR_arr.flatten())
    write_result_to_tsv(f'{output_path}/LR_count.tsv',output_df,theta_LR_arr.flatten())
    if config.inital_theta in ['hybrid_unique','hybrid']:
        # theta_SR_arr /= theta_SR_arr.sum()
        # theta_LR_arr /= theta_LR_arr.sum()
        if config.eps_strategy == 'add_eps':
            theta_SR_arr += config.pseudo_count_SR
            theta_LR_arr += config.pseudo_count_LR
            theta_SR_arr = theta_SR_arr / theta_SR_arr.sum()
            theta_LR_arr = theta_LR_arr / theta_LR_arr.sum()
            theta_arr = (theta_SR_arr+theta_LR_arr)/2
        elif config.eps_strategy == 'add_eps_small':
            theta_arr = (theta_SR_arr/theta_SR_arr.sum() + theta_LR_arr/theta_LR_arr.sum())/2
            theta_SR_arr[theta_arr<config.inital_theta_eps] += config.pseudo_count_SR
            theta_LR_arr[theta_arr<config.inital_theta_eps] += config.pseudo_count_LR
            theta_SR_arr = theta_SR_arr / theta_SR_arr.sum()
            theta_LR_arr = theta_LR_arr / theta_LR_arr.sum()
            theta_arr = (theta_SR_arr+theta_LR_arr)/2
    else:
        # theta_SR_arr /= theta_SR_arr.sum()
        # theta_LR_arr /= theta_LR_arr.sum()
        if config.eps_strategy == 'add_eps_small':
            theta_SR_arr[theta_SR_arr/theta_SR_arr.sum()<config.inital_theta_eps] += config.pseudo_count_SR
            theta_LR_arr[theta_LR_arr/theta_LR_arr.sum()<config.inital_theta_eps] += config.pseudo_count_LR
        elif config.eps_strategy == 'add_eps':
            theta_SR_arr += config.pseudo_count_SR
            theta_LR_arr += config.pseudo_count_LR
        theta_SR_arr = theta_SR_arr / theta_SR_arr.sum()
        theta_LR_arr = theta_LR_arr / theta_LR_arr.sum()
        if config.inital_theta in ['SR','SR_unique']:
            theta_arr = theta_SR_arr
        elif config.inital_theta in ['LR','LR_unique']:
            theta_arr = theta_LR_arr
        elif config.inital_theta == 'uniform':
            theta_arr = np.ones(shape=theta_SR_arr.shape)
            theta_arr = theta_arr/theta_arr.sum()
        elif config.inital_theta == 'random':
            theta_arr = np.random.dirichlet(np.ones(shape=theta_SR_arr.flatten().shape))
            if config.eps_strategy == 'add_eps_small':
                count_arr = theta_arr * num_LRs
                count_arr[theta_arr<config.inital_theta_eps] += config.pseudo_count_LR
                theta_arr = count_arr
            elif config.eps_strategy == 'add_eps':
                theta_arr = theta_arr * num_LRs
                theta_arr += config.inital_theta_eps
            theta_arr = np.expand_dims(theta_arr,0)
            theta_arr = theta_arr/theta_arr.sum()
    # np.savez_compressed(f'{output_path}/initial_theta',theta=theta_arr)
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    write_result_to_tsv(f'{output_path}/initial_theta.tsv',output_df,theta_arr.flatten())
    TPM_arr = (theta_arr.flatten()/theta_arr.flatten().sum())*1e6
    write_result_to_tsv(f'{output_path}/initial_theta_TPM.tsv',output_df,TPM_arr)
    print('Using {} as initial theta'.format(config.inital_theta))
    
    # theta_arr = theta_arr/theta_arr.sum()
    EM_manager(threads,eff_len_arr,theta_arr,output_path,SR_num_batches_dict,LR_num_batches_dict,output_df)
    
