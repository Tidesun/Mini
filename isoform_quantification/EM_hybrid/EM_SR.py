from pathlib import Path
import pandas as pd
import pickle
import numpy as np
import multiprocessing as mp
import glob
import os
from EM_hybrid.prepare_hits import prepare_hits
import config
def E_step_cal(hits_df,isoform_df):
    hits_df['len_theta_product'] = isoform_df.loc[hits_df.index.droplevel(0),['len_theta_product']].set_index(hits_df.index)
    q = hits_df['ant'] * hits_df['len_theta_product']
    q_sum = q.groupby('read').sum()
    q_sum[q_sum == 0] = 1
    q = q/q_sum
    isoform_q_df = q.groupby('isoform').sum()
    return isoform_q_df
def E_step_MT(args):
    os.nice(10)
#     MIN_PROB = 1e-100
    worker_id,pipe,queue,output_path = args
    all_hits_dict_paths = list(glob.glob(f'{output_path}/temp/hits_dict/{worker_id}_*'))
    p_out,p_in = pipe
#     p_out.close()
    while True:
        msg = p_out.recv()
        if msg == 'kill':
            break
        else:
            isoform_df_fpath = msg
            with open(isoform_df_fpath,'rb') as f:
                isoform_df = pickle.load(f)
            for fpath in all_hits_dict_paths:
                batch_id = fpath.split('/')[-1].split('_')[1]
                with open(fpath,'rb') as f:
                    hits_df = pickle.load(f)
                isoform_q_df = E_step_cal(hits_df,isoform_df)
                del hits_df
                out_fpath = f'{output_path}/temp/EM_isoform_q/{worker_id}_{batch_id}'
                with open(out_fpath,'wb') as f:
                    pickle.dump(isoform_q_df,f)
                queue.put(out_fpath)
            queue.put('done')
    return
def M_step(isoform_q_df,isoform_df):
    ss = isoform_df['eff_len'] / ((isoform_df['eff_len']*isoform_df['theta']).sum())
    new_theta_df = isoform_q_df / (isoform_q_df.sum() * ss)
    new_theta_df = new_theta_df/new_theta_df.sum()
    return new_theta_df
def EM_listener(watcher_args):
    threads,eff_len_dict,theta_df,output_path,queue,all_pipes = watcher_args
    min_diff = 1e-3
    num_iters = config.EM_SR_num_iters
    Path(f'{output_path}/temp/EM_isoform_q/').mkdir(exist_ok=True,parents=True)
    # with open(f'{output_path}/EM_log.txt','w') as logger:
    for i in range(num_iters):
        # define isoform_df
        isoform_df = pd.DataFrame({'eff_len':eff_len_dict,'theta':theta_df})
        isoform_df = isoform_df.fillna(0)
        isoform_df['len_theta_product'] = isoform_df['eff_len'] * isoform_df['theta'] / ((isoform_df['eff_len'] * isoform_df['theta']).sum())
        isoform_df_fpath = f'{output_path}/temp/EM_isoform_df'
        with open(isoform_df_fpath,'wb') as f:
            pickle.dump(isoform_df,f)
        for (p_out,p_in) in all_pipes:
            p_in.send(isoform_df_fpath)
#             print('Send isoform df')
        isoform_q_df = pd.Series(0,index=isoform_df.index)
        num_workers_done = 0
        while True:
            msg = queue.get()
            if msg == 'done':
                num_workers_done += 1
            else:
                isoform_q_df_path = msg
                with open(isoform_q_df_path,'rb') as f:
                    new_isoform_q_df = pickle.load(f)
                isoform_q_df = isoform_q_df.add(new_isoform_q_df,fill_value=0)
            if num_workers_done == threads:
                break

        new_theta_df = M_step(isoform_q_df,isoform_df)
        diff = np.abs(theta_df - new_theta_df)/new_theta_df
        diff = diff[new_theta_df>1e-7]
        diff[new_theta_df==0] = 0
        if i % 1 == 0:
            print(f'Iteration:{i}')
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max(),flush=True)
            theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')

        if diff[diff > min_diff].shape[0] == 0:
            print(f'Iteration:{i}')
        #     print('Sum:')
        #     print(Q.sum())
            print('bChange:')
            print(diff.max(),flush=True)
            theta_df.to_csv(f'{output_path}/EM_iterations/Iter_{i}_theta.tsv',sep='\t')
            break
        theta_df = new_theta_df.copy()
    # finished EM
    for (p_out,p_in) in all_pipes:
        p_in.send('kill')
    return theta_df
def callback_error(result):
    print('ERR:', result)
def EM_manager(threads,eff_len_dict,theta_df,output_path):
    pool = mp.Pool(threads+1)
    manager = mp.Manager()
    queue = manager.Queue()
    futures = []
    all_pipes = []
    for worker_id in range(threads):
        pipe = mp.Pipe()
        args = worker_id,pipe,queue,output_path
        futures.append(pool.apply_async(E_step_MT,(args,),error_callback=callback_error))
        all_pipes.append(pipe)
    watcher_args = threads,eff_len_dict,theta_df,output_path,queue,all_pipes
    watcher = pool.apply_async(EM_listener,(watcher_args,),error_callback=callback_error)
    for future in futures:
        future.get()
    theta_df = watcher.get()
    for (p_out,p_in) in all_pipes:
        p_in.close()
        p_out.close()
    pool.close()
    pool.join()
    return theta_df
def EM_algo_SR(SR_sam,output_path,threads):
    theta_df,eff_len_dict = prepare_hits(SR_sam,output_path,threads)
    Path(f'{output_path}/EM_iterations/').mkdir(exist_ok=True,parents=True)
    theta_df = EM_manager(threads,eff_len_dict,theta_df,output_path)
    TPM_df = (theta_df/theta_df.sum())*1e6
    eff_len_df = pd.Series(eff_len_dict)
    out_df = pd.DataFrame({'TPM':TPM_df,'Effective_length':eff_len_df})
    out_df.index.name = 'Isoform'
    out_df.to_csv(f'{output_path}/EM_expression.out',sep='\t')