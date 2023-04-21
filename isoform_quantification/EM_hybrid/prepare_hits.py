import pysam
from pathlib import Path
import pandas as pd
import pickle
import numpy as np
import multiprocessing as mp
import glob
import config
import shutil
from EM_hybrid.cal_eff_len import get_eff_len_dict
def dump_hits_dict(all_fragment_lengths,hits_dict,worker_id,batch_id,output_path):
    num_reads_dict = {}
    for read in hits_dict:
        if len(hits_dict[read]) == 1:
            if hits_dict[read][0]['isoform'] not in num_reads_dict:
                num_reads_dict[hits_dict[read][0]['isoform']] = 0
            num_reads_dict[hits_dict[read][0]['isoform']] += 1
        else:
            for hit in hits_dict[read]:
                if hit['isoform'] not in num_reads_dict:
                    num_reads_dict[hit['isoform']] = 0
                num_reads_dict[hit['isoform']] += 1/len(hits_dict[read])
    rows = []
    for read in hits_dict:
        for mapping in hits_dict[read]:
            rows.append([read,mapping['isoform'],mapping['fragment_length']])
    hits_df = pd.DataFrame(rows,columns=['read','isoform','fragment_length']).set_index(['read','isoform'])
    with open(f'{output_path}/temp/hits_dict/{worker_id}_{batch_id}','wb') as f:
        pickle.dump(hits_df,f)
    num_reads_df = pd.Series(num_reads_dict)
    with open(f'{output_path}/temp/fragment_lengths/{worker_id}_{batch_id}','wb') as f:
        pickle.dump((all_fragment_lengths,num_reads_df),f)
    frag_lengths_path = f'{output_path}/temp/fragment_lengths/{worker_id}_{batch_id}'
    return frag_lengths_path
def get_hits_dict(args):
    os.nice(10)
    worker_id,alignment_file_path,start_pos,end_pos,output_path,queue = args
    hits_dict = {}
    all_fragment_lengths = []
    previous_read_name = None
    batch_id = 0
    buffer_size = 1e6
    with pysam.AlignmentFile(alignment_file_path, "r") as f:
        f.seek(start_pos)
        for read in f:
            if previous_read_name is None:
                previous_read_name = read.query_name
            elif previous_read_name != read.query_name:
                previous_read_name = read.query_name
                if len(hits_dict) == buffer_size:
                    frag_lengths_path = dump_hits_dict(all_fragment_lengths,hits_dict,worker_id,batch_id,output_path)
                    queue.put(frag_lengths_path)
                    all_fragment_lengths = []
                    hits_dict = {}
                    batch_id += 1
                if f.tell() >= end_pos:
                    break
            if read.template_length > 0:
                if read.query_name not in hits_dict:
                    hits_dict[read.query_name] = []
                if '|' in read.reference_name:
                    isoform_name = read.reference_name.split('|')[0]
                else:
                    isoform_name = read.reference_name
                hits_dict[read.query_name].append({'fragment_length':read.template_length,'isoform':isoform_name})
                all_fragment_lengths.append(read.template_length)
    if len(hits_dict) > 0:
        frag_lengths_path = dump_hits_dict(all_fragment_lengths,hits_dict,worker_id,batch_id,output_path)
        queue.put(frag_lengths_path)
        all_fragment_lengths = []
        hits_dict = {}
        batch_id += 1
    queue.put('done')
    
    print(f'Get_hits_dict: Worker {worker_id} done with {batch_id} batches!',flush=True)
    return
import os
def get_aln_line_marker(alignment_file_path,threads):
    '''
    Split the sam file into THREADS chunks
    !Split by read
    '''
    file_stats = os.stat(alignment_file_path)
    total_bytes = file_stats.st_size
    chunksize, extra = divmod(total_bytes, threads)
    if extra:
        chunksize += 1
    byte_marker = []
    for i in range(threads):
        with open(alignment_file_path,'r') as f:
            if i == 0:
                start_offset = 0
                for line in f:
                    if line[0] != '@':
                        break
                    start_offset += len(line)
            else:
                f.seek(i*chunksize)
                f.readline()
                previous_read_name = None
                previous_byte_pos = f.tell()
                while True:
                    line = f.readline()
                    read_name = line.split('\t')[0]
                    if previous_read_name is None:
                        previous_read_name = read_name
                    elif previous_read_name != read_name:
                        start_offset = previous_byte_pos
                        break
                    previous_byte_pos = f.tell()
            byte_marker.append(start_offset)
    byte_marker.append(total_bytes)
    return byte_marker
def get_all_hits_dict_listener(queue,threads,alignment_file_path):
    all_fragment_lengths = []
    list_of_num_reads_df = []
    num_workers_done = 0
    while True:
        msg = queue.get()
        if msg == 'kill':
            break
        elif msg == 'done':
            num_workers_done += 1
            if num_workers_done == threads:
                mean_f_len = np.mean(all_fragment_lengths)
                std_f_len = np.std(all_fragment_lengths)
                eff_len_dict = get_eff_len_dict(alignment_file_path,mean_f_len,std_f_len,config.eff_len_option)
                eff_len_df = pd.Series(eff_len_dict)
                num_reads_df = pd.Series(0,index=eff_len_df.index)
                for single_thread_num_reads_df in list_of_num_reads_df:
                    num_reads_df = num_reads_df.add(single_thread_num_reads_df,fill_value=0)
                theta_df = num_reads_df/eff_len_df
                theta_df = theta_df/theta_df.sum()
                break
        else:
            frag_lengths_path = msg
            with open(frag_lengths_path,'rb') as f:
                single_thread_fragment_lengths,num_reads_df  = pickle.load(f)
            all_fragment_lengths += single_thread_fragment_lengths
            list_of_num_reads_df.append(num_reads_df)
    return theta_df,mean_f_len,std_f_len,eff_len_dict
def get_all_hits_dict(alignment_file_path,byte_marker,threads,output_path):
    pool = mp.Pool(threads+1)
    manager = mp.Manager()
    queue = manager.Queue()    
    watcher = pool.apply_async(get_all_hits_dict_listener, args=(queue,threads,alignment_file_path))
    futures = []
    for i in range(threads):
        start_pos,end_os = byte_marker[i],byte_marker[i+1]
        args = i,alignment_file_path,start_pos,end_os,output_path,queue
        futures.append(pool.apply_async(get_hits_dict,(args,)))
    for future in futures:
        future.get()
    queue.put('kill')
    theta_df,mean_f_len,std_f_len,eff_len_dict = watcher.get()
    pool.close()
    pool.join()
    return theta_df,mean_f_len,std_f_len,eff_len_dict
def get_ant(row,mean_f_len,std_f_len,eff_len_dict):
    fragment_length = row['fragment_length']
    if row.name[1] in eff_len_dict:
        eff_length = eff_len_dict[row.name[1]]
    else:
        eff_length = 1
    # assert eff_length - fragment_length + 1 > 0
    ant = 1/(std_f_len * np.sqrt(np.pi)) * np.e**(-1/2*((fragment_length-mean_f_len)/std_f_len)**2) / eff_length
    return ant
def get_all_ant(args):
    worker_id,eff_len_dict,mean_f_len,std_f_len,output_path = args
    num_batches_processed = 0
    for fpath in glob.glob(f'{output_path}/temp/hits_dict/{worker_id}_*'):
        with open(fpath,'rb') as f:
            hits_df = pickle.load(f)
        hits_df['ant'] = hits_df.apply(lambda row:get_ant(row,mean_f_len,std_f_len,eff_len_dict),axis=1)
        hits_df = hits_df.drop(columns=['fragment_length'])
        with open(fpath,'wb') as f:
            pickle.dump(hits_df,f)
        num_batches_processed += 1
    print(f'Get_all_ant: Worker {worker_id} done with {num_batches_processed} batches!',flush=True)
def get_ant_all_workers(eff_len_dict,mean_f_len,std_f_len,threads,output_path):
    pool = mp.Pool(threads)
    futures = []
    for i in range(threads):
        args = i,eff_len_dict,mean_f_len,std_f_len,output_path
        futures.append(pool.apply_async(get_all_ant,(args,)))
    for future in futures:
        future.get()
    pool.close()
    pool.join()
def prepare_hits(SR_sam,output_path,threads):
    Path(f'{output_path}/temp/').mkdir(exist_ok=True,parents=True)
    print('Sorting sam file by read name...',flush=True)
    pysam.sort(SR_sam,'-n','-@',str(threads),'-o',f'{output_path}/temp/SR.sam')
    alignment_file_path = f'{output_path}/temp/SR.sam'
    print('Done',flush=True)
    print('Getting short reads info...',flush=True)
    byte_marker = get_aln_line_marker(alignment_file_path,threads)
    Path(f'{output_path}/temp/hits_dict/').mkdir(exist_ok=True,parents=True)
    Path(f'{output_path}/temp/fragment_lengths/').mkdir(exist_ok=True,parents=True)
    theta_df,mean_f_len,std_f_len,eff_len_dict =  get_all_hits_dict(alignment_file_path,byte_marker,threads,output_path)
    print('Get all hits dict done',flush=True)
    get_ant_all_workers(eff_len_dict,mean_f_len,std_f_len,threads,output_path)
    print('Done',flush=True)
    try:
        shutil.rmtree(f'{output_path}/temp/fragment_lengths/')
    except:
        pass
    try:
        Path(alignment_file_path).unlink()
    except:
        pass
    return theta_df,eff_len_dict