import pysam
from pathlib import Path
import pandas as pd
import dill as pickle
import numpy as np
import multiprocessing as mp
def get_hits_dict(args):
    worker_id,alignment_file_path,start_pos,num_reads,output_path = args
    hits_dict = {}
    all_fragment_lengths = []
    num_reads_processed = 0
    previous_read_name = None
    with pysam.AlignmentFile(alignment_file_path, "r") as f:
        f.seek(start_pos)
        for read in f:
            if previous_read_name is None:
                previous_read_name = read.query_name
            elif previous_read_name != read.query_name:
                previous_read_name = read.query_name
                num_reads_processed += 1
            if num_reads_processed > num_reads:
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
    num_reads_df = pd.Series(num_reads_dict)
    with open(f'{output_path}/temp/hits_dict/{worker_id}','wb') as f:
        pickle.dump(hits_dict,f)
    print(f'Worker {worker_id} Done!',flush=True)
    return all_fragment_lengths,num_reads_df
def get_aln_line_marker(alignment_file_path,threads):
    previous_read_name = None
    with pysam.AlignmentFile(alignment_file_path, "r") as f:
        line_offset = []
        for line in f:
            read_name = line.query_name
            if previous_read_name is None or read_name != previous_read_name:
                line_offset.append(f.tell())
                previous_read_name = read_name
    num_reads = len(line_offset)
    chunksize, extra = divmod(num_reads, threads)
    if extra:
        chunksize += 1
    aln_line_marker = []
    for i in range(threads):
        aln_line_marker.append((line_offset[i*chunksize],chunksize,i))
    return aln_line_marker
def get_all_hits_dict(eff_len_dict,alignment_file_path,aln_line_marker,threads,output_path):
    pool = mp.Pool(threads)
    futures = []
    for i in range(threads):
        start_pos,num_reads,worker_id = aln_line_marker[i]
        args = worker_id,alignment_file_path,start_pos,num_reads,output_path
        futures.append(pool.apply_async(get_hits_dict,(args,)))
    all_fragment_lengths = []
    list_of_num_reads_df = []
    for future in futures:
        single_thread_fragment_lengths,num_reads_df =future.get()
        all_fragment_lengths += single_thread_fragment_lengths
        list_of_num_reads_df.append(num_reads_df) 
    pool.close()
    pool.join()
    mean_f_len = np.mean(all_fragment_lengths)
    std_f_len = np.std(all_fragment_lengths)
    eff_len_df = pd.Series(eff_len_dict)
    num_reads_df = pd.Series(0,index=eff_len_df.index)
    for single_thread_num_reads_df in list_of_num_reads_df:
        num_reads_df = num_reads_df.add(single_thread_num_reads_df,fill_value=0)
    theta_df = num_reads_df/eff_len_df
    theta_df = theta_df/theta_df.sum()
    return theta_df,mean_f_len,std_f_len
def get_ant(fragment_length,eff_length,mean_f_len,std_f_len):
    # assert eff_length - fragment_length + 1 > 0
    ant = 1/(std_f_len * np.sqrt(np.pi)) * np.e**(-1/2*((fragment_length-mean_f_len)/std_f_len)**2) / eff_length
    return ant
def get_all_ant(args):
    worker_id,eff_len_dict,mean_f_len,std_f_len,output_path = args
    with open(f'{output_path}/temp/hits_dict/{worker_id}','rb') as f:
        hits_dict = pickle.load(f)
    for read in hits_dict:
        for hit in hits_dict[read]:
            eff_length = eff_len_dict[hit['isoform']]
            hit['ant'] = get_ant(hit['fragment_length'],eff_length,mean_f_len,std_f_len)
    with open(f'{output_path}/temp/hits_dict/{worker_id}','wb') as f:
        pickle.dump(hits_dict,f)
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
def prepare_hits(SR_sam,eff_len_dict,output_path,threads):
    Path(f'{output_path}/temp/').mkdir(exist_ok=True,parents=True)
    print('Sorting sam file by read name...',flush=True)
    pysam.sort(SR_sam,'-n','-@',str(threads),'-o',f'{output_path}/temp/SR.sam')
    alignment_file_path = f'{output_path}/temp/SR.sam'
    print('Done',flush=True)
    print('Getting short reads info...',flush=True)
    aln_line_marker = get_aln_line_marker(alignment_file_path,threads)
    Path(f'{output_path}/temp/hits_dict/').mkdir(exist_ok=True,parents=True)
    theta_df,mean_f_len,std_f_len =  get_all_hits_dict(eff_len_dict,alignment_file_path,aln_line_marker,threads,output_path)
    get_ant_all_workers(eff_len_dict,mean_f_len,std_f_len,threads,output_path)
    print('Done',flush=True)
    return theta_df


    
