import dill as pickle
from pathlib import Path
def prepare_MT(reads_isoform_info,read_len_dict,threads,output_path):
    Path(f'{output_path}/temp/').mkdir(exist_ok=True,parents=True)
    chunksize, extra = divmod(len(reads_isoform_info), threads)
    if extra:
        chunksize += 1
    single_thread_reads_isoform_info = {}
    single_thread_read_len_dict = {}
    ct = 0
    worker_id = 0
    for read in reads_isoform_info:
        single_thread_reads_isoform_info[read] = reads_isoform_info[read]
        single_thread_read_len_dict[read] = read_len_dict[read]
        ct += 1
        if ct == chunksize:
            with open(f'{output_path}/temp/{worker_id}','wb') as f:
                pickle.dump((single_thread_reads_isoform_info,single_thread_read_len_dict),f)
            worker_id += 1
            ct = 0
            single_thread_reads_isoform_info = {}
            single_thread_read_len_dict = {}
    if ct != 0:
        with open(f'{output_path}/temp/{worker_id}','wb') as f:
            pickle.dump((single_thread_reads_isoform_info,single_thread_read_len_dict),f)