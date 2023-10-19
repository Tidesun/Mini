import multiprocessing as mp
from EM_hybrid.get_reads_isoform_cond_prob import get_cond_prob_MT_LIQA_modified,get_Sm_dict
from EM_hybrid.get_reads_isoform_info import get_reads_isoform_info,get_read_len_dist
from EM_hybrid.prepare_MT import prepare_MT
import pickle
import pandas as pd
import multiprocessing as mp
import numpy as np
from pathlib import Path
import gc
import glob
import shutil

def get_isoform_df(isoform_len_df,expression_dict):
    theta_df = pd.Series(expression_dict) 
    theta_df = theta_df/theta_df.sum()
    isoform_df = pd.DataFrame({'isoform_len':isoform_len_df,'theta':theta_df})
    isoform_df = isoform_df.fillna(0)
    return isoform_df,theta_df
def prepare_LR(isoform_len_df,isoform_index_dict,isoform_index_series,threads,output_path):
    expression_dict = {}
    read_len_dist_dict = {}
    for fpath in glob.glob(f'{output_path}/temp/LR_alignments/dist_*_*'):
        with open(fpath,'rb') as f:
            [new_read_len_dist_dict,new_expression_dict] = pickle.load(f)
        for isoform in new_expression_dict:
            if isoform not in expression_dict:
                expression_dict[isoform] = 0
            expression_dict[isoform] += new_expression_dict[isoform]
        for read_len in new_read_len_dist_dict:
            if read_len not in read_len_dist_dict:
                read_len_dist_dict[read_len] = 0
            read_len_dist_dict[read_len] += new_read_len_dist_dict[read_len]
    read_len_dist = pd.Series(read_len_dist_dict).to_frame()
    read_len_dist.columns = ['PDF']
    read_len_dist = read_len_dist.sort_index(ascending=True)
    isoform_df,theta_df = get_isoform_df(isoform_len_df,expression_dict)
    del isoform_len_df
    del expression_dict
    gc.collect()
    Sm_dict = get_Sm_dict(read_len_dist,isoform_df)
<<<<<<< HEAD
    get_cond_prob_MT_LIQA_modified(threads,output_path,isoform_df,read_len_dist,Sm_dict)
    return theta_df,isoform_df
=======
    num_batches_dict = get_cond_prob_MT_LIQA_modified(threads,output_path,isoform_df,isoform_index_dict,read_len_dist,Sm_dict)
    theta_arr = np.expand_dims(pd.concat([isoform_index_series, theta_df], axis=1).fillna(0).sort_values('Index')[0].values,axis=0)
    return theta_arr,isoform_df,num_batches_dict
>>>>>>> 24929069bf9997b145e21733b499aeb1f08cef25
# def EM_algo_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_mapping,LR_gene_regions_dict,threads,output_path,EM_choice):
#     theta_df,isoform_df = prepare_LR(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_mapping,LR_gene_regions_dict,threads,output_path,EM_choice)
#     final_theta_df = EM_algo(threads,theta_df,output_path)
#     TPM_df = (final_theta_df/final_theta_df.sum())*1e6
#     TPM_df.name = 'TPM'
#     TPM_df.index.name = 'Isoform'
#     TPM_df = TPM_df.to_frame().join(isoform_df).fillna(0)
#     TPM_df = TPM_df['TPM']
#     TPM_df.to_csv(f'{output_path}/EM_expression.out',sep='\t')
