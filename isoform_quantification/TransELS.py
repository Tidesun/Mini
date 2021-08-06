from collections import defaultdict
import concurrent.futures
import datetime
import pysam
import time
from itertools import repeat
import dill as pickle
import numpy as np
from numpy import linalg as LA
from qpsolvers import solve_qp
from util import get_long_read_M_dist,get_filtered_out_long_read_M_dist,get_very_short_isoforms
from construct_feature_matrix import generate_all_feature_matrix_short_read,generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment
from parse_alignment_main import parse_alignment
from generate_output import generate_TransELS_output,generate_TrEESR_output
from quantification import quantification
def TransELS(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,region_expression_calculation_method,alpha,beta,P,filtering,threads=1,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=0):
    start_time = time.time()
    print('Preprocessing...')
    LR_gene_read_min_len_dict = None
    # LR_gene_read_min_len_dict = get_long_read_gene_distribution(ref_file_path,long_read_alignment_file_path)
    # print(LR_gene_read_min_len_dict)
    print('Start parsing annoation...')
    if short_read_alignment_file_path is not None:
        with pysam.AlignmentFile(short_read_alignment_file_path, "r") as samfile:
            for read in samfile.head(1):
                READ_LEN = read.infer_query_length()
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,LR_gene_read_min_len_dict)
    gene_regions_points_list,gene_range,gene_interval_tree_dict = process_annotation_for_alignment(gene_exons_dict,gene_points_dict)
    end_time_1 = time.time()
    print('Done in %.3f s'%(end_time_1-start_time))
    print('Mapping short read to regions...')  
    short_read_gene_regions_read_count,SR_read_len,num_SRs = parse_alignment(short_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,gene_isoforms_length_dict,False,False,threads)
    SR_read_len = READ_LEN
    end_time_2 = time.time()
    print('Mapped {} short reads'.format(num_SRs))
    print('Done in %.3f s'%(end_time_2-end_time_1))
    print('Mapping long read to regions...')
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
    end_time_3 = time.time()
    print('Mapped {} long reads'.format(num_LRs))
    print('Done in %.3f s'%(end_time_3-end_time_2))
    print('Constructing matrix and calculating condition number...')
    # unique_dist_df,multi_dist_df = get_long_read_M_dist(long_read_gene_regions_read_count,LR_gene_regions_dict)
    # filtered_unique_dist_df,filtered_multi_dist_df = get_filtered_out_long_read_M_dist(output_path,filtered_gene_regions_read_length,LR_gene_regions_dict)
    
    short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(gene_isoforms_dict,SR_gene_regions_dict,short_read_gene_regions_read_count,SR_read_len,SR_genes_regions_len_dict,num_SRs,region_expression_calculation_method)
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,LR_gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,LR_genes_regions_len_dict,num_LRs,total_long_read_length,region_expression_calculation_method)
    raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict = defaultdict(dict),defaultdict(dict),defaultdict(dict)
    raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            raw_gene_num_exon_dict[chr_name][gene_name] = len(raw_gene_exons_dict[chr_name][gene_name])
            gene_num_exon_dict[chr_name][gene_name] = len(gene_exons_dict[chr_name][gene_name])
            gene_num_isoform_dict[chr_name][gene_name] = len(gene_isoforms_dict[chr_name][gene_name])
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name]:
                raw_isoform_num_exon_dict[isoform_name] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                isoform_length_dict[isoform_name] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                num_isoforms_dict[isoform_name] =  len(raw_isoform_exons_dict[chr_name][gene_name])
    info_dict_list = [raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict,raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict]
    # get_very_short_isoforms(output_path,filtered_gene_regions_read_length,LR_gene_regions_dict,isoform_length_dict)
    # generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list)
    # unique_dist_df.to_csv('{}/lr_M_unique_dist.tsv'.format(output_path),sep='\t',index=False)
    # multi_dist_df.to_csv('{}/lr_M_multi_dist.tsv'.format(output_path),sep='\t',index=False)
    # filtered_unique_dist_df.to_csv('{}/filtered_out_lr_M_unique_dist.tsv'.format(output_path),sep='\t',index=False)
    # filtered_multi_dist_df.to_csv('{}/filtered_out_lr_M_multi_dist.tsv'.format(output_path),sep='\t',index=False)
    # with open('{}/lr.pkl'.format(output_path),'wb') as f:
    #     pickle.dump((gene_isoforms_length_dict,long_read_gene_regions_read_length,long_read_gene_regions_read_count,LR_gene_regions_dict,filtered_gene_regions_read_length),f)

    end_time_4 = time.time()
    print('Done in %.3f s'%(end_time_4-end_time_3))
    gene_isoform_tpm_expression_dict,list_of_all_genes_chrs = quantification(short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,alpha,beta,P)
    
    end_time_5 = time.time()
    # import dill as pickle
    # rep_name = output_path.split('/')[-2]
    # # rep_name = 1
    # pickle.dump((short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_points_dict,SR_gene_regions_dict,LR_gene_regions_dict,gene_isoform_expression_dict,gene_isoform_tpm_expression_dict),open('/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/human_simulation/jobs/hybrid_simulation/validation/quantif_pkl/{}.p'.format(rep_name),'wb'))
    print('Done in %.3f s'%(end_time_5-end_time_4))
    generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict)
