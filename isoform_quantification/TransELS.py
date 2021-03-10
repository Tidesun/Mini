from collections import defaultdict
import concurrent.futures
import datetime
from itertools import repeat

import numpy as np
from numpy import linalg as LA
from qpsolvers import solve_qp

from construct_feature_matrix import generate_all_feature_matrix_short_read,generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment
from parse_alignment_main import parse_alignment
from generate_output import generate_TransELS_output


def normalize_expression(gene_isoform_expression_dict,total_num_reads):
    gene_isoform_tpm_expression_dict = defaultdict(dict)
    isoform_expression_sum = 0
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name].sum()
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            gene_isoform_tpm_expression_dict[chr_name][gene_name] = gene_isoform_expression_dict[chr_name][gene_name] * 1e6 / isoform_expression_sum
    return gene_isoform_tpm_expression_dict
# def estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,alpha,beta):
#     num_isoforms = SR_isoform_region_matrix.shape[1]
#     Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta * np.identity(num_isoforms)
#     c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T)
#     lb = np.zeros(num_isoforms)
#     isoform_expression = solve_qp(Q, c, lb = lb)
#     if ((isoform_expression+1e-5<0).any()):
#         raise ValueError('Obtain negative value for isoform expression')
#     isoform_expression[isoform_expression<0] = 0
#     def div0( a, b ):
#         with np.errstate(divide='ignore', invalid='ignore'):
#             c = np.true_divide(a, b)
#             c[ ~ np.isfinite( c )] = 0
#         return c
#     isoform_expression = div0(isoform_expression,isoform_lengths)
#     return isoform_expression   

def estimate_isoform_expression_grid_search_iteration(args,params):
    (SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P) = args
    alpha,beta = params['alpha'], params['beta']
    num_isoforms = SR_isoform_region_matrix.shape[1]
    Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta * np.identity(num_isoforms)
    c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T)
    lb = np.zeros(num_isoforms)
    G = - np.concatenate((SR_isoform_region_matrix[SR_region_read_count_matrix>0,:], LR_isoform_region_matrix[LR_region_read_count_matrix>0,:]), axis=0)
    h = - np.ones(G.shape[0])/(1/P)
    isoform_expression = solve_qp(Q, c,G,h, lb = lb)
    # if ((isoform_expression+1e-6<0).any()):
    #     print(alpha)
    #     print(beta)
    #     print(isoform_expression)
    #     raise ValueError('Obtain negative value for isoform expression')
    isoform_expression[isoform_expression<0] = 0
    target = (1.0 - alpha) * LA.norm(SR_region_read_count_matrix - np.matmul(SR_isoform_region_matrix,isoform_expression)) + alpha * LA.norm(LR_region_read_count_matrix - np.matmul(LR_isoform_region_matrix,isoform_expression))
    return isoform_expression,target


def estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,alpha,beta,P):
    num_isoforms = SR_isoform_region_matrix.shape[1]
    if ((SR_region_read_count_matrix<=0).all() & (LR_region_read_count_matrix<=0).all()):
        return np.zeros(num_isoforms)
    min_target = float('inf')
    best_params = None
    best_isoform_expression = np.zeros(num_isoforms)
    alpha_selections = [i/20 for i in range(20)]
    beta_selections = [10**(-i) for i in range(1,10)]
    params_grid =[{'alpha':alpha,'beta':beta} for alpha in alpha_selections for beta in beta_selections]
    args = (SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P)
    for params in params_grid:
        isoform_expression,target = estimate_isoform_expression_grid_search_iteration(args,params)
        if (target < min_target):
                min_target = target
                best_isoform_expression = isoform_expression
                best_params = params
    # with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #     for params, (isoform_expression,target) in zip(params_grid, executor.map(estimate_isoform_expression_grid_search_iteration,repeat(args),params_grid,chunksize=100)):
    #         if (target < min_target):
    #             min_target = target
    #             best_isoform_expression = isoform_expression
    #             best_params = params
        # future_grid_search = {executor.submit(estimate_isoform_expression_grid_search_iteration, args,params): params for params in params_grid}
        # for future in concurrent.futures.as_completed(future_grid_search):
        #     params = future_grid_search[future]
        #     try:
        #         isoform_expression,target = future.result()
        #         if (target < min_target):
        #             min_target = target
        #             best_isoform_expression = isoform_expression
        #             best_params = params
        #     except Exception as exc:
        #         print('%r generated an exception: %s' % (params, exc))
    return best_isoform_expression 

def estimate_isoform_expression_single_chr(per_chr_short_read_gene_matrix_dict,per_chr_long_read_gene_matrix_dict,per_chr_gene_isoforms_length_dict,alpha,beta,P):
    gene_list = set(per_chr_short_read_gene_matrix_dict.keys()).intersection(set(per_chr_long_read_gene_matrix_dict.keys()))
    per_chr_gene_isoform_expression_dict = {}
    for gene_name in gene_list:
        isoform_lengths = np.zeros((len(per_chr_gene_isoforms_length_dict[gene_name])))
        isoform_names_indics = per_chr_short_read_gene_matrix_dict[gene_name]['isoform_names_indics']
        for isoform_name in isoform_names_indics:
            isoform_lengths[isoform_names_indics[isoform_name]] = per_chr_gene_isoforms_length_dict[gene_name][isoform_name]
            try:
                per_chr_gene_isoform_expression_dict[gene_name] = estimate_isoform_expression(
                                        per_chr_short_read_gene_matrix_dict[gene_name]['isoform_region_matrix'],
                                        per_chr_short_read_gene_matrix_dict[gene_name]['region_read_count_matrix'],
                                        per_chr_long_read_gene_matrix_dict[gene_name]['isoform_region_matrix'],
                                        per_chr_long_read_gene_matrix_dict[gene_name]['region_read_count_matrix'],isoform_lengths,alpha,beta,P)
            except ValueError as e:
                print("Error encountered for gene %s:"%(gene_name) + e)
    return per_chr_gene_isoform_expression_dict
def TransELS(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,region_expression_calculation_method,alpha,beta,P,threads=1,READ_LEN=79,READ_JUNC_MIN_MAP_LEN=10):
    start_time = datetime.datetime.now()
    print('Start parsing annoation...')
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    gene_regions_points_list,gene_range = process_annotation_for_alignment(gene_points_dict,gene_regions_dict)
    end_time_1 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_1-start_time).total_seconds()))
    print('Mapping short read to regions...')  
    short_read_gene_regions_read_count,num_SRs = parse_alignment(short_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict,long_read = False)
    end_time_2 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_2-end_time_1).total_seconds()))
    print('Mapping long read to regions...')  
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict, long_read = True)
    end_time_3 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_3-end_time_2).total_seconds()))
    print('Constructing matrix and calculating condition number...')
    short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(gene_isoforms_dict,gene_regions_dict,short_read_gene_regions_read_count,genes_regions_len_dict,num_SRs,region_expression_calculation_method)
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,genes_regions_len_dict,num_LRs,total_long_read_length,region_expression_calculation_method)
    chr_list = set(short_read_gene_matrix_dict.keys()).intersection(set(long_read_gene_matrix_dict.keys()))
    gene_isoform_expression_dict = defaultdict(dict)
    end_time_4 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_4-end_time_3).total_seconds()))
    print('Calculating the isoform expression...')
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        future_estimate = {executor.submit(estimate_isoform_expression_single_chr,short_read_gene_matrix_dict[chr_name],long_read_gene_matrix_dict[chr_name],gene_isoforms_length_dict[chr_name],alpha,beta,P): chr_name for chr_name in chr_list}
        for future in concurrent.futures.as_completed(future_estimate):
            chr_name = future_estimate[future]
            try:
                gene_isoform_expression_dict[chr_name] = future.result()
            except Exception as exc:
                print('%r generated an exception: chr %s' % (chr_name, exc))
    # for chr_name in chr_list:
    #     gene_isoform_expression_dict[chr_name] = estimate_isoform_expression_single_chr(short_read_gene_matrix_dict[chr_name],long_read_gene_matrix_dict[chr_name],gene_isoforms_length_dict[chr_name],alpha,beta,P)
    gene_isoform_tpm_expression_dict = normalize_expression(gene_isoform_expression_dict,num_SRs+num_LRs)
    end_time_5 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_5-end_time_4).total_seconds()))
    generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict)
    return short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict