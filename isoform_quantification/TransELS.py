from collections import defaultdict
import concurrent.futures
import datetime
import time
from itertools import repeat
import pickle
import numpy as np
from numpy import linalg as LA
from qpsolvers import solve_qp

from construct_feature_matrix import generate_all_feature_matrix_short_read,generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment
from parse_alignment_main import parse_alignment
from generate_output import generate_TransELS_output
from get_long_read_gene_distribution import get_long_read_gene_distribution

def adjust_isoform_expression_by_gene_expression(gene_isoform_expression_dict,gene_isoforms_length_dict,short_read_gene_matrix_dict,SR_read_len,long_read_gene_matrix_dict):
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            if gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'].sum() == 0:
                continue
            # isoform_eff_length_arr = np.array([gene_isoforms_length_dict[chr_name][gene_name][isoform] - SR_read_len + 1 for isoform in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']])
            # isoform_eff_length_arr[isoform_eff_length_arr < 1] = 1
            # adjusted_isoform_expression = gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] * isoform_eff_length_arr
            # gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] = adjusted_isoform_expression/adjusted_isoform_expression.sum() * short_read_gene_matrix_dict[chr_name][gene_name]['num_SRs_mapped_gene'] / isoform_eff_length_arr
            gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] = gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] / gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'].sum() * long_read_gene_matrix_dict[chr_name][gene_name]['num_LRs_mapped_gene']
    return gene_isoform_expression_dict
def normalize_expression(gene_isoform_expression_dict,gene_isoforms_length_dict,short_read_gene_matrix_dict,total_num_reads):
    gene_isoform_tpm_expression_dict = defaultdict(lambda: defaultdict(dict))
    isoform_expression_sum = 0
    SR_L2_sum,SR_perfect_sum,LR_L2_sum,LR_perfect_sum = 0,0,0,0
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'].sum()
            
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            # gene_isoform_length_arr = np.array([gene_isoforms_length_dict[chr_name][gene_name][isoform] for isoform in short_read_gene_matrix_dict[chr_name][gene_name]['']])
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'] = gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] * 1e6 / isoform_expression_sum
            # gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm_by_gene'] = gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] * 1e6 / gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'].sum() 
            # gene_isoform_tpm_expression_dict[chr_name][gene_name]['perfect_tpm_by_gene'] = gene_isoform_expression_dict[chr_name][gene_name]['perfect_isoform_expression'] * 1e6 / gene_isoform_expression_dict[chr_name][gene_name]['perfect_isoform_expression'].sum()
            # gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm_by_gene_sr'] = num_mapped_SRs[chr_name][gene_name] * gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'] / (gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression']  * gene_isoform_length_arr)
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
    #l2 norm
    Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta * np.identity(num_isoforms)
    c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T)
    lb = np.zeros(num_isoforms)
    G = - np.concatenate((SR_isoform_region_matrix[SR_region_read_count_matrix>0,:], LR_isoform_region_matrix[LR_region_read_count_matrix>0,:]), axis=0)
    h = - np.ones(G.shape[0])/(1/P)
    # l1 norm
    # Q = 2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix)
    # c = -2 * (1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T) + beta* np.ones(num_isoforms)
    # lb = np.zeros(num_isoforms)
    # G = - np.concatenate((SR_isoform_region_matrix[SR_region_read_count_matrix>0,:], LR_isoform_region_matrix[LR_region_read_count_matrix>0,:]), axis=0)
    # h = - np.ones(G.shape[0])/(1/P)
    isoform_expression = solve_qp(Q, c,G,h, lb = lb)
    # isoform_expression = solve_qp(Q, c,lb = lb)
    # if ((isoform_expression+1e-6<0).any()):
    #     raise ValueError('Obtain negative value for isoform expression')
    isoform_expression[isoform_expression<0] = 0
    target = (1.0 - alpha) * LA.norm(SR_region_read_count_matrix - np.matmul(SR_isoform_region_matrix,isoform_expression)) + alpha * LA.norm(LR_region_read_count_matrix - np.matmul(LR_isoform_region_matrix,isoform_expression))
    # 
    perfect_isoform_expression = np.matmul(LA.inv((1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix)+ beta * np.identity(num_isoforms)), (1 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix) )
    return isoform_expression,perfect_isoform_expression,target


def estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,alpha,beta,P):
    num_isoforms = SR_isoform_region_matrix.shape[1]
    if ((SR_region_read_count_matrix<=0).all() and (LR_region_read_count_matrix<=0).all()):
        return np.zeros(num_isoforms),np.zeros(num_isoforms)
    min_target = float('inf')
    best_params = None
    best_isoform_expression = np.zeros(num_isoforms)
    best_perfect_isoform_expression = np.zeros(num_isoforms)
    if (alpha == 'adaptive'):
        alpha_selections = [i/20 for i in range(20)]
    else:
        alpha_selections = [alpha]
    if (beta == 'adaptive'):
        beta_selections = [10**(-i) for i in range(1,10)]
    else:
        beta_selections = [beta]
    params_grid =[{'alpha':alpha,'beta':beta} for alpha in alpha_selections for beta in beta_selections]
    args = (SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P)
    for params in params_grid:
        isoform_expression,perfect_isoform_expression, target = estimate_isoform_expression_grid_search_iteration(args,params)
        if (target < min_target):
                min_target = target
                best_isoform_expression = isoform_expression
                best_perfect_isoform_expression = perfect_isoform_expression
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
    return best_isoform_expression,best_perfect_isoform_expression
def estimate_isoform_expression_single_gene(args):
    (short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,alpha,beta,P) = args
    SR_isoform_region_matrix = short_read_gene_matrix_dict['isoform_region_matrix']
    SR_region_read_count_matrix = short_read_gene_matrix_dict['region_abund_matrix']
    LR_isoform_region_matrix = long_read_gene_matrix_dict['isoform_region_matrix']
    LR_region_read_count_matrix = long_read_gene_matrix_dict['region_abund_matrix']
    isoform_lengths = np.zeros((len(gene_isoforms_length_dict)))
    isoform_names_indics = short_read_gene_matrix_dict['isoform_names_indics']
    for isoform_name in isoform_names_indics:
        isoform_lengths[isoform_names_indics[isoform_name]] = gene_isoforms_length_dict[isoform_name]
    return estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,alpha,beta,P)
# def estimate_isoform_expression_single_chr(per_chr_short_read_gene_matrix_dict,per_chr_long_read_gene_matrix_dict,per_chr_gene_isoforms_length_dict,alpha,beta,P):
#     gene_list = set(per_chr_short_read_gene_matrix_dict.keys()).intersection(set(per_chr_long_read_gene_matrix_dict.keys()))
#     per_chr_gene_isoform_expression_dict = {}
#     for gene_name in gene_list:
#         isoform_lengths = np.zeros((len(per_chr_gene_isoforms_length_dict[gene_name])))
#         isoform_names_indics = per_chr_short_read_gene_matrix_dict[gene_name]['isoform_names_indics']
#         for isoform_name in isoform_names_indics:
#             isoform_lengths[isoform_names_indics[isoform_name]] = per_chr_gene_isoforms_length_dict[gene_name][isoform_name]
#         try:
#             per_chr_gene_isoform_expression_dict[gene_name] = {}
#             per_chr_gene_isoform_expression_dict[gene_name]['isoform_expression'],per_chr_gene_isoform_expression_dict[gene_name]['perfect_isoform_expression'] = estimate_isoform_expression(
#                                     per_chr_short_read_gene_matrix_dict[gene_name]['isoform_region_matrix'],
#                                     per_chr_short_read_gene_matrix_dict[gene_name]['region_abund_matrix'],
#                                     per_chr_long_read_gene_matrix_dict[gene_name]['isoform_region_matrix'],
#                                     per_chr_long_read_gene_matrix_dict[gene_name]['region_abund_matrix'],isoform_lengths,alpha,beta,P)
#         except ValueError as e:
#             print("Error encountered for gene {} :{}".format(gene_name,e))
#     return per_chr_gene_isoform_expression_dict
def TransELS(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,region_expression_calculation_method,alpha,beta,P,threads=1,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=10):
    start_time = time.time()
    print('Preprocessing...')
    LR_gene_read_min_len_dict = None
    # LR_gene_read_min_len_dict = get_long_read_gene_distribution(ref_file_path,long_read_alignment_file_path)
    # print(LR_gene_read_min_len_dict)
    print('Start parsing annoation...')
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,_ = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,LR_gene_read_min_len_dict)
    gene_regions_points_list,gene_range,gene_interval_tree_dict = process_annotation_for_alignment(gene_exons_dict,gene_points_dict)
    end_time_1 = time.time()
    print('Done in %.3f s'%(end_time_1-start_time))
    print('Mapping short read to regions...')  
    short_read_gene_regions_read_count,SR_read_len,num_SRs = parse_alignment(short_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,False,threads)
    end_time_2 = time.time()
    print('Mapped {} short reads'.format(num_SRs))
    print('Done in %.3f s'%(end_time_2-end_time_1))
    print('Mapping long read to regions...')
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict, True,threads)
    end_time_3 = time.time()
    print('Mapped {} long reads'.format(num_LRs))
    print('Done in %.3f s'%(end_time_3-end_time_2))
    print('Constructing matrix and calculating condition number...')
    short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(gene_isoforms_dict,SR_gene_regions_dict,short_read_gene_regions_read_count,SR_read_len,SR_genes_regions_len_dict,num_SRs,region_expression_calculation_method)
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,LR_gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,LR_genes_regions_len_dict,num_LRs,total_long_read_length,region_expression_calculation_method)
    gene_isoform_expression_dict = defaultdict(lambda:defaultdict(dict))
    end_time_4 = time.time()
    print('Done in %.3f s'%(end_time_4-end_time_3))
    print('Calculating the isoform expression...')
    list_of_all_genes_chrs = []
    pickle.dump((short_read_gene_matrix_dict,long_read_gene_matrix_dict),open('/fs/ess/scratch/PCON0009/haoran/TransELS/a.pkl','wb'))
    for chr_name in long_read_gene_matrix_dict:
        if chr_name in short_read_gene_matrix_dict:
            for gene_name in long_read_gene_matrix_dict[chr_name]:
                if gene_name in short_read_gene_matrix_dict[chr_name]:
                    list_of_all_genes_chrs.append((gene_name,chr_name))
    list_of_args = [(short_read_gene_matrix_dict[chr_name][gene_name],long_read_gene_matrix_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],alpha,beta,P) for gene_name,chr_name in list_of_all_genes_chrs]
    # if threads == 1:
    if True:
        for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, [estimate_isoform_expression_single_gene(args) for args in list_of_args]):
            try:
                gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'],gene_isoform_expression_dict[chr_name][gene_name]['perfect_isoform_expression'] = result
            except Exception as e:
                print(e)
                raise e
    else:
        chunksize, extra = divmod(len(list_of_all_genes_chrs), threads)
        if extra:
            chunksize += 1
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, executor.map(estimate_isoform_expression_single_gene,list_of_args,chunksize=chunksize)):
                try:
                    gene_isoform_expression_dict[chr_name][gene_name]['isoform_expression'],gene_isoform_expression_dict[chr_name][gene_name]['perfect_isoform_expression'] = result
                except Exception as e:
                    print(e)
                    raise e

    # for chr_name in chr_list:
    #     gene_isoform_expression_dict[chr_name] = estimate_isoform_expression_single_chr(short_read_gene_matrix_dict[chr_name],long_read_gene_matrix_dict[chr_name],gene_isoforms_length_dict[chr_name],alpha,beta,P)
    # gene_isoform_expression_dict = adjust_isoform_expression_by_gene_expression(gene_isoform_expression_dict,gene_isoforms_length_dict,short_read_gene_matrix_dict,SR_read_len,long_read_gene_matrix_dict)   
    gene_isoform_tpm_expression_dict = normalize_expression(gene_isoform_expression_dict,gene_isoforms_length_dict,short_read_gene_matrix_dict,num_SRs+num_LRs)
    end_time_5 = time.time()
    # import dill as pickle
    # rep_name = output_path.split('/')[-2]
    # # rep_name = 1
    # pickle.dump((short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_points_dict,SR_gene_regions_dict,LR_gene_regions_dict,gene_isoform_expression_dict,gene_isoform_tpm_expression_dict),open('/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/human_simulation/jobs/hybrid_simulation/validation/quantif_pkl/{}.p'.format(rep_name),'wb'))
    print('Done in %.3f s'%(end_time_5-end_time_4))
    generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict)