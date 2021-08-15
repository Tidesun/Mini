from collections import defaultdict
import numpy as np
from numpy import linalg as LA
from qpsolvers import solve_qp
from predict_params import load_model,predict_params
def normalize_expression(gene_isoform_expression_dict):
    gene_isoform_tpm_expression_dict = defaultdict(lambda: defaultdict(dict))
    SR_isoform_expression_sum = 0
    LR_isoform_expression_sum = 0
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_expected_counts'] = gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts']
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_expected_counts'] = gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression']
            SR_isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'].sum()
            LR_isoform_expression_sum += gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'].sum()

    if SR_isoform_expression_sum == 0:
        SR_isoform_expression_sum = 1
    if LR_isoform_expression_sum == 0:
        LR_isoform_expression_sum = 1 
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_tpm'] = gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] * 1e6 / SR_isoform_expression_sum
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_tpm'] = gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'] * 1e6 / LR_isoform_expression_sum
            alpha = gene_isoform_expression_dict[chr_name][gene_name]['alpha']
            gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'] = (1 - alpha) * gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_tpm'] + alpha * gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_tpm']
            if gene_name == 'ENSDARG00000114503':
                print(gene_isoform_tpm_expression_dict[chr_name][gene_name])
    return gene_isoform_tpm_expression_dict 
def estimate_isoform_expression_grid_search_iteration(args,params):
    (SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P) = args
    if SR_region_read_count_matrix.sum() != 0:
        SR_region_read_count_matrix = SR_region_read_count_matrix / SR_region_read_count_matrix.sum()
    if LR_region_read_count_matrix.sum() != 0:
        LR_region_read_count_matrix = LR_region_read_count_matrix / LR_region_read_count_matrix.sum()
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
    # isoform_expression = solve_qp(Q, c,G,h, lb = lb)
    isoform_expression = solve_qp(Q, c,lb = lb)
    # if ((isoform_expression+1e-6<0).any()):
    #     raise ValueError('Obtain negative value for isoform expression')
    isoform_expression[isoform_expression<0] = 0
    target = (1.0 - alpha) * LA.norm(SR_region_read_count_matrix - np.matmul(SR_isoform_region_matrix,isoform_expression)) + alpha * LA.norm(LR_region_read_count_matrix - np.matmul(LR_isoform_region_matrix,isoform_expression))
    # 
    perfect_isoform_expression = np.matmul(LA.inv((1.0 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix)+ beta * np.identity(num_isoforms)), (1 - alpha) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix) + alpha * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix) )
    return isoform_expression,perfect_isoform_expression,target


def estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,SR_gene_counts,alpha,beta,P,model):
    num_isoforms = SR_isoform_region_matrix.shape[1]
    if ((SR_region_read_count_matrix<=0).all() and (LR_region_read_count_matrix<=0).all()):
        return np.zeros(num_isoforms),np.zeros(num_isoforms),np.zeros(num_isoforms),0.5
    if (alpha == 'adaptive' or beta == 'adaptive'):
        pred_alpha,pred_beta = predict_params(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,model)
    if (alpha == 'adaptive'):
        # alpha_selections = [i/20 for i in range(20)]
        selected_alpha = pred_alpha
    else:
        selected_alpha = alpha
    if (beta == 'adaptive'):
        # beta_selections = [10**(-i) for i in range(1,10)]
        selected_beta = pred_beta
    else:
        selected_beta = beta
    args = (SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,P)
    params = {'alpha':selected_alpha,'beta':selected_beta}
    # for params in params_grid:
    isoform_expression,_, _ = estimate_isoform_expression_grid_search_iteration(args,params)
    if isoform_expression.sum() != 0:
        isoform_expression = isoform_expression/isoform_expression.sum()
    SR_expression = isoform_lengths * isoform_expression
    if SR_expression.sum() != 0:
        SR_expression = SR_expression/SR_expression.sum()
    SR_expected_counts = SR_gene_counts * SR_expression
    SR_isoform_expression = SR_expected_counts / isoform_lengths

    LR_isoform_expression = LR_region_read_count_matrix.sum() * isoform_expression
    return SR_isoform_expression,SR_expected_counts,LR_isoform_expression,params['alpha']
def estimate_isoform_expression_single_gene(args):
    (short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,alpha,beta,P,model) = args
    SR_isoform_region_matrix = short_read_gene_matrix_dict['isoform_region_matrix']
    SR_region_read_count_matrix = short_read_gene_matrix_dict['region_abund_matrix']
    SR_gene_counts = short_read_gene_matrix_dict['num_SRs_mapped_gene']
    LR_isoform_region_matrix = long_read_gene_matrix_dict['isoform_region_matrix']
    LR_region_read_count_matrix = long_read_gene_matrix_dict['region_abund_matrix']
    isoform_lengths = np.zeros((len(gene_isoforms_length_dict)))
    isoform_names_indics = short_read_gene_matrix_dict['isoform_names_indics']
    for isoform_name in isoform_names_indics:
        isoform_lengths[isoform_names_indics[isoform_name]] = gene_isoforms_length_dict[isoform_name]
    return estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,SR_gene_counts,alpha,beta,P,model)
def quantification(short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,alpha,beta,P):
    print('Calculating the isoform expression...')
    gene_isoform_expression_dict = defaultdict(lambda:defaultdict(dict))
    if (alpha == 'adaptive' or beta == 'adaptive'):
        model = load_model('model_10.pt')
    else:
        model = None
    list_of_all_genes_chrs = []
    for chr_name in long_read_gene_matrix_dict:
        if chr_name in short_read_gene_matrix_dict:
            for gene_name in long_read_gene_matrix_dict[chr_name]:
                if gene_name in short_read_gene_matrix_dict[chr_name]:
                    list_of_all_genes_chrs.append((gene_name,chr_name))
    list_of_args = [(short_read_gene_matrix_dict[chr_name][gene_name],long_read_gene_matrix_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],alpha,beta,P,model) for gene_name,chr_name in list_of_all_genes_chrs]
    for (gene_name,chr_name), args in zip(list_of_all_genes_chrs, list_of_args):
        try:
            result = estimate_isoform_expression_single_gene(args)
            gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'],gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts'],gene_isoform_expression_dict[chr_name][gene_name]['LR_isoform_expression'],gene_isoform_expression_dict[chr_name][gene_name]['alpha'] = result
        except Exception as e:
            print(e)
    gene_isoform_tpm_expression_dict = normalize_expression(gene_isoform_expression_dict)
    return gene_isoform_tpm_expression_dict,list_of_all_genes_chrs