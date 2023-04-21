import numpy as np
from numpy import linalg as LA
from util import check_region_type,cal_inner_region_len
import config
def construct_index(region_names,isoform_names):
    # indexing the region names and isoform names
    region_names_indics = {x:i for i,x in enumerate(region_names)}
    isoform_names_indics = {x:i for i,x in enumerate(isoform_names)}
    return region_names_indics,isoform_names_indics

def construct_isoform_region_matrix(isoform_region_dict,region_names_indics,isoform_names_indics,normalize_A=True):
    isoform_region_matrix = np.zeros((len(region_names_indics),len(isoform_names_indics)))
    for region_name in isoform_region_dict:
        for isoform_name in isoform_region_dict[region_name]:
            isoform_region_matrix[region_names_indics[region_name],isoform_names_indics[isoform_name]] = 1
    if normalize_A:
        sum_A = isoform_region_matrix.sum(axis=0)
        sum_A[sum_A==0] = 1
        isoform_region_matrix = isoform_region_matrix/sum_A
    return isoform_region_matrix
def calculate_eff_length(region_len_dict,SR_read_len):
    region_eff_length_dict = {}
    for region_name in region_len_dict:
        region_len = region_len_dict[region_name]
        if check_region_type(region_name) == 'one_junction':
            region_eff_length = SR_read_len - config.READ_JUNC_MIN_MAP_LEN
        elif check_region_type(region_name) == 'two_exons':
            region_eff_length = SR_read_len - 1
        elif check_region_type(region_name) == 'one_exon':
            region_eff_length = region_len - SR_read_len + 1 if SR_read_len < region_len else 1
        # else:
        #     try:
        #         inner_region_len = cal_inner_region_len(region_name,region_len_dict)
        #     except:
        #         print(region_len_dict)
        #         print(region_name)
        #     region_eff_length = SR_read_len - config.READ_JUNC_MIN_MAP_LEN
        elif check_region_type(region_name) == 'exons':
            inner_region_len = cal_inner_region_len(region_name,region_len_dict)
            region_eff_length =  SR_read_len - inner_region_len - 1
        elif check_region_type(region_name) == 'junctions':
            inner_region_len = cal_inner_region_len(region_name,region_len_dict)
            region_eff_length = SR_read_len - inner_region_len - config.READ_JUNC_MIN_MAP_LEN
        if region_eff_length <= 0 :
            region_eff_length = 1
        region_eff_length_dict[region_name] = region_eff_length
    return region_eff_length_dict
def construct_region_abundance_matrix_short_read(region_read_count_dict,region_eff_length_dict,region_names_indics,num_SRs):
    region_read_count_matrix = np.zeros((len(region_names_indics)))
    for region_name in region_read_count_dict:
        # if (region_expression_calculation_method == 'coverage'):
        #     region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] * 150 / region_len_dict[region_name]
        # elif (region_expression_calculation_method == 'div_read_length'):
        #     region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / num_SRs
        # elif (region_expression_calculation_method == 'original'):
            # region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / (region_len_dict[region_name] * num_SRs)
        region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / region_eff_length_dict[region_name]
            # region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name]
        # else:
        #     raise Exception('Invalid region expression calculation option!')
    return region_read_count_matrix
        
#     region_tpm_matrix = np.zeros((len(region_names_indics)))
#     region_fpkm_matrix = np.zeros((len(region_names_indics)))
#     for region_name in region_read_count_dict:
#         if (is_long_read):
#             region_tpm_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / region_len_dict[region_name]
#         else:
#             region_tpm_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / region_len_dict[region_name]
# #             region_tpm_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / (region_len_dict[region_name] - 150 + 1)
#         region_fpkm_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / region_len_dict[region_name]
#     #TODO handle 0 read count for whole gene
#     if (not sum(region_tpm_matrix) == 0):
#         region_tpm_matrix = region_tpm_matrix * 1e6 / sum(region_tpm_matrix)
#         region_fpkm_matrix = region_fpkm_matrix * 1e9 / region_fpkm_matrix.shape[0]
#     return region_tpm_matrix,region_fpkm_matrix
def divide_by_zero(a,b):
    if b == 0:
        return float('inf')
    else:
        return a/b
def get_condition_number(isoform_region_matrix):
    # Calculate K value
    multiply_transpose_matrix = isoform_region_matrix.T.dot(isoform_region_matrix)
    singular_values = LA.svd(multiply_transpose_matrix,compute_uv=False)
    rank = LA.matrix_rank(multiply_transpose_matrix)

    svd_val_max = np.sqrt(singular_values[0])
    svd_val_pos_min = np.sqrt(singular_values[singular_values > 0].min())
    svd_val_min = 0
    if (rank == multiply_transpose_matrix.shape[0]):
        svd_val_min = np.sqrt(singular_values[-1])
        # full rank
        kvalue = (svd_val_max - svd_val_min)/svd_val_max
    else:
        # not full rank
        kvalue =  (svd_val_max/svd_val_pos_min)
    
    # Calculate condition number
    regular_condition_number = divide_by_zero(svd_val_max,svd_val_min)

    # Calculate generalized condition number
    generalized_condition_number = svd_val_max/svd_val_pos_min
    if (rank == multiply_transpose_matrix.shape[0]):
        assert regular_condition_number == generalized_condition_number

    singular_value_product = svd_val_max * svd_val_pos_min
    return kvalue,regular_condition_number,generalized_condition_number,singular_value_product
def calculate_condition_number(region_isoform_dict,isoform_names,normalize_A):
    region_names = region_isoform_dict.keys()
    (region_names_indics,isoform_names_indics) = construct_index(region_names,isoform_names)
    isoform_region_matrix = construct_isoform_region_matrix(region_isoform_dict,region_names_indics,isoform_names_indics,normalize_A)
    try:
        condition_numbers = get_condition_number(isoform_region_matrix)
    except:
        # print(isoform_region_matrix)
        # print(region_names_indics)
        # print(isoform_names_indics)
        condition_numbers = (np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan'))
    matrix_dict = {'isoform_region_matrix':isoform_region_matrix,'condition_number':condition_numbers,
                   'region_names_indics':region_names_indics,'isoform_names_indics':isoform_names_indics}
    return matrix_dict
def cal_weight_multi_exon_region(short_read_gene_matrix_dict,region_len_dict,SR_read_len):
    isoform_region_matrix = short_read_gene_matrix_dict['isoform_region_matrix'].copy()
    for isoform_name,isoform_index in short_read_gene_matrix_dict['isoform_names_indics'].items():
        for region,region_index in short_read_gene_matrix_dict['region_names_indics'].items():
            if isoform_region_matrix[region_index][isoform_index] != 0:
                if check_region_type(region) in ['junctions','exons']:
                    inner_region_len = cal_inner_region_len(region,region_len_dict)
                    weight = (SR_read_len - inner_region_len - 1)/SR_read_len
                    if weight <= 0 or weight > 1:
                        weight = 1
                    isoform_region_matrix[region_index][isoform_index] = weight
    return isoform_region_matrix

# def filter_regions(regions_dict,long_read = False):
#     filtered_regions_dict = {}
#     for region_name in regions_dict:
#         try:
#             points = [int(p) for p in region_name.replace('P','').replace(':','-').split('-')]
#         except:
#             print(region_name)
#         if (not long_read):
#             if check_region_type(region_name) in ['two_exons','one_junction','one_exon']:
#                 filtered_regions_dict[region_name] = regions_dict[region_name]
#         else:
#             if check_region_type(region_name) in ['two_exons','one_junction','one_exon','others']:
#                 filtered_regions_dict[region_name] = regions_dict[region_name]

#     return filtered_regions_dict

def calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons):
    gene_matrix_dict = dict()
    for chr_name in gene_isoforms_dict:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_isoforms_dict[chr_name]:
            isoform_names = gene_isoforms_dict[chr_name][gene_name]
            # for short read only allow exon and exon-exon junction
            region_isoform_dict = gene_regions_dict[chr_name][gene_name]
            # if (not allow_multi_exons):
            #     region_isoform_dict = filter_regions(gene_regions_dict[chr_name][gene_name],long_read=False)
            # else:
            #     region_isoform_dict = filter_regions(gene_regions_dict[chr_name][gene_name],long_read=True)
            gene_matrix_dict[chr_name][gene_name] = calculate_condition_number(region_isoform_dict,isoform_names,config.normalize_sr_A)
    return gene_matrix_dict
def generate_all_feature_matrix_short_read(gene_isoforms_dict,gene_regions_dict,gene_regions_read_count,SR_read_len,gene_region_len_dict,num_SRs,normalize_A=True):
    gene_matrix_dict = dict()
    for chr_name in gene_isoforms_dict:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_isoforms_dict[chr_name]:
            isoform_names = gene_isoforms_dict[chr_name][gene_name]
            # for short read only allow exon and exon-exon junction
            # region_isoform_dict = filter_regions(gene_regions_dict[chr_name][gene_name],long_read=False)
            # region_read_count_dict = filter_regions(gene_regions_read_count[chr_name][gene_name],long_read=False)
            # region_len_dict = filter_regions(gene_region_len_dict[chr_name][gene_name],long_read=False)
            region_isoform_dict = {}
            for region in gene_regions_read_count[chr_name][gene_name]:
                region_isoform_dict[region] = gene_regions_dict[chr_name][gene_name][region]
            # region_isoform_dict = gene_regions_dict[chr_name][gene_name]
            region_read_count_dict = gene_regions_read_count[chr_name][gene_name]
            region_len_dict = gene_region_len_dict[chr_name][gene_name]
    
            matrix_dict = calculate_condition_number(region_isoform_dict,isoform_names,config.normalize_sr_A)
            matrix_dict['region_eff_length_dict'] = calculate_eff_length(region_len_dict,SR_read_len)
            if config.multi_exon_region_weight == 'minus_inner_region':
                matrix_dict['isoform_region_matrix'] = cal_weight_multi_exon_region(matrix_dict,region_len_dict,SR_read_len)
                if config.normalize_sr_A:
                    sum_A = matrix_dict['isoform_region_matrix'].sum(axis=0)
                    sum_A[sum_A==0] = 1
                    matrix_dict['isoform_region_matrix'] = matrix_dict['isoform_region_matrix']/sum_A
                matrix_dict['condition_number'] = get_condition_number(matrix_dict['isoform_region_matrix'])
            matrix_dict['region_abund_matrix'] = construct_region_abundance_matrix_short_read(region_read_count_dict,matrix_dict['region_eff_length_dict'],matrix_dict['region_names_indics'],num_SRs)
            num_SRs_mapped_gene = 0
            for region in region_read_count_dict:
                num_SRs_mapped_gene += region_read_count_dict[region]
            matrix_dict['num_SRs_mapped_gene'] = num_SRs_mapped_gene
            gene_matrix_dict[chr_name][gene_name] = matrix_dict

    return gene_matrix_dict
def is_multi_isoform_region(matrix_dict,region):
    index = matrix_dict['region_names_indics'][region]
    A = matrix_dict['isoform_region_matrix'].copy()
    A[A != 0] = 1
    sum_A = A.sum(axis=1)
    return sum_A[index] > 1