import numpy as np
from numpy import linalg as LA
def construct_index(region_names,isoform_names):
    # indexing the region names and isoform names
    region_names_indics = {x:i for i,x in enumerate(region_names)}
    isoform_names_indics = {x:i for i,x in enumerate(isoform_names)}
    return region_names_indics,isoform_names_indics

def construct_isoform_region_matrix(isoform_region_dict,region_names_indics,isoform_names_indics):
    isoform_region_matrix = np.zeros((len(region_names_indics),len(isoform_names_indics)))
    for region_name in isoform_region_dict:
        for isoform_name in isoform_region_dict[region_name]:
            isoform_region_matrix[region_names_indics[region_name],isoform_names_indics[isoform_name]] = 1
    return isoform_region_matrix

def construct_region_read_count_matrix(region_read_count_dict,region_len_dict,region_names_indics,is_long_read):
    region_read_count_matrix = np.zeros((len(region_names_indics)))
    for region_name in region_read_count_dict:
        if (is_long_read):
            region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name]
        else:
            region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name] / region_len_dict[region_name]
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

# def get_condition_number(isoform_region_matrix):
#     _,singular_values,_ = LA.svd(isoform_region_matrix)
#     svd_val_max = singular_values[0]
#     rank = LA.matrix_rank(isoform_region_matrix)
#     if (rank == min(isoform_region_matrix.shape[0],isoform_region_matrix.shape[1])):
#         # full rank
#         svd_val_min = singular_values[-1]
#         kvalue = (svd_val_max - svd_val_min)/svd_val_max
#         regular_condition_number = generalized_condition_number = svd_val_min/svd_val_max
#     else:
#         # not full rank
#         svd_val_pos_min = singular_values[singular_values > 0].min()
#         kvalue = generalized_condition_number =  (svd_val_max/svd_val_pos_min)
#         regular_condition_number = float("inf")
# #     condition_number = LA.cond(isoform_region_matrix)
#     return kvalue,regular_condition_number,generalized_condition_number
def get_condition_number(isoform_region_matrix):
    rank = LA.matrix_rank(isoform_region_matrix)
    singular_values = LA.svd(isoform_region_matrix.T.dot(isoform_region_matrix)+0.01*np.identity(isoform_region_matrix.shape[1]),compute_uv=False)
    svd_val_max = singular_values[0]
    svd_val_min = singular_values[-1]
    kvalue = svd_val_max/svd_val_min
    singular_values = LA.svd(isoform_region_matrix.T.dot(isoform_region_matrix),compute_uv=False)
    svd_val_max = singular_values[0]
    svd_val_pos_min = singular_values[singular_values > 0].min()
    regular_condition_number = generalized_condition_number = svd_val_max/svd_val_pos_min
    return kvalue,regular_condition_number,generalized_condition_number
def calculate_condition_number(region_isoform_dict,isoform_names):
    region_names = region_isoform_dict.keys()
    (region_names_indics,isoform_names_indics) = construct_index(region_names,isoform_names)
    isoform_region_matrix = construct_isoform_region_matrix(region_isoform_dict,region_names_indics,isoform_names_indics)
    condition_numbers = get_condition_number(isoform_region_matrix)
    matrix_dict = {'isoform_region_matrix':isoform_region_matrix,'condition_number':condition_numbers,
                   'region_names_indics':region_names_indics,'isoform_names_indics':isoform_names_indics}
    return matrix_dict


def generate_feature_matrix(regions_isoform_dict,region_read_count_dict,region_len_dict,isoform_names,is_long_read):
    matrix_dict = calculate_condition_number(regions_isoform_dict,isoform_names)
    matrix_dict['region_read_count_matrix'] = construct_region_read_count_matrix(region_read_count_dict,region_len_dict,matrix_dict['region_names_indics'],is_long_read)
    return matrix_dict

def filter_multi_exons_regions(regions_dict):
    filtered_regions_dict = {}
    for region_name in regions_dict:
        points = [int(p) for p in region_name.replace('P','').replace(':','-').split('-')]
        if len(points) == 2 and points[1] - points[0] == 1:
            filtered_regions_dict[region_name] = regions_dict[region_name]
    return filtered_regions_dict

def calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons):
    gene_matrix_dict = dict()
    for chr_name in gene_isoforms_dict:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_isoforms_dict[chr_name]:
            isoform_names = gene_isoforms_dict[chr_name][gene_name]
            # for short read only allow exon and exon-exon junction
            if (not allow_multi_exons):
                region_isoform_dict = filter_multi_exons_regions(gene_regions_dict[chr_name][gene_name])
            else:
                region_isoform_dict = gene_regions_dict[chr_name][gene_name]
            gene_matrix_dict[chr_name][gene_name] = calculate_condition_number(region_isoform_dict,isoform_names)
    return gene_matrix_dict

def generate_all_feature_matrix(gene_isoforms_dict,gene_regions_dict,gene_regions_read_count,gene_region_len_dict,allow_multi_exons):
    gene_matrix_dict = dict()
    for chr_name in gene_isoforms_dict:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_isoforms_dict[chr_name]:
            isoform_names = gene_isoforms_dict[chr_name][gene_name]
            # for short read only allow exon and exon-exon junction
            if (not allow_multi_exons):
                region_isoform_dict = filter_multi_exons_regions(gene_regions_dict[chr_name][gene_name])
                region_read_count_dict = filter_multi_exons_regions(gene_regions_read_count[chr_name][gene_name])
                region_len_dict = filter_multi_exons_regions(gene_region_len_dict[chr_name][gene_name])
                gene_matrix_dict[chr_name][gene_name] = generate_feature_matrix(region_isoform_dict,region_read_count_dict,region_len_dict,isoform_names,False)
            else:
                region_isoform_dict = gene_regions_dict[chr_name][gene_name]
                region_read_count_dict = gene_regions_read_count[chr_name][gene_name]
                region_len_dict = gene_region_len_dict[chr_name][gene_name]
                gene_matrix_dict[chr_name][gene_name] = generate_feature_matrix(region_isoform_dict,region_read_count_dict,region_len_dict,isoform_names,True)
    return gene_matrix_dict