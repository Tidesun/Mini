from pathlib import Path
from kde_weight_calculation import create_new_matrix_dict,cal_isoform_region_weight
import numpy as np
import pickle
from construct_feature_matrix import calculate_condition_number,is_multi_isoform_region,get_condition_number
import config
def construct_region_abundance_matrix_long_read(region_read_length,region_read_count_dict,region_len_dict,region_names_indics,num_LRs,total_long_read_lengths):
    region_read_count_matrix = np.zeros((len(region_names_indics)))
    for region_name in region_read_count_dict:
        # if (region_expression_calculation_method == 'coverage'):
            # region_read_count_matrix[region_names_indics[region_name]] = sum(region_read_length[region_name]) / region_len_dict[region_name]
        # elif (region_expression_calculation_method == 'div_read_length'):
        #     region_read_count_matrix[region_names_indics[region_name]] = sum(region_read_length[region_name]) / total_long_read_lengths
        # elif (region_expression_calculation_method == 'original'):
        region_read_count_matrix[region_names_indics[region_name]] = region_read_count_dict[region_name]
        # else:
        #     raise Exception('Invalid region expression calculation option!')
    return region_read_count_matrix

def generate_all_feature_matrix_long_read(gene_isoforms_dict,gene_regions_dict,gene_regions_read_count,gene_regions_read_length,gene_region_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,num_LRs,total_long_read_length,READ_JUNC_MIN_MAP_LEN,output_dir,threads,normalize_A=True):
    # with open('/users/PCON0009/haoranli/_projects/complex_profile_quant/weight_dict.pkl','rb') as f:
    #     long_reads_isoform_region_weight_matrix_dict = pickle.load(f)
    with open('/fs/ess/scratch/PCON0009/haoran/weight_calculation/weight_dict.pkl','rb') as f:
        long_reads_isoform_region_weight_matrix_dict = pickle.load(f)
            
    # long_reads_isoform_region_weight_matrix_dict = cal_isoform_region_weight(gene_regions_dict,gene_region_len_dict,gene_isoforms_length_dict,READ_JUNC_MIN_MAP_LEN,output_dir,threads)
    # long_reads_isoform_region_weight_matrix_dict = {}
    gene_matrix_dict = dict()
    for chr_name in gene_regions_read_count:
        gene_matrix_dict[chr_name] = dict()
        for gene_name in gene_regions_read_count[chr_name]:
            isoform_names = gene_isoforms_dict[chr_name][gene_name]

            # region_isoform_dict = filter_regions(gene_regions_dict[chr_name][gene_name],long_read=True)
            # region_read_count_dict = filter_regions(gene_regions_read_count[chr_name][gene_name],long_read=True)
            # region_len_dict = filter_regions(gene_region_len_dict[chr_name][gene_name],long_read=True)
            # region_read_length = filter_regions(gene_regions_read_length[chr_name][gene_name],long_read=True)
            region_isoform_dict = {}
            for region in gene_regions_read_count[chr_name][gene_name]:
                region_isoform_dict[region] = gene_regions_dict[chr_name][gene_name][region]


            # region_isoform_dict = gene_regions_dict[chr_name][gene_name]
            region_read_count_dict = gene_regions_read_count[chr_name][gene_name]
            region_len_dict = gene_region_len_dict[chr_name][gene_name]
            region_read_length = gene_regions_read_length[chr_name][gene_name]

            matrix_dict = calculate_condition_number(region_isoform_dict,isoform_names,False)
            matrix_dict['region_abund_matrix'] = construct_region_abundance_matrix_long_read(region_read_length,region_read_count_dict,region_len_dict,matrix_dict['region_names_indics'],num_LRs,total_long_read_length)
            num_LRs_mapped_gene = 0
            for region in region_read_count_dict:
                num_LRs_mapped_gene += region_read_count_dict[region]
            matrix_dict['num_LRs_mapped_gene'] = num_LRs_mapped_gene
            matrix_dict['num_exons'] = {}
            for isoform in isoform_names:
                matrix_dict['num_exons'][isoform] = len(raw_isoform_exons_dict[chr_name][gene_name][isoform]['start_pos'])
            matrix_dict['multi_isoforms_count'],matrix_dict['unique_isoforms_count'] = 0, 0
            for region in matrix_dict['region_names_indics']:
                index = matrix_dict['region_names_indics'][region]
                count = matrix_dict['region_abund_matrix'][index]
                if is_multi_isoform_region(matrix_dict,region):
                    matrix_dict['multi_isoforms_count'] += count
                else:
                    matrix_dict['unique_isoforms_count'] += count
            if config.use_weight_matrix:
                if chr_name in long_reads_isoform_region_weight_matrix_dict:
                    if gene_name in long_reads_isoform_region_weight_matrix_dict[chr_name]:
                        if len(long_reads_isoform_region_weight_matrix_dict[chr_name][gene_name]) > 0:
                            matrix_dict = create_new_matrix_dict(matrix_dict,long_reads_isoform_region_weight_matrix_dict[chr_name][gene_name])
            # matrix_dict['condition_number']  = get_condition_number(matrix_dict['isoform_region_matrix'])
            gene_matrix_dict[chr_name][gene_name] = matrix_dict

    return gene_matrix_dict