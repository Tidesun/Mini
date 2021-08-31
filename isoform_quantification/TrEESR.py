from collections import defaultdict
from construct_feature_matrix import calculate_all_condition_number
from parse_annotation_main import parse_reference_annotation
from generate_output import generate_TrEESR_output
from construct_feature_matrix import generate_all_feature_matrix_short_read,generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation,process_annotation_for_alignment
from parse_alignment_main import parse_alignment
import datetime
def TrEESR(ref_file_path,output_path,long_read_alignment_file_path,threads,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=0):
    start_time = datetime.datetime.now()
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict = \
        parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,None)
    end_time_1 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_1-start_time).total_seconds()),flush=True)
    print('Calculating the condition number...',flush=True)
    short_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,SR_gene_regions_dict,allow_multi_exons=False)
    gene_regions_points_list,gene_range,gene_interval_tree_dict = process_annotation_for_alignment(gene_exons_dict,gene_points_dict)
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,True,threads)
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict,LR_gene_regions_dict,long_read_gene_regions_read_count,long_read_gene_regions_read_length,LR_genes_regions_len_dict,num_LRs,total_long_read_length,'original')
    # long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,LR_gene_regions_dict,allow_multi_exons=True)
    end_time_2 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_2-end_time_1).total_seconds()),flush=True)
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
    generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list)
    return short_read_gene_matrix_dict,long_read_gene_matrix_dict
def get_kvalues_dict(ref_file_path,threads,READ_LEN=150,READ_JUNC_MIN_MAP_LEN=10):
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=True)
    return long_read_gene_matrix_dict