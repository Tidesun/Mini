from construct_feature_matrix import calculate_all_condition_number
from parse_annotation_main import parse_reference_annotation
from generate_output import generate_TrEESR_output
import datetime

def TrEESR(ref_file_path,output_path,READ_LEN=79,READ_JUNC_MIN_MAP_LEN=10):
    start_time = datetime.datetime.now()
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    end_time_1 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_1-start_time).total_seconds()))
    print('Calculating the condition number...')
    short_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=False)
    long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=True)
    end_time_2 = datetime.datetime.now()
    print('Done in %.3f s'%((end_time_2-end_time_1).total_seconds()))
    generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict)
    return short_read_gene_matrix_dict,long_read_gene_matrix_dict
def get_kvalues_dict(ref_file_path,READ_LEN=79,READ_JUNC_MIN_MAP_LEN=10):
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    long_read_gene_matrix_dict = calculate_all_condition_number(gene_isoforms_dict,gene_regions_dict,allow_multi_exons=True)
    return long_read_gene_matrix_dict