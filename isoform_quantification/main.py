from parse_annotation import parse_annotation
from construct_feature_matrix import generate_all_feature_matrix,calculate_all_condition_number,filter_multi_exons_regions
from parse_alignment import map_read,map_long_read
from operator import itemgetter, attrgetter
from collections import defaultdict

import numpy as np
from numpy import linalg as LA
from qpsolvers import solve_qp

import datetime
import traceback
import argparse
import time
from pathlib import Path
def parse_reference_annotation(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    print('Parsing the reference annotation...')
    start = time.time()
    [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict] = parse_annotation(ref_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN)
    for chr_name in gene_regions_dict:
        for gene_name in gene_regions_dict[chr_name].copy():
            if len(filter_multi_exons_regions(gene_regions_dict[chr_name][gene_name])) == 0:
                del gene_points_dict[chr_name][gene_name],gene_isoforms_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],raw_isoform_exons_dict[chr_name][gene_name]
    print('Done in {}s !'.format(time.time()-start))
    return gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict

def process_annotation_for_alignment(gene_points_dict,gene_regions_dict):
    print('Process the reference annotation for read_count...')
    # process the reference annotation for read_count TODO: a universal format
    gene_regions_points_list,gene_range = dict(),dict()
    for chr_name in gene_points_dict:
        gene_regions_points_list[chr_name],gene_range[chr_name] = dict(),list()
        for gene_name in gene_points_dict[chr_name]:
            points = list(gene_points_dict[chr_name][gene_name].keys())
            gene_regions_points_list[chr_name][gene_name] = points
            gene_range[chr_name].append([gene_name,points[0],points[-1]])
    return gene_regions_points_list,gene_range

def parse_alignment(alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict,long_read):
    read_names = set()
    gene_regions_read_count = defaultdict(lambda: defaultdict(dict))
    for chr in gene_regions_dict:
        for gene in gene_regions_dict[chr]:
            for region in gene_regions_dict[chr][gene]:
                gene_regions_read_count[chr][gene][region] = 0
    # Create sorted end and start positions
    start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST = dict(),dict(),dict(),dict(),list(gene_range.keys())
    CHR_LIST = list(gene_range.keys())
    for chr in CHR_LIST:     
        # Sort based on start position
        temp_list = sorted(gene_range[chr], key=itemgetter(1))
        start_pos_list[chr] = [temp_list[j][1] for j in range(len(temp_list))]
        start_gname_list[chr] = [temp_list[j][0] for j in range(len(temp_list))]
        # Sort based on end position
        temp_list = sorted(gene_range[chr], key=itemgetter(2))
        end_pos_list[chr] = [temp_list[j][2] for j in range(len(temp_list))]
        end_gname_list[chr] = [temp_list[j][0] for j in range(len(temp_list))]
                       
    # for each read in the alignment
    with open(alignment_file_path, 'r') as f:
        for line in f.readlines():
            try:
                is_legal_alignment = True
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[1] == 4):
                    continue
                illegal_ops = ['S','H','P','=','X']
                for op in illegal_ops:
                    if op in fields[5]:
                        is_legal_alignment = False
                if (is_legal_alignment):
                    if (not long_read):
                        mapping = map_read(line, gene_regions_read_count, gene_regions_points_list, 
                            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
                    else:
                        mapping = map_long_read(line, gene_regions_read_count, gene_regions_points_list, 
                            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
                    read_names.add(fields[0])
            except Exception as e:
                tb = traceback.format_exc()
                continue
                raise Exception('Failed to on ' + line, tb)
    return gene_regions_read_count,len(read_names)
def generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    with open(output_path+"/kvalues.out",'w') as f:
        f.write('Gene\tChr\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\t\n')
        for chr_name in short_read_gene_matrix_dict:
            for gene_name in short_read_gene_matrix_dict[chr_name]:
                SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name,chr_name,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
    
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
def estimate_isoform_expression(SR_isoform_region_matrix,SR_region_read_count_matrix,LR_isoform_region_matrix,LR_region_read_count_matrix,isoform_lengths,weight,beta):
    num_isoforms = SR_isoform_region_matrix.shape[1]
    Q = 2 * (1.0 - weight) * np.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * weight * np.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta * np.identity(num_isoforms)
    c = -2 * (1.0 - weight) * np.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * weight * np.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T)
    lb = np.zeros(num_isoforms)
    isoform_expression = solve_qp(Q, c, lb = lb)
    if ((isoform_expression+1e-5<0).any()):
        raise ValueError('Obtain negative value for isoform expression')
    isoform_expression[isoform_expression<0] = 0
    def div0( a, b ):
        with np.errstate(divide='ignore', invalid='ignore'):
            c = np.true_divide(a, b)
            c[ ~ np.isfinite( c )] = 0
        return c
    isoform_expression = div0(isoform_expression,isoform_lengths)
    return isoform_expression

def generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,gene_isoform_fpkm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    with open(output_path+"/expression_gene.out",'w') as f_gene:
        with open(output_path+"/expression_isoform.out",'w') as f_isoform:
            f_gene.write('Gene\tChr\tTPM\tFPKM\n')
            f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tFPKM\n')
            for chr_name in short_read_gene_matrix_dict:
                for gene_name in short_read_gene_matrix_dict[chr_name]:
                    tpm_sum = fpkm_sum = 0
                    for isoform_name in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']:
                        start_pos = min(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                        end_pos = max(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'])
                        isoform_len = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                        isoform_index = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics'][isoform_name]
                        tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]
                        fpkm = gene_isoform_fpkm_expression_dict[chr_name][gene_name][isoform_index]
                        tpm_sum += tpm
                        fpkm_sum += fpkm
                        f_isoform.write('%s\t%s\t%s\t%d\t%d\t%d\t%f\t%f\n'%(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,fpkm))
                    f_gene.write('%s\t%s\t%f\t%f\n'%(gene_name,chr_name,tpm_sum,fpkm_sum))
def normalize_expression(gene_isoform_expression_dict,total_num_reads,isoform_expression_sum):
    gene_isoform_tpm_expression_dict = defaultdict(dict)
    gene_isoform_fpkm_expression_dict = defaultdict(dict)
    for chr_name in gene_isoform_expression_dict:
        for gene_name in gene_isoform_expression_dict[chr_name]:
            gene_isoform_tpm_expression_dict[chr_name][gene_name] = gene_isoform_expression_dict[chr_name][gene_name] * 1e6 / isoform_expression_sum
            gene_isoform_fpkm_expression_dict[chr_name][gene_name] = gene_isoform_expression_dict[chr_name][gene_name] * 1e9 / total_num_reads
            
    return gene_isoform_tpm_expression_dict,gene_isoform_fpkm_expression_dict
            
def TransELS(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,weight = 0.01,beta = 0.01,READ_LEN=79,READ_JUNC_MIN_MAP_LEN=10):
    gene_points_dict,gene_isoforms_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = parse_reference_annotation(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    gene_regions_points_list,gene_range = process_annotation_for_alignment(gene_points_dict,gene_regions_dict)
    print('Mapping short read to regions...')  
    short_read_gene_regions_read_count,num_SRs = parse_alignment(short_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict,long_read = False)
    print('Mapping long read to regions...')  
    long_read_gene_regions_read_count,num_LRs = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict, long_read = True)
    print('Constructing matrix and calculating condition number...')
    short_read_gene_matrix_dict = generate_all_feature_matrix(gene_isoforms_dict,gene_regions_dict,short_read_gene_regions_read_count,genes_regions_len_dict,allow_multi_exons=False)
    long_read_gene_matrix_dict = generate_all_feature_matrix(gene_isoforms_dict,gene_regions_dict,long_read_gene_regions_read_count,genes_regions_len_dict,allow_multi_exons=True)
    chr_list = set(short_read_gene_matrix_dict.keys()).intersection(set(long_read_gene_matrix_dict.keys()))
    gene_isoform_expression_dict = defaultdict(dict)
    print('Calculating the isoform expression...')
    isoform_expression_sum = 0
    for chr_name in chr_list:
        gene_list = set(short_read_gene_matrix_dict[chr_name].keys()).intersection(set(long_read_gene_matrix_dict[chr_name].keys()))
        for gene_name in gene_list:
            isoform_lengths = np.zeros((len(gene_isoforms_length_dict[chr_name][gene_name])))
            isoform_names_indics = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']
            for isoform_name in isoform_names_indics:
                isoform_lengths[isoform_names_indics[isoform_name]] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                try:
                    gene_isoform_expression_dict[chr_name][gene_name] = estimate_isoform_expression(
                                            short_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix'],
                                            short_read_gene_matrix_dict[chr_name][gene_name]['region_read_count_matrix'],
                                            long_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix'],
                                            long_read_gene_matrix_dict[chr_name][gene_name]['region_read_count_matrix'],isoform_lengths,weight,beta)
                    isoform_expression_sum += sum(gene_isoform_expression_dict[chr_name][gene_name])
                except ValueError as e:
                    print("Error encountered for gene %s chr %s:"%((chr_name,gene_name)) + e)
    gene_isoform_tpm_expression_dict,gene_isoform_fpkm_expression_dict = normalize_expression(gene_isoform_expression_dict,num_SRs+num_LRs,isoform_expression_sum)
    generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,gene_isoform_fpkm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict)
    return short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,gene_isoform_fpkm_expression_dict
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_TrEESR = subparsers.add_parser('TrEESR', help='TrEESR')
    parser_TransELS = subparsers.add_parser('TransELS', help='TransELS')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for TrEESR')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for TrEESR')
    requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=True)
    requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    optional_TransELS.add_argument('--alpha',type=float,default=0.1, help="Alpha")
    optional_TransELS.add_argument('--beta',type=float, default=0.1,help="Beta")
#     optional.add_argument(
#     '-h',
#     '--help',
#     action='help',
#     default=argparse.SUPPRESS,
#     help='show this help message and exit')
    
    args = parser.parse_args()
    if args.subparser_name == 'TrEESR':
        print('Using TrEESR')
        TrEESR(args.gtf_annotation_path,args.output_path)
    else:
        print('Using TransELS')
        TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,args.alpha,args.beta)
if __name__ == "__main__":
    parse_arguments()
