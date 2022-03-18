from collections import defaultdict
import pysam
import pandas as pd
import time
import dill as pickle
import shutil
from pathlib import Path
from util import get_long_read_M_dist, get_filtered_out_long_read_M_dist, get_very_short_isoforms
from construct_feature_matrix import generate_all_feature_matrix_short_read
from construct_long_reads_feature_matrix import generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation, process_annotation_for_alignment
from parse_alignment_main import parse_alignment
from generate_output import generate_TransELS_output
from quantification import quantification
from SR_external_quantification import SR_external_quantification
def infer_read_len(short_read_alignment_file_path):
    with pysam.AlignmentFile(short_read_alignment_file_path, "r") as samfile:
            for read in samfile:
                READ_LEN = read.infer_query_length()
                if READ_LEN is not None:
                    break
    return READ_LEN
def parse(ref_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN, short_read_alignment_file_path, threads):
    print('Start parsing annotation...')
    start_time = time.time()
    if short_read_alignment_file_path is not None:
        READ_LEN = infer_read_len(short_read_alignment_file_path)
    gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict = parse_reference_annotation(
        ref_file_path, threads, READ_LEN, READ_JUNC_MIN_MAP_LEN, 'read_length')
    _, gene_range, gene_interval_tree_dict = process_annotation_for_alignment(
        gene_exons_dict, gene_points_dict)
    end_time = time.time()
    print('Done in %.3f s' % (end_time-start_time), flush=True)
    return gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict, gene_range, gene_interval_tree_dict
def map_short_reads(short_read_alignment_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN, gene_isoforms_dict, gene_points_dict, gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, output_path, multi_mapping_filtering, threads):
    print('Mapping short read to regions...', flush=True)
    start_time = time.time()
    if short_read_alignment_file_path is not None:
        if multi_mapping_filtering == 'unique_only':
            pysam.view('-F', '2816', '-q', '10', '-@', f'{threads}', '-h', '-o',
                       f'{output_path}/temp_sr.sam', short_read_alignment_file_path, catch_stdout=False)
            short_read_alignment_file_path = f'{output_path}/temp_sr.sam'
        elif multi_mapping_filtering == 'best':
            pysam.view('-F', '2816', '-@', f'{threads}', '-h', '-o',
                       f'{output_path}/temp_sr.sam', short_read_alignment_file_path, catch_stdout=False)
            short_read_alignment_file_path = f'{output_path}/temp_sr.sam'
    short_read_gene_regions_read_count, _, num_SRs = parse_alignment(short_read_alignment_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN, gene_points_dict,
                                                                               gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, False, False, threads)
    SR_read_len = READ_LEN
    print('Mapped {} short reads with read length {}'.format(
        num_SRs, SR_read_len), flush=True)
    end_time = time.time()
    print('Constructing matrix and calculating condition number...', flush=True)
    short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(
        gene_isoforms_dict, SR_gene_regions_dict, short_read_gene_regions_read_count, SR_read_len, SR_genes_regions_len_dict, num_SRs)
    print('Done in %.3f s' % (end_time-start_time), flush=True)
    return short_read_gene_matrix_dict, SR_read_len
def map_long_reads(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict):
    print('Mapping long read to regions...',flush=True)
    start_time = time.time()
    if multi_mapping_filtering == 'unique_only':
        pysam.view('-F','2820','-q','10','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    elif multi_mapping_filtering == 'best':
        pysam.view('-F','2820','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length = parse_alignment(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
    try:
        Path(f'{output_path}/temp_sr.sam').unlink()
    except:
        pass
    try:
        Path(f'{output_path}/temp_lr.sam').unlink()
    except:
        pass
    print('Mapped {} long reads'.format(num_LRs,flush=True))
    print('Constructing matrix and calculating condition number...',flush=True)
    long_read_gene_matrix_dict = generate_all_feature_matrix_long_read(gene_isoforms_dict, LR_gene_regions_dict, long_read_gene_regions_read_count, long_read_gene_regions_read_length,
                                                                       LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, num_LRs, total_long_read_length, READ_JUNC_MIN_MAP_LEN, output_path, threads)
    end_time = time.time()
    print('Done in %.3f s'%(end_time-start_time),flush=True)
    return long_read_gene_matrix_dict
def get_info_dict_list(gene_isoforms_dict,gene_exons_dict,raw_gene_exons_dict,raw_isoform_exons_dict,gene_isoforms_length_dict):
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
    return info_dict_list
def generate_training_dict(list_of_all_genes_chrs,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,output_path):
    training_dict = {}
    for gene_name,chr_name in list_of_all_genes_chrs:
        if chr_name not in training_dict:
            training_dict[chr_name] = {}
        single_dict = {}
        single_dict['sr_A'] = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix']
        single_dict['sr_b'] = short_read_gene_matrix_dict[chr_name][gene_name]['region_abund_matrix']
        single_dict['lr_A'] = long_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix']
        single_dict['lr_b'] = long_read_gene_matrix_dict[chr_name][gene_name]['region_abund_matrix']
        single_dict['SR_tpm'] = gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_tpm']
        single_dict['LR_tpm'] = gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_tpm']
        single_dict['isoform_names_indics'] = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']
        training_dict[chr_name][gene_name] = single_dict
    with open(f'{output_path}/training.pkl','wb') as f:
        pickle.dump(training_dict,f)
    print('DONE')
def TransELS(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,alpha,beta,P,filtering,multi_mapping_filtering='best',SR_quantification_option='Mili',SR_fastq_list=[],reference_genome='',training=False,DL_model='',assign_unique_mapping_option='',threads=1,READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15):
    print(alpha)
    Path(output_path).mkdir(parents=True, exist_ok=True)
    print('Preprocessing...',flush=True)
    _,gene_points_dict,gene_isoforms_dict,\
        SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,\
            gene_isoforms_length_dict,raw_isoform_exons_dict,_,\
                same_structure_isoform_dict,removed_gene_isoform_dict,gene_range,gene_interval_tree_dict = \
                    parse(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,short_read_alignment_file_path,threads)
    short_read_gene_matrix_dict,SR_read_len = map_short_reads(short_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,gene_isoforms_length_dict,output_path,multi_mapping_filtering,threads)
    long_read_gene_matrix_dict = map_long_reads(long_read_alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict)
    

    # get_very_short_isoforms(output_path,filtered_gene_regions_read_length,LR_gene_regions_dict,isoform_length_dict)
    # generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list)
    # unique_dist_df.to_csv('{}/lr_M_unique_dist.tsv'.format(output_path),sep='\t',index=False)
    # multi_dist_df.to_csv('{}/lr_M_multi_dist.tsv'.format(output_path),sep='\t',index=False)
    # filtered_unique_dist_df.to_csv('{}/filtered_out_lr_M_unique_dist.tsv'.format(output_path),sep='\t',index=False)
    # filtered_multi_dist_df.to_csv('{}/filtered_out_lr_M_multi_dist.tsv'.format(output_path),sep='\t',index=False)
    # with open('{}/lr.pkl'.format(output_path),'wb') as f:
    #     pickle.dump((gene_isoforms_length_dict,long_read_gene_regions_read_length,long_read_gene_regions_read_count,LR_gene_regions_dict,filtered_gene_regions_read_length),f)
    
    # info_dict_list = get_info_dict_list(gene_isoforms_dict,gene_exons_dict,raw_gene_exons_dict,raw_isoform_exons_dict,gene_isoforms_length_dict)
    print('Start quantification...',flush=True)
    start_time = time.time()
    SR_gene_isoform_expression_dict = None
    if SR_quantification_option != 'Mili':
        if short_read_alignment_file_path is not None:
            ref_genome = reference_genome
            SR_gene_isoform_expression_dict = SR_external_quantification(short_read_gene_matrix_dict,gene_isoforms_length_dict,SR_quantification_option,SR_fastq_list,SR_read_len,ref_file_path,ref_genome,output_path,threads)
    import dill as pickle
    # pickle.dump((long_read_gene_matrix_dict,gene_points_dict,LR_gene_regions_dict,short_read_gene_matrix_dict),open(f'{output_path}/dict.pkl','wb'))
    gene_isoform_tpm_expression_dict,list_of_all_genes_chrs = quantification(short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoforms_length_dict,SR_gene_isoform_expression_dict,SR_quantification_option,DL_model,alpha,beta,P,assign_unique_mapping_option)
    end_time = time.time()
    print('Done in %.3f s'%(end_time-start_time),flush=True)
    # import dill as pickle
    # rep_name = output_path.split('/')[-2]
    # # rep_name = 1
    with open(f'{output_path}/dict.pkl','wb') as f:
        pickle.dump([long_read_gene_matrix_dict,gene_points_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict],f)
    print('Generating output...',flush=True)
    if training:
        generate_training_dict(list_of_all_genes_chrs,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,output_path)
    generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict,same_structure_isoform_dict,removed_gene_isoform_dict)
    try:
        shutil.rmtree(f'{output_path}/temp/')
    except:
        pass
    print('Done in %.3f s'%(end_time-start_time),flush=True)