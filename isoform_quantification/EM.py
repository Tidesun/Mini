from collections import defaultdict
import pysam
import time
import gc
from pathlib import Path
from construct_feature_matrix import generate_all_feature_matrix_short_read
from construct_long_reads_feature_matrix import generate_all_feature_matrix_long_read
from parse_annotation_main import parse_reference_annotation, process_annotation_for_alignment
from parse_alignment_main import parse_alignment
from EM_libraries.EM_algo import EM_algo_main
from EM_theta_iterative_libraries.EM_algo import EM_algo_theta_iter_main
from EM_kde_libraries.EM_algo import EM_algo_kde_main
from EM_kde_score_libraries.EM_algo import EM_algo_kde_score_main
from EM_SR.EM_SR import EM_algo_SR
import config
def infer_read_len(short_read_alignment_file_path):
    READ_LEN = 150
    sr_sam_valid = False
    with pysam.AlignmentFile(short_read_alignment_file_path, "r") as samfile:
            for read in samfile:
                READ_LEN = read.infer_query_length()
                if READ_LEN is not None:
                    sr_sam_valid = True
                    break
    if not sr_sam_valid:
        print('The SR sam seems invalid. Please check!')
        exit()
    return READ_LEN
def parse(ref_file_path, READ_JUNC_MIN_MAP_LEN, short_read_alignment_file_path, threads):
    print('Start parsing annotation...')
    start_time = time.time()
    if short_read_alignment_file_path is not None:
        READ_LEN = infer_read_len(short_read_alignment_file_path)
    else:
        READ_LEN = 150
    gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict = parse_reference_annotation(
        ref_file_path, threads, READ_LEN, READ_JUNC_MIN_MAP_LEN, 'real_data')
    _, gene_range, gene_interval_tree_dict = process_annotation_for_alignment(
        gene_exons_dict, gene_points_dict)
    end_time = time.time()
    print('Done in %.3f s' % (end_time-start_time), flush=True)
    return gene_exons_dict, gene_points_dict, gene_isoforms_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, LR_gene_regions_dict, LR_genes_regions_len_dict, gene_isoforms_length_dict, raw_isoform_exons_dict, raw_gene_exons_dict, same_structure_isoform_dict, removed_gene_isoform_dict, gene_range, gene_interval_tree_dict
def map_short_reads(short_read_alignment_file_path, READ_JUNC_MIN_MAP_LEN, gene_isoforms_dict, gene_points_dict, gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, output_path, multi_mapping_filtering, threads):
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
    short_read_gene_regions_read_count, SR_read_len, num_SRs = parse_alignment(short_read_alignment_file_path, READ_JUNC_MIN_MAP_LEN, gene_points_dict,
                                                                               gene_range, gene_interval_tree_dict, SR_gene_regions_dict, SR_genes_regions_len_dict, gene_isoforms_length_dict, False, False, threads)
    config.READ_LEN = SR_read_len
    print('Mapped {} short reads with read length {}'.format(
        num_SRs, SR_read_len), flush=True)
    end_time = time.time()
    print('Constructing matrix and calculating condition number...', flush=True)
    short_read_gene_matrix_dict = generate_all_feature_matrix_short_read(
        gene_isoforms_dict, SR_gene_regions_dict, short_read_gene_regions_read_count, SR_read_len, SR_genes_regions_len_dict, num_SRs)
    print('Done in %.3f s' % (end_time-start_time), flush=True)
    return short_read_gene_matrix_dict, SR_read_len
def map_long_reads(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict):
    print('Mapping long read to regions...',flush=True)
    start_time = time.time()
    if multi_mapping_filtering == 'unique_only':
        pysam.view('-F','2820','-q','10','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    elif multi_mapping_filtering == 'best':
        pysam.view('-F','2820','-@',f'{threads}','-h','-o',f'{output_path}/temp_lr.sam',long_read_alignment_file_path,catch_stdout=False)
        long_read_alignment_file_path = f'{output_path}/temp_lr.sam'
    long_read_gene_regions_read_count,long_read_gene_regions_read_length,total_long_read_length,num_LRs,filtered_gene_regions_read_length,gene_regions_read_pos = parse_alignment(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict, True,filtering,threads)
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
    return long_read_gene_matrix_dict,gene_regions_read_pos,long_read_gene_regions_read_length
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
import re
def parse_for_EM_algo(annotation):
    gene_exon_dict = {}
    gene_isoform_dict = {}
    isoform_exon_dict = {}
    strand_dict = {}
    with open(annotation,'r') as f:
        for line in f:
            if line.lstrip()[0] == "#":
                continue
            fields = line.split('\t')
            if (fields[2] != 'exon'):
                continue
            strand = fields[6]
            chr_name = fields[0]
            gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
            isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            if gene_name not in gene_exon_dict:
                gene_exon_dict[gene_name] = []
                gene_isoform_dict[gene_name] = set()
            if isoform_name not in isoform_exon_dict:
                isoform_exon_dict[isoform_name] = []
            gene_exon_dict[gene_name].append([start_pos,end_pos])
            gene_isoform_dict[gene_name].add(isoform_name)
            isoform_exon_dict[isoform_name].append([start_pos,end_pos])
            strand_dict[gene_name] = strand
    for isoform in isoform_exon_dict:
        isoform_exon_dict[isoform] = sorted(isoform_exon_dict[isoform],key=lambda x:(x[0],x[1]))
    isoform_len_dict = {}
    for isoform in isoform_exon_dict:
        isoform_len_dict[isoform] = 0
        for exon in isoform_exon_dict[isoform]:
            isoform_len_dict[isoform] += exon[1] - exon[0] + 1
#     isoform_len_set = set(isoform_len_dict.values())
    return isoform_len_dict,isoform_exon_dict,strand_dict
def EM(ref_file_path,short_read_alignment_file_path,long_read_alignment_file_path,output_path,alpha,beta,P,filtering,multi_mapping_filtering='best',SR_quantification_option='Mili',SR_fastq_list=[],reference_genome='',training=False,DL_model='',assign_unique_mapping_option='',threads=1,READ_LEN=0,READ_JUNC_MIN_MAP_LEN=15,EM_choice='LIQA_modified',iter_theta='True'):
    print(alpha)
    Path(output_path).mkdir(parents=True, exist_ok=True)
    print('Preprocessing...',flush=True)
    _,gene_points_dict,gene_isoforms_dict,\
        _,_,LR_gene_regions_dict,LR_genes_regions_len_dict,\
            gene_isoforms_length_dict,raw_isoform_exons_dict,_,\
                _,_,gene_range,gene_interval_tree_dict = \
                    parse(ref_file_path,READ_JUNC_MIN_MAP_LEN,short_read_alignment_file_path,threads)
    _,gene_regions_read_pos,_ = map_long_reads(long_read_alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_isoforms_dict,gene_points_dict,gene_range,gene_interval_tree_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,filtering,output_path,multi_mapping_filtering,threads,raw_isoform_exons_dict)
    isoform_len_dict,isoform_exon_dict,strand_dict = parse_for_EM_algo(ref_file_path)
    del gene_points_dict
    del gene_isoforms_dict
    del LR_genes_regions_len_dict
    del gene_isoforms_length_dict
    del raw_isoform_exons_dict
    del gene_range
    del gene_interval_tree_dict
    gc.collect()
    print('Start quantification...',flush=True)
    start_time = time.time()
    if EM_choice == 'kde':
        EM_algo_kde_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice,config.kde_path,set(isoform_len_dict.values()))
    elif EM_choice == 'kde_score':
        EM_algo_kde_score_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice,config.kde_path,set(isoform_len_dict.values()))
    else:
        if iter_theta == 'True':
            EM_algo_theta_iter_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice)
        else:
            EM_algo_main(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_pos,LR_gene_regions_dict,threads,output_path,EM_choice)
    end_time = time.time()
    print('Done in %.3f s'%(end_time-start_time),flush=True)
def EM_SR(eff_len_file,short_read_alignment_file_path,output_path,threads):
    start_time = time.time()
    # isoform_len_dict,_,_ = parse_for_EM_algo(ref_file_path)
    EM_algo_SR(short_read_alignment_file_path,eff_len_file,output_path,threads)
    end_time = time.time()
    print('Done in %.3f s'%(end_time-start_time),flush=True)
