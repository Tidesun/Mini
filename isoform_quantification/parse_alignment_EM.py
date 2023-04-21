from parse_alignment import map_read,parse_read_line
from patch_mp import patch_mp_connection_bpo_17560
from parse_annotation_main import check_valid_region
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import numpy as np
import time
import random
import multiprocessing as mp
import os
from util import check_region_type
import config
import pickle
from pathlib import Path
import pandas as pd
import shutil

# from memory_profiler import profile
# def parse_alignment_iteration(alignment_file_path,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
#                     start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#                     READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,map_f,line_nums):
def debuginfoStr(info):
    print(info,flush=True)
    with open('/proc/self/status') as f:
        memusage = f.read().split('VmRSS:')[1].split('\n')[0][:-3]
    mem = int(memusage.strip())/1024
    print('Mem consumption: '+str(mem),flush=True)
def get_reads_isoform_info(output_path,gene_regions_read_mapping,LR_gene_regions_dict):
    with open(f'{output_path}/temp/LR_alignments_dict/isoform_dict','rb') as f:
        [isoform_len_dict,isoform_exon_dict,strand_dict] = pickle.load(f)
    reads_isoform_info = {}
    expression_dict = {}
    unique_mapping_expression_dict = {}
    for rname in gene_regions_read_mapping:
        for gname in gene_regions_read_mapping[rname]:
            for region in gene_regions_read_mapping[rname][gname]:
                for read_mapping in gene_regions_read_mapping[rname][gname][region]:
                    read_start,read_end = read_mapping['read_pos']
                    read_name = read_mapping['read_name']
                    if read_name not in reads_isoform_info:
                        reads_isoform_info[read_name] = {}
                    for isoform in LR_gene_regions_dict[rname][gname][region]: 
                        # assign non-zero expression isoform
                        if isoform not in expression_dict:
                            expression_dict[isoform] = 0
                            unique_mapping_expression_dict[isoform] = 0
                        isoform_len = isoform_len_dict[isoform]
                        start_offset = 0
                        for [exon_start,exon_end] in isoform_exon_dict[isoform]:
                            if exon_end > read_start:
                                start_offset += read_start - exon_start + 1
                                break
                            else:
                                start_offset += exon_end - exon_start + 1
                        end_offset = isoform_len_dict[isoform] - start_offset - read_mapping['read_length']
                        if start_offset < 0:
                            start_offset = 0
                        if end_offset < 0:
                            end_offset = 0
                        if strand_dict[gname]  == '+':
                            reads_isoform_info[read_name][isoform] = [isoform_len,start_offset,end_offset]
                        else:
                            reads_isoform_info[read_name][isoform] = [isoform_len,end_offset,start_offset]
                    if len(LR_gene_regions_dict[rname][gname][region]) == 1:
                         # use unique mapping reads as initial expression
                        isoform = next(iter(LR_gene_regions_dict[rname][gname][region]))
                        expression_dict[isoform] += 1
                        unique_mapping_expression_dict[isoform] += 1
                    else:
                        num_isoforms = len(LR_gene_regions_dict[rname][gname][region])
                        for isoform in LR_gene_regions_dict[rname][gname][region]:
                            expression_dict[isoform] += 1/num_isoforms
    return reads_isoform_info,expression_dict,unique_mapping_expression_dict
def get_read_len_dist(reads_isoform_info):
    read_len_dict = {}
    read_len_dist_dict = {}
    for read in reads_isoform_info:
        [isoform_len,f_end,t_end] = next(iter(reads_isoform_info[read].values()))
        read_len = isoform_len -  f_end - t_end
        read_len_dict[read] = read_len
        if read_len not in read_len_dist_dict:
            read_len_dist_dict[read_len] = 0
        read_len_dist_dict[read_len] += 1
    # read_len_dist = pd.Series(read_len_dist_dict).to_frame()
    # read_len_dist.columns = ['PDF']
    # read_len_dist = read_len_dist.sort_index(ascending=True)
    return read_len_dict,read_len_dist_dict
def parse_alignment_iteration(alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,CHR_LIST,output_path,aln_line_marker):
    os.nice(10)
    Path(f'{output_path}/temp/LR_alignments/').mkdir(exist_ok=True,parents=True)
    start_file_pos,end_file_pos,worker_id = aln_line_marker
    batch_id = 0
    points_dict,interval_tree_dict, gene_regions_dict,\
        start_pos_list, start_gname_list, end_pos_list, end_gname_list = {},{},{},{},{},{},{}
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_pos = {}
        aln_file.seek(start_file_pos)
        max_buffer_size = 1e5
        buffer_size = 0
        while True:
            line = aln_file.readline()
            if aln_file.tell() >= end_file_pos:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'):
                    continue
                aln_line = parse_read_line(line)
                [_, _, rname, _] = aln_line
                if rname not in points_dict:
                    with open(f'{output_path}/temp/LR_alignments_dict/{rname}','rb') as f:
                        [points_dict[rname],interval_tree_dict[rname],gene_regions_dict[rname],start_pos_list[rname],end_pos_list[rname],start_gname_list[rname],end_gname_list[rname]] = pickle.load(f)
                mapping = map_f(points_dict,interval_tree_dict, gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
                if (mapping['read_mapped']):
                    for mapping_area in [random.choice(mapping['mapping_area'])]:
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                        if rname not in local_gene_regions_read_pos:
                            local_gene_regions_read_pos[rname] = {}
                        if gname not in local_gene_regions_read_pos[rname]:
                            local_gene_regions_read_pos[rname][gname] = {}
                        if region_name not in local_gene_regions_read_pos[rname][gname]:
                            local_gene_regions_read_pos[rname][gname][region_name] = []
                        # if long_read:
                        local_gene_regions_read_pos[rname][gname][region_name].append(mapping)
                    buffer_size += 1
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
            if buffer_size > max_buffer_size:
                reads_isoform_info,expression_dict,unique_mapping_expression_dict = get_reads_isoform_info(output_path,local_gene_regions_read_pos,gene_regions_dict)
                read_len_dict,read_len_dist_dict = get_read_len_dist(reads_isoform_info)
                with open(f'{output_path}/temp/LR_alignments/{worker_id}_{batch_id}','wb') as f:
                    pickle.dump([reads_isoform_info,read_len_dict,read_len_dist_dict,expression_dict,unique_mapping_expression_dict],f)
                batch_id += 1
                del reads_isoform_info
                del expression_dict
                del unique_mapping_expression_dict
                local_gene_regions_read_pos = {}
                buffer_size = 0
        if buffer_size > 0:
            reads_isoform_info,expression_dict,unique_mapping_expression_dict =get_reads_isoform_info(output_path,local_gene_regions_read_pos,gene_regions_dict)
            read_len_dict,read_len_dist_dict = get_read_len_dist(reads_isoform_info)
            with open(f'{output_path}/temp/LR_alignments/{worker_id}_{batch_id}','wb') as f:
                pickle.dump([reads_isoform_info,read_len_dict,read_len_dist_dict,expression_dict,unique_mapping_expression_dict],f)
    return 
# @profile
# def get_aln_line_marker(alignment_file_path,threads):
#     with open(alignment_file_path, 'r') as aln_file:
#         line_offset = []
#         offset = 0
#         for line in aln_file:
#             if line[0] != '@':
#                 line_offset.append(offset)
#             offset += len(line)
#     num_aln_lines = len(line_offset)
#     chunksize, extra = divmod(num_aln_lines, threads)
#     if extra:
#         chunksize += 1
#     aln_line_marker = []
#     for i in range(threads):
#         aln_line_marker.append((line_offset[i*chunksize],chunksize,i))
#     return aln_line_marker
def get_aln_line_marker(alignment_file_path,threads):
    '''
    Split the sam file into THREADS chunks
    !Split by read
    '''
    file_stats = os.stat(alignment_file_path)
    total_bytes = file_stats.st_size
    chunksize, extra = divmod(total_bytes, threads)
    if extra:
        chunksize += 1
    byte_marker = []
    with open(alignment_file_path,'r') as f:
        for i in range(threads):
            if i == 0:
                start_offset = 0
                for line in f:
                    if line[0] != '@':
                        break
                    start_offset += len(line)
            else:
                f.seek(i*chunksize)
                f.readline()
                start_offset = f.tell()
            byte_marker.append(start_offset)
    byte_marker.append(total_bytes)
    aln_line_marker = []
    for i in range(threads):
         aln_line_marker.append((byte_marker[i],byte_marker[i+1],i))
    return aln_line_marker
def parse_alignment_EM(alignment_file_path,READ_JUNC_MIN_MAP_LEN,output_path,threads,CHR_LIST):
    print(threads)
    patch_mp_connection_bpo_17560()
    # Create sorted end and start positions
    map_f = map_read
    # debuginfoStr('Before MP mem usage')
    pool = mp.Pool(threads)
    # partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path)
    # temp_queue = manager.Queue()    
    # watcher = pool.apply_async(mapping_listener, args=(temp_queue,gene_regions_read_count,gene_regions_read_length,gene_regions_read_pos))
    partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,CHR_LIST,output_path)
    futures = []
    aln_line_marker = get_aln_line_marker(alignment_file_path,threads)
    for marker in aln_line_marker:
        futures.append(pool.apply_async(partial_read_alignment,(marker,)))
    # with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #     futures = [executor.submit(partial_read_alignment,marker) for marker in aln_line_marker]
    #     concurrent.futures.wait(futures)
    for future in futures:
        future.get()
    # temp_queue.put('kill')
    # gene_regions_read_count,gene_regions_read_length,num_mapped_lines,gene_regions_read_pos = watcher.get()
    pool.close()
    pool.join()
    try:
        shutil.rmtree(f'{output_path}/temp/LR_alignments_dict/')
    except:
        pass

        
   