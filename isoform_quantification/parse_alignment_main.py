from parse_alignment import map_read,parse_read_line
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import concurrent.futures
import numpy as np
import time
import random
import multiprocessing as mp
import os
# from memory_profiler import profile
# def parse_alignment_iteration(alignment_file_path,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
#                     start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#                     READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,map_f,line_nums):
def parse_alignment_iteration(alignment_file_path,READ_LEN, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,aln_line_marker):
    os.nice(10)
    start_file_pos,num_lines = aln_line_marker
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_count = {}
        local_gene_regions_read_length = {}
        aln_file.seek(start_file_pos)
        line_num_ct = 0
        max_buffer_size = 1e3
        buffer_size = 0
        for line in aln_file:
            if line_num_ct >= num_lines:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'):
                    continue
                aln_line = parse_read_line(line, READ_LEN)
                mapping = map_f(points_dict,interval_tree_dict, filtered_gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
                if (mapping['read_mapped']):
                    for mapping_area in [random.choice(mapping['mapping_area'])]:
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                        if rname not in local_gene_regions_read_count:
                            local_gene_regions_read_count[rname],local_gene_regions_read_length[rname] = {},{}
                        if gname not in local_gene_regions_read_count[rname]:
                            local_gene_regions_read_count[rname][gname],local_gene_regions_read_length[rname][gname] = {},{}
                        if region_name not in local_gene_regions_read_count[rname][gname]:
                            local_gene_regions_read_count[rname][gname][region_name],local_gene_regions_read_length[rname][gname][region_name] = 0,[]
                        local_gene_regions_read_count[rname][gname][region_name] += 1 
                        local_gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
                    buffer_size += 1
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
            if buffer_size > max_buffer_size:
                temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length))
                local_gene_regions_read_count,local_gene_regions_read_length = {},{}
                buffer_size = 0
            line_num_ct += 1
        if buffer_size > 0:
            temp_queue.put((local_gene_regions_read_count,local_gene_regions_read_length))
    return 
def mapping_listener(temp_queue,gene_regions_read_count,gene_regions_read_length):
    num_mapped_lines = 0
    num_lines = 0
    while True:
        msg = temp_queue.get()
        if msg == 'kill':
            break
        else:
            local_gene_regions_read_count,local_gene_regions_read_length = msg
            for chr in local_gene_regions_read_count:
                for gene in local_gene_regions_read_count[chr]:
                    for region in local_gene_regions_read_count[chr][gene]:
                        num_mapped_lines += local_gene_regions_read_count[chr][gene][region]
                        gene_regions_read_count[chr][gene][region] += local_gene_regions_read_count[chr][gene][region]
                        gene_regions_read_length[chr][gene][region] += local_gene_regions_read_length[chr][gene][region]

            # for mapping in local_all_mappings:
            #     num_lines += 1
            #     if len(mapping['gene_candidates'])>0:
            #         num_mapped_to_gene += 1
            #     if (mapping['read_mapped']):
            #         num_mapped_lines += 1
            #         for mapping_area in [random.choice(mapping['mapping_area'])]:
            #             rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
            #             if region_name in gene_regions_read_count[rname][gname]:
            #                 gene_regions_read_count[rname][gname][region_name] += 1 
            #                 gene_regions_read_length[rname][gname][region_name].append(mapping['read_length'])
            #         read_lens.append(mapping['read_length'])
            #         read_names.update(local_read_names)
    return gene_regions_read_count,gene_regions_read_length,num_mapped_lines

# @profile
def parse_alignment(alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,gene_regions_dict,genes_regions_len_dict,long_read,threads):
    start_t = time.time()
    manager = mp.Manager()
    gene_regions_read_count = {}
    gene_regions_read_length ={}
    global filtered_gene_regions_dict
    filtered_gene_regions_dict = defaultdict(lambda: defaultdict(dict))
    for chr in gene_regions_dict:
        gene_regions_read_count[chr],gene_regions_read_length[chr] = {},{}
        for gene in gene_regions_dict[chr]:
            gene_regions_read_count[chr][gene],gene_regions_read_length[chr][gene] = {},{}
            # if (not long_read):
            #     per_gene_regions_dict = filter_regions(gene_regions_dict[chr][gene],long_read=False)
            # else:
            #     per_gene_regions_dict =  filter_regions(gene_regions_dict[chr][gene],long_read=True)
            per_gene_regions_dict =  gene_regions_dict[chr][gene]
            for region in per_gene_regions_dict:
                gene_regions_read_count[chr][gene][region] = 0
                gene_regions_read_length[chr][gene][region] = []
                filtered_gene_regions_dict[chr][gene][region] = True
    if alignment_file_path == None:
        return gene_regions_read_count,150,0
    # Create sorted end and start positions
    global start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST
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
    global points_dict,interval_tree_dict
    points_dict,interval_tree_dict = gene_points_dict,gene_interval_tree_dict
    map_f = map_read
    with open(alignment_file_path, 'r') as aln_file:
        line_offset = []
        offset = 0
        for line in aln_file:
            if line[0] != '@':
                line_offset.append(offset)
            offset += len(line)
    num_aln_lines = len(line_offset)
    pool = mp.Pool(threads+1)
    # partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path)
    chunksize, extra = divmod(num_aln_lines, threads)
    if extra:
        chunksize += 1
    aln_line_marker = []
    # def split_max_size(num):
    #     if num <= 2**32 -1:
    #         return [num]
    #     if num > 2**32 -1:
    #         return [2**32-1] + split_max_size(num-2**32+1)
    for i in range(threads):
        # if (line_offset[i*chunksize] > 2**31-1):
        #     new_line_offset = split_max_size(line_offset[i*chunksize])
        #     aln_line_marker.append((new_line_offset,chunksize))
        # else:
        aln_line_marker.append((line_offset[i*chunksize],chunksize))
    temp_queue = manager.Queue()    
    watcher = pool.apply_async(mapping_listener, args=(temp_queue,gene_regions_read_count,gene_regions_read_length))
    partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path,READ_LEN, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue)
    futures = []
    for marker in aln_line_marker:
        futures.append(pool.apply_async(partial_read_alignment,(marker,)))
    # with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #     futures = [executor.submit(partial_read_alignment,marker) for marker in aln_line_marker]
    #     concurrent.futures.wait(futures)
    for future in futures:
        future.get()
    temp_queue.put('kill')
    gene_regions_read_count,gene_regions_read_length,num_mapped_lines = watcher.get()
    pool.close()
    pool.join()
    read_lens = []
    for chr in gene_regions_read_length:
        for gene in gene_regions_read_length[chr]:
            for region in gene_regions_read_length[chr][gene]:
                read_lens += gene_regions_read_length[chr][gene][region]
    # MAX_BUFF_LINES_IN_MEM = 1e5
    # with open(alignment_file_path, 'r') as aln_file:
    #     aln_lines = []
    #     for line in aln_file:
    #          # parse SAM files
    #         try:
    #             if line[0] == '@':
    #                 continue
    #             fields = line.split('\t')
    #             if (fields[2] == '*'):
    #                 continue
    #             [read_name, read_start_pos, rname, read_len_list] = parse_read_line(line, READ_LEN)
    #             aln_lines.append([read_name, read_start_pos, rname, read_len_list])
    #             read_names.add(fields[0])
    #         except Exception as e:
    #             tb = traceback.format_exc()
    #             print(Exception('Failed to on ' + line, tb))
    #             continue
    #         aln_lines = list(range(MAX_BUFF_LINES_IN_MEM+1))
    #         # if not exceed buffer
    #         if (len(aln_lines) < MAX_BUFF_LINES_IN_MEM):
    #             continue
    #         else: 
    #             print('Buffer exceeded.Starting mapping...')
    #             # if exceed buffer
    #             if (long_read):
    #                 map_f = map_long_read
    #             else:
    #                 map_f = map_read   
    #             partial_map_read = partial(map_f,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
    #                 start_pos_list, start_gname_list, end_pos_list, end_gname_list,
    #                 READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)      
    #             if threads > 1:
    #                 # parallel
    #                 chunksize, extra = divmod(len(aln_lines), threads)
    #                 if extra:
    #                     chunksize += 1
    #                 with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
    #                     all_mappings = executor.map(partial_map_read, aln_lines,chunksize=chunksize)
    #                     # all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
    #             else:
    #                 # serial
    #                 all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
    #             for mapping in all_mappings:
    #                 if (mapping['read_mapped']):
    #                     for mapping_area in [random.choice(mapping['mapping_area'])]:
    #                         rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
    #                         if region_name in gene_regions_read_count[rname][gname]:
    #                             gene_regions_read_count[rname][gname][region_name] += 1
    #                     read_lens.append(mapping['read_length'])
    #             # for testing
    #             # pickle.dump(all_mappings,mapping_pkl_file)
    #             aln_lines = []
    #             print('Done.')
    #     if (len(aln_lines)>0):
    #         if (long_read):
    #             map_f = map_long_read
    #         else:
    #             map_f = map_read   
    #         partial_map_read = partial(map_f,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
    #             start_pos_list, start_gname_list, end_pos_list, end_gname_list,
    #             READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)     
    #         all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
    #         for mapping in all_mappings:
    #             if (mapping['read_mapped']):
    #                 for mapping_area in [random.choice(mapping['mapping_area'])]:
    #                     rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
    #                     if region_name in gene_regions_read_count[rname][gname]:
    #                         gene_regions_read_count[rname][gname][region_name] += 1
    #                 read_lens.append(mapping['read_length'])
    #         # for testing
    #         # pickle.dump(all_mappings,mapping_pkl_file)
    #         aln_lines = []
    # mapping_pkl_file.close()
        
    if (long_read):
        for chr in gene_regions_read_length.copy():
            for gene in gene_regions_read_length[chr].copy():
                region_lens = []
                for region in gene_regions_read_length[chr][gene]:
                    if gene_regions_read_count[chr][gene][region] > 0:
                        region_lens.append(gene_regions_read_length[chr][gene][region])
                if len(region_lens) == 0:
                    del gene_regions_read_length[chr][gene]
                    del gene_regions_read_count[chr][gene]
                else:
                    min_region_len = min(region_lens)
                    for region in gene_regions_read_length[chr][gene].copy():
                        if gene_regions_read_length[chr][gene][region] < min_region_len:
                            del gene_regions_read_length[chr][gene][region]
                            del gene_regions_read_count[chr][gene][region]
                    if len(gene_regions_read_length[chr][gene]) == 0:
                        del gene_regions_read_length[chr][gene]
                        del gene_regions_read_count[chr][gene]
            if len(gene_regions_read_length[chr]) == 0:
                del gene_regions_read_length[chr]
                del gene_regions_read_count[chr]
        return gene_regions_read_count,gene_regions_read_length,sum(read_lens),num_mapped_lines
    else:
        SR_read_len = 150
        return gene_regions_read_count,SR_read_len,num_mapped_lines