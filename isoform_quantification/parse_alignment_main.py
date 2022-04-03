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
def parse_alignment_iteration(alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,long_read,aln_line_marker):
    os.nice(10)
    start_file_pos,num_lines = aln_line_marker
    with open(alignment_file_path, 'r') as aln_file:
        local_gene_regions_read_count = {}
        local_gene_regions_read_length = {}
        aln_file.seek(start_file_pos)
        line_num_ct = 0
        max_buffer_size = 1e2
        buffer_size = 0
        for line in aln_file:
            line_num_ct += 1
            if line_num_ct > num_lines:
                break
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[2] == '*'):
                    continue
                aln_line = parse_read_line(line)
                mapping = map_f(points_dict,interval_tree_dict, filtered_gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_JUNC_MIN_MAP_LEN, CHR_LIST,aln_line)
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
                        # if long_read:
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
            for rname in local_gene_regions_read_count:
                for gname in local_gene_regions_read_count[rname]:
                    for region in local_gene_regions_read_count[rname][gname]:
                        num_mapped_lines += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_count[rname][gname][region] += local_gene_regions_read_count[rname][gname][region]
                        gene_regions_read_length[rname][gname][region] += local_gene_regions_read_length[rname][gname][region]

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
def get_aln_line_marker(alignment_file_path,threads):
    with open(alignment_file_path, 'r') as aln_file:
        line_offset = []
        offset = 0
        for line in aln_file:
            if line[0] != '@':
                line_offset.append(offset)
            offset += len(line)
    num_aln_lines = len(line_offset)
    chunksize, extra = divmod(num_aln_lines, threads)
    if extra:
        chunksize += 1
    aln_line_marker = []
    for i in range(threads):
        aln_line_marker.append((line_offset[i*chunksize],chunksize))
    return aln_line_marker

def parse_alignment(alignment_file_path,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,gene_regions_dict,genes_regions_len_dict,gene_isoforms_length_dict,long_read,filtering,threads):
    patch_mp_connection_bpo_17560()
    start_t = time.time()
    manager = mp.Manager()
    gene_regions_read_count = {}
    gene_regions_read_length ={}
    global filtered_gene_regions_dict
    filtered_gene_regions_dict = defaultdict(lambda: defaultdict(dict))
    for rname in gene_regions_dict:
        gene_regions_read_count[rname],gene_regions_read_length[rname] = {},{}
        for gname in gene_regions_dict[rname]:
            gene_regions_read_count[rname][gname],gene_regions_read_length[rname][gname] = {},{}
            # if (not long_read):
            #     per_gene_regions_dict = filter_regions(gene_regions_dict[rname][gname],long_read=False)
            # else:
            #     per_gene_regions_dict =  filter_regions(gene_regions_dict[rname][gname],long_read=True)
            per_gene_regions_dict =  gene_regions_dict[rname][gname]
            for region in per_gene_regions_dict:
                gene_regions_read_count[rname][gname][region] = 0
                gene_regions_read_length[rname][gname][region] = []
                filtered_gene_regions_dict[rname][gname][region] = True
    if alignment_file_path == None:
        return gene_regions_read_count,150,0
    # Create sorted end and start positions
    global start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST
    start_pos_list,end_pos_list,start_gname_list,end_gname_list,CHR_LIST = dict(),dict(),dict(),dict(),list(gene_range.keys())
    CHR_LIST = list(gene_range.keys())
    for rname in CHR_LIST:     
        # Sort based on start position
        temp_list = sorted(gene_range[rname], key=itemgetter(1))
        start_pos_list[rname] = [temp_list[j][1] for j in range(len(temp_list))]
        start_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
        # Sort based on end position
        temp_list = sorted(gene_range[rname], key=itemgetter(2))
        end_pos_list[rname] = [temp_list[j][2] for j in range(len(temp_list))]
        end_gname_list[rname] = [temp_list[j][0] for j in range(len(temp_list))]
    global points_dict,interval_tree_dict
    points_dict,interval_tree_dict = gene_points_dict,gene_interval_tree_dict
    map_f = map_read
    debuginfoStr('Before MP mem usage')
    pool = mp.Pool(threads+1)
    # partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path)
    temp_queue = manager.Queue()    
    watcher = pool.apply_async(mapping_listener, args=(temp_queue,gene_regions_read_count,gene_regions_read_length))
    partial_read_alignment = partial(parse_alignment_iteration,alignment_file_path, READ_JUNC_MIN_MAP_LEN,map_f,temp_queue,long_read)
    futures = []
    aln_line_marker = get_aln_line_marker(alignment_file_path,threads)
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
    for rname in gene_regions_read_length:
        for gname in gene_regions_read_length[rname]:
            for region in gene_regions_read_length[rname][gname]:
                read_lens += gene_regions_read_length[rname][gname][region]
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
        num_long_reads = 0
        gene_full_length_region_dict = defaultdict(lambda:{})
        isoform_max_num_exons_dict = {}
        for rname in gene_regions_dict:
            for gname in gene_regions_dict[rname]:
                regions_set = set()
                isoform_region_dict = defaultdict(lambda:set())
                for region in gene_regions_dict[rname][gname]:
                    for isoform in gene_regions_dict[rname][gname][region]:
                        isoform_region_dict[isoform].add(region)
                for isoform in isoform_region_dict:
                    max_region_exon_num = 0
                    longest_region = ''
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if max_region_exon_num < region_exon_num:
                            max_region_exon_num = region_exon_num
                            longest_region = region
                    isoform_max_num_exons_dict[isoform] = max_region_exon_num
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if region_exon_num == max_region_exon_num:
                            regions_set.add(region)
                gene_full_length_region_dict[rname][gname] = regions_set
        filtered_gene_regions_read_length = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:[])))
        for rname in gene_regions_read_length.copy():
            for gname in gene_regions_read_length[rname].copy():    
                # region_lens = []
                for region in gene_regions_read_length[rname][gname].copy():
                    region_exon_num = region.count(':')
                    read_length_list = []
                    for read_length in gene_regions_read_length[rname][gname][region]:
                        is_valid_read = True
                        if (filtering):
                            for isoform in gene_regions_dict[rname][gname][region]:
                                if (isoform_max_num_exons_dict[isoform] - region_exon_num  > 7) and (read_length / gene_isoforms_length_dict[rname][gname][isoform] <= 0.2):
                                        is_valid_read = False
                                        filtered_gene_regions_read_length[rname][gname][region].append(read_length)
                                        break
                        if (is_valid_read):
                            read_length_list.append(read_length)
                    gene_regions_read_length[rname][gname][region] = read_length_list
                    gene_regions_read_count[rname][gname][region] = len(read_length_list)
                    num_long_reads += len(read_length_list)

                    if gene_regions_read_count[rname][gname][region] == 0:
                        if region not in gene_full_length_region_dict[rname][gname]:
                            del gene_regions_read_length[rname][gname][region]
                            del gene_regions_read_count[rname][gname][region]
        return gene_regions_read_count,gene_regions_read_length,sum(read_lens),num_long_reads,filtered_gene_regions_read_length
    else:
        all_read_len = []
        for rname in gene_regions_read_count.copy():
            for gname in gene_regions_read_count[rname].copy():    
                # region_lens = []
                for region in gene_regions_read_count[rname][gname].copy():
                    if gene_regions_read_count[rname][gname][region] == 0:
                        if config.sr_region_selection == 'real_data':
                            if config.keep_sr_exon_region:
                                if check_region_type(region) != 'one_exon':
                                    del gene_regions_read_count[rname][gname][region]
                            else:
                                del gene_regions_read_count[rname][gname][region]
                    else:
                        all_read_len += gene_regions_read_length[rname][gname][region]
        SR_read_len = sum(all_read_len) / len(all_read_len)
        return gene_regions_read_count,SR_read_len,num_mapped_lines