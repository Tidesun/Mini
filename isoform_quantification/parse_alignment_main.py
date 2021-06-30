from parse_alignment import map_read,map_long_read,parse_read_line
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from functools import partial
import concurrent.futures
import numpy as np
import time
import random
# from memory_profiler import profile

# @profile
def parse_alignment(alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_points_dict,gene_range,gene_interval_tree_dict,gene_regions_dict,genes_regions_len_dict,long_read,threads):
    start_t = time.time()
    read_names = set()
    total_read_length = 0
    gene_regions_read_count = defaultdict(lambda: defaultdict(dict))
    gene_regions_read_length = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda : [])))
    filtered_gene_regions_dict = defaultdict(lambda: defaultdict(dict))
    for chr in gene_regions_dict:
        for gene in gene_regions_dict[chr]:
            # if (not long_read):
            #     per_gene_regions_dict = filter_regions(gene_regions_dict[chr][gene],long_read=False)
            # else:
            #     per_gene_regions_dict =  filter_regions(gene_regions_dict[chr][gene],long_read=True)
            per_gene_regions_dict =  gene_regions_dict[chr][gene]
            for region in per_gene_regions_dict:
                gene_regions_read_count[chr][gene][region] = 0
                filtered_gene_regions_dict[chr][gene][region] = True
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
    read_lens = []
    aln_lines = []
    # for each read in the alignment
    ## for testing
    # import pickle
    if (long_read):
        mapping_pkl_path = '/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/human_simulation/jobs/hybrid_simulation/validation/jobs/GM/lr_only/mapping.pkl'
    else:
        mapping_pkl_path = '/fs/project/PCON0009/Au-scratch2/haoran/quantification_evaluation/human_simulation/jobs/hybrid_simulation/validation/jobs/GM/sr_only/mapping.pkl'
    mapping_pkl_file =  open(mapping_pkl_path,'ab')
    mapping_pkl_file.truncate(0)
    MAX_BUFF_LINES_IN_MEM = 1e5
    with open(alignment_file_path, 'r') as aln_file:
        aln_lines = []
        for line in aln_file:
             # parse SAM files
            try:
                if line[0] == '@':
                    continue
                fields = line.split('\t')
                if (fields[1] == '4'):
                    continue
                [read_name, read_start_pos, rname, read_len_list] = parse_read_line(line, READ_LEN)
                aln_lines.append([read_name, read_start_pos, rname, read_len_list])
                read_names.add(fields[0])
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
            # if not exceed buffer
            if (len(aln_lines) < MAX_BUFF_LINES_IN_MEM):
                continue
            else: 
                print('Buffer exceeded.Starting mapping...')
                # if exceed buffer
                if (long_read):
                    map_f = map_long_read
                else:
                    map_f = map_read   
                partial_map_read = partial(map_f,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
                    start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                    READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)      
                if threads > 1:
                    # parallel
                    chunksize, extra = divmod(len(aln_lines), threads)
                    if extra:
                        chunksize += 1
                    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
                        all_mappings = [executor.map(partial_map_read, aln_lines,chunksize=chunksize)]
                        # all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
                else:
                    # serial
                    all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
                for mapping in all_mappings:
                    if (mapping['read_mapped']):
                        for mapping_area in [random.choice(mapping['mapping_area'])]:
                            rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                            if region_name in gene_regions_read_count[rname][gname]:
                                gene_regions_read_count[rname][gname][region_name] += 1
                        read_lens.append(mapping['read_length'])
                # for testing
                # pickle.dump(all_mappings,mapping_pkl_file)
                aln_lines = []
                print('Done.')
        if (len(aln_lines)>0):
            if (long_read):
                map_f = map_long_read
            else:
                map_f = map_read   
            partial_map_read = partial(map_f,gene_points_dict,gene_interval_tree_dict, filtered_gene_regions_dict,
                start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)     
            all_mappings = [partial_map_read(aln_line) for aln_line in aln_lines]
            for mapping in all_mappings:
                if (mapping['read_mapped']):
                    for mapping_area in [random.choice(mapping['mapping_area'])]:
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                        if region_name in gene_regions_read_count[rname][gname]:
                            gene_regions_read_count[rname][gname][region_name] += 1
                    read_lens.append(mapping['read_length'])
            # for testing
            # pickle.dump(all_mappings,mapping_pkl_file)
            aln_lines = []
    mapping_pkl_file.close()
        
    if (long_read):
        return gene_regions_read_count,gene_regions_read_length,sum(read_lens),len(read_names)
    else:
        SR_read_len = 150
        print('Short read length is {} bp'.format(SR_read_len))
        return gene_regions_read_count,SR_read_len,len(read_names)