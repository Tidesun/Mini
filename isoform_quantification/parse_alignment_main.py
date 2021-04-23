from parse_alignment import map_read,map_long_read,parse_read_line
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
from construct_feature_matrix import filter_regions
from functools import partial
import concurrent.futures
import numpy as np
import time
def parse_alignment(alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict,genes_regions_len_dict,long_read,threads):
    start_t = time.time()
    read_names = set()
    total_read_length = 0
    gene_regions_read_count = defaultdict(lambda: defaultdict(dict))
    gene_regions_read_length = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda : [])))
    for chr in gene_regions_dict:
        for gene in gene_regions_dict[chr]:
            if (not long_read):
                per_gene_regions_dict = filter_regions(gene_regions_dict[chr][gene],long_read=False)
            else:
                per_gene_regions_dict =  filter_regions(gene_regions_dict[chr][gene],long_read=True)
            for region in per_gene_regions_dict:
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
    read_lens = []
    aln_lines = []
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
                    [read_name, read_start_pos, rname, read_len_list] = parse_read_line(line, READ_LEN)
                    aln_lines.append([read_name, read_start_pos, rname, read_len_list])
                    read_names.add(fields[0])
            except Exception as e:
                tb = traceback.format_exc()
                print(Exception('Failed to on ' + line, tb))
                continue
    if (long_read):
        for parsed_line in aln_lines:
            mapping,gene_regions_read_count,gene_regions_read_length = map_long_read(parsed_line, gene_regions_read_count,gene_regions_read_length, gene_regions_points_list, 
            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
            if ((mapping['read_mapped']) & (fields[0] not in read_names)):
                total_read_length += sum(mapping['read_len'])
    else:
        # Determine the chunksize
        chunksize, extra = divmod(len(aln_lines), threads)
        if extra:
            chunksize += 1
        print(len(aln_lines))
        print(chunksize)
        print(threads)
        # for parsed_line in aln_lines:
        #     mapping = map_read(parsed_line,*args)
        #     if (mapping['read_mapped']):
        #         mapping_area = mapping['mapping_area']
        #         rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
        #         if region_name in gene_regions_read_count[rname][gname]:
        #             gene_regions_read_count[rname][gname][region_name] += 1
        #             read_lens.append(mapping['read_length'])
        # for parsed_line in aln_lines:
        #     mapping = helper_map_read(parsed_line)
        #     if (mapping['read_mapped']):
        #         mapping_area = mapping['mapping_area']
        #         rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
        #         if region_name in gene_regions_read_count[rname][gname]:
        #             gene_regions_read_count[rname][gname][region_name] += 1
        #             read_lens.append(mapping['read_length'])
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            partial_map_read = partial(map_read,gene_regions_points_list, genes_regions_len_dict,
            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
            for mapping in executor.map(partial_map_read, aln_lines,chunksize=chunksize):
                try:
                    if (mapping['read_mapped']):
                        mapping_area = mapping['mapping_area']
                        rname,gname,region_name = mapping_area['chr_name'],mapping_area['gene_name'],mapping_area['region_name']
                        if region_name in gene_regions_read_count[rname][gname]:
                            gene_regions_read_count[rname][gname][region_name] += 1
                            read_lens.append(mapping['read_length'])
                except Exception as e:
                    print(e)
    if (long_read):
        return gene_regions_read_count,gene_regions_read_length,total_read_length,len(read_names)
    else:
        SR_read_len = np.median(read_lens)
        print('Short read length is {} bp'.format(SR_read_len))
        return gene_regions_read_count,SR_read_len,len(read_names)