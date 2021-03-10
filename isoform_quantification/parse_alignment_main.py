from parse_alignment import map_read,map_long_read
from collections import defaultdict
import traceback
from operator import itemgetter, attrgetter
def parse_alignment(alignment_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN,gene_regions_points_list,gene_range,gene_regions_dict,long_read):
    read_names = set()
    total_read_length = 0
    gene_regions_read_count = defaultdict(lambda: defaultdict(dict))
    gene_regions_read_length = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda : [])))
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
                        mapping,gene_regions_read_count = map_read(line, gene_regions_read_count, gene_regions_points_list, 
                            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
                    else:
                        mapping,gene_regions_read_count,gene_regions_read_length = map_long_read(line, gene_regions_read_count,gene_regions_read_length, gene_regions_points_list, 
                            start_pos_list, start_gname_list, end_pos_list, end_gname_list,
                            READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST)
                        if ((mapping['read_mapped']) & (fields[0] not in read_names)):
                            total_read_length += sum(mapping['read_len'])

                    read_names.add(fields[0])
            except Exception as e:
                tb = traceback.format_exc()
                continue
                raise Exception('Failed to on ' + line, tb)
    if (long_read):
        return gene_regions_read_count,gene_regions_read_length,total_read_length,len(read_names)
    else:
        return gene_regions_read_count,len(read_names)