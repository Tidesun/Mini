from parse_annotation import parse_annotation
from collections import defaultdict
def check_region_type(region_name):
    if ((region_name.count(':') == 1) and ('-' not in region_name)):
        return 'one_exon'
    elif ((region_name.count(':') == 2) and ('-' not in region_name)):
        return 'two_exons'
    elif (region_name.count('-') == 1):
        return 'one_junction'
    elif ((region_name.count(':') > 2) and ('-' not in region_name)):
        return 'exons'
    else:
        return 'junctions'
def calculate_exon_min_read_mapped_length(exon_region_name,point_dict,exon_position):
    assert '-' not in exon_region_name
    points = exon_region_name.split(':')
    if len(points) == 2:
        return 1
    if exon_position == 'left':
        return point_dict[points[-1]] - point_dict[points[1]] + 1 + 1
    elif exon_position == 'right':
        return point_dict[points[-2]] - point_dict[points[0]] + 1 + 1
    elif exon_position == 'center':
        return point_dict[points[-2]] - point_dict[points[1]] + 2 + 1       
# def filter_regions(gene_regions_dict,gene_points_dict,genes_regions_len_dict,READ_JUNC_MIN_MAP_LEN,min_read_len,max_read_len=None,LR_gene_read_min_len_dict=None):
#     new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
#     new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
#     for chr_name in gene_regions_dict:
#         for gene_name in gene_regions_dict[chr_name]:
#             if LR_gene_read_min_len_dict is not None:
#                 if (chr_name not in LR_gene_read_min_len_dict) or (gene_name not in LR_gene_read_min_len_dict[chr_name]):
#                     continue
#                 min_read_len = LR_gene_read_min_len_dict[chr_name][gene_name]
#             point_dict = {}
#             for coord in gene_points_dict[chr_name][gene_name]:
#                 point_dict['P{}'.format(gene_points_dict[chr_name][gene_name][coord])] = int(coord)
#             for region_name in gene_regions_dict[chr_name][gene_name]:
#                 is_region_valid = False
#                 if (genes_regions_len_dict[chr_name][gene_name][region_name] > min_read_len):              
#                     if (max_read_len is not None):
#                         if check_region_type(region_name) in ['one_exon','two_exons']:
#                             is_region_valid = True
#                         elif check_region_type(region_name) == 'one_junction':
#                             exons = region_name.split('-')
#                             for i,exon in zip(range(len(exons)),exons):
#                                 points = exon.split(':')
#                                 if (i == 0):
#                                     first_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
#                                 elif i == len(exons) - 1:
#                                     last_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
#                             if (first_exon_length > READ_JUNC_MIN_MAP_LEN) and (last_exon_length > READ_JUNC_MIN_MAP_LEN):
#                                 if region_name.count(':') == 2:
#                                     is_region_valid = True
#                                 else:
#                                     [exon_1,exon_2] = region_name.split('-')
#                                     if 2 * READ_JUNC_MIN_MAP_LEN - 2 + calculate_exon_min_read_mapped_length(exon_1,point_dict,exon_position='left') + calculate_exon_min_read_mapped_length(exon_2,point_dict,exon_position='right') <= max_read_len:
#                                         is_region_valid = True
#                         elif check_region_type(region_name) == 'exons':
#                             if calculate_exon_min_read_mapped_length(region_name,point_dict,exon_position='center') <= max_read_len:
#                                 is_region_valid = True
#                         elif check_region_type(region_name) == 'junctions':
#                             exons = region_name.split('-')
#                             exon_length = 0
#                             for i,exon in zip(range(len(exons)),exons):
#                                 points = exon.split(':')
#                                 if (i == 0):
#                                     first_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
#                                 elif i == len(exons) - 1:
#                                     last_exon_length = point_dict[points[-1]] - point_dict[points[0]] + 1
#                                 else:
#                                     exon_length += point_dict[points[-1]] - point_dict[points[0]] + 1
#                             if (first_exon_length > READ_JUNC_MIN_MAP_LEN) and (last_exon_length > READ_JUNC_MIN_MAP_LEN):
#                                 if exon_length + 2 * READ_JUNC_MIN_MAP_LEN - 2 + calculate_exon_min_read_mapped_length(exons[0],point_dict,exon_position='left') + calculate_exon_min_read_mapped_length(exons[-1],point_dict,exon_position='right') <= max_read_len:
#                                     is_region_valid = True
#                     else:
#                         is_region_valid = True
#                 if (is_region_valid):
#                     new_gene_regions_dict[chr_name][gene_name][region_name] = gene_regions_dict[chr_name][gene_name][region_name]
#                     new_genes_regions_len_dict[chr_name][gene_name][region_name] = genes_regions_len_dict[chr_name][gene_name][region_name]
#     return new_gene_regions_dict,new_genes_regions_len_dict
def filter_regions(gene_regions_dict,genes_regions_len_dict):
    new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
    new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
    for rname in gene_regions_dict:
        for gname in gene_regions_dict[rname]:
            regions_set = set()
            for region in gene_regions_dict[rname][gname]:
                if check_region_type(region) in ['one_exon','one_junction']:
                    regions_set.add(region)
            for region_name in regions_set:
                new_gene_regions_dict[rname][gname][region_name] = gene_regions_dict[rname][gname][region_name]
                new_genes_regions_len_dict[rname][gname][region_name] = genes_regions_len_dict[rname][gname][region_name]
    return new_gene_regions_dict,new_genes_regions_len_dict
# def filter_long_read_regions(gene_regions_dict,genes_regions_len_dict):
#     new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
#     new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
#     for rname in gene_regions_dict:
#         for gname in gene_regions_dict[rname]:
#             regions_set = set()
#             isoform_region_dict = defaultdict(lambda:set())
#             for region in gene_regions_dict[rname][gname]:
#                 for isoform in gene_regions_dict[rname][gname][region]:
#                     isoform_region_dict[isoform].add(region)
#             for isoform in isoform_region_dict:
#                 max_region_exon_num = 0
#                 longest_region = ''
#                 allowed_regions_set = set()
#                 for region in isoform_region_dict[isoform]:
#                     region_exon_num = region.count(':')
#                     if max_region_exon_num < region_exon_num:
#                         max_region_exon_num = region_exon_num
#                         longest_region = region
#                     if region_exon_num >= 3:
#                         allowed_regions_set.add(region)
#                 if max_region_exon_num > 3:
#                     regions_set = regions_set.union(allowed_regions_set)
#                 else:
#                     if longest_region != '':
#                         regions_set.add(longest_region)
#             for region_name in regions_set:
#                 new_gene_regions_dict[rname][gname][region_name] = gene_regions_dict[rname][gname][region_name]
#                 new_genes_regions_len_dict[rname][gname][region_name] = genes_regions_len_dict[rname][gname][region_name]
#     return new_gene_regions_dict,new_genes_regions_len_dict

def filter_long_read_regions(gene_regions_dict,genes_regions_len_dict):
    new_gene_regions_dict = defaultdict(lambda:defaultdict(dict))
    new_genes_regions_len_dict = defaultdict(lambda:defaultdict(dict))
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
                if max_region_exon_num <= 1:
                    if longest_region != '':
                        regions_set.add(longest_region)
                else:
                    for region in isoform_region_dict[isoform]:
                        region_exon_num = region.count(':')
                        if (region_exon_num >= max_region_exon_num - 1):
                            regions_set.add(region)
            for region_name in regions_set:
                new_gene_regions_dict[rname][gname][region_name] = gene_regions_dict[rname][gname][region_name]
                new_genes_regions_len_dict[rname][gname][region_name] = genes_regions_len_dict[rname][gname][region_name]
    return new_gene_regions_dict,new_genes_regions_len_dict

                                
                        


def parse_reference_annotation(ref_file_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN,LR_gene_read_min_len_dict):
    [gene_exons_dict, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        _, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict] = parse_annotation(ref_file_path, threads,READ_LEN, READ_JUNC_MIN_MAP_LEN)
    # SR_gene_regions_dict,SR_genes_regions_len_dict = filter_regions(gene_regions_dict,gene_points_dict,genes_regions_len_dict,READ_JUNC_MIN_MAP_LEN,150,150,None)
    SR_gene_regions_dict,SR_genes_regions_len_dict = filter_regions(gene_regions_dict,genes_regions_len_dict)
    LR_gene_regions_dict,LR_genes_regions_len_dict = gene_regions_dict,genes_regions_len_dict
    # LR_gene_regions_dict,LR_genes_regions_len_dict = filter_long_read_regions(gene_regions_dict,genes_regions_len_dict)
    for chr_name in gene_points_dict:
        for gene_name in gene_points_dict[chr_name].copy():
            if (len(SR_gene_regions_dict[chr_name][gene_name]) == 0 or len(LR_gene_regions_dict[chr_name][gene_name]) == 0):
                for dic in [gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict]:
                    if chr_name in dic and gene_name in dic[chr_name]:
                        del dic[chr_name][gene_name]
    return gene_exons_dict,gene_points_dict,gene_isoforms_dict,SR_gene_regions_dict,SR_genes_regions_len_dict,LR_gene_regions_dict,LR_genes_regions_len_dict,gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict
from intervaltree import IntervalTree
def process_annotation_for_alignment(gene_exons_dict,gene_points_dict):
    print('Process the reference annotation for read_count...')
    # process the reference annotation for read_count TODO: a universal format
    gene_regions_points_list,gene_range,gene_interval_tree_dict = dict(),dict(),dict()
    for chr_name in gene_points_dict:
        gene_regions_points_list[chr_name],gene_range[chr_name],gene_interval_tree_dict[chr_name] = dict(),list(),dict()
        for gene_name in gene_points_dict[chr_name]:
            points = list(gene_points_dict[chr_name][gene_name].keys())
            gene_regions_points_list[chr_name][gene_name] = points
            gene_range[chr_name].append([gene_name,points[0],points[-1]])
            gene_interval_tree_dict[chr_name][gene_name] = IntervalTree()
            for [start_pos,end_pos,_] in gene_exons_dict[chr_name][gene_name]:
                # Interval tree exclude end position
                gene_interval_tree_dict[chr_name][gene_name].addi(start_pos, end_pos+1, None)
    return gene_regions_points_list,gene_range,gene_interval_tree_dict

