from parse_annotation import parse_annotation
from construct_feature_matrix import filter_multi_exons_regions
import time
def parse_reference_annotation(ref_file_path,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    print('Parsing the reference annotation...')
    start = time.time()
    [_, gene_points_dict, gene_isoforms_dict, genes_regions_len_dict,
        _, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict] = parse_annotation(ref_file_path, READ_LEN, READ_JUNC_MIN_MAP_LEN)
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