import pandas as pd
import numpy as np
def get_reads_isoform_info(isoform_len_dict,isoform_exon_dict,strand_dict,gene_regions_read_mapping,LR_gene_regions_dict):
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
    multi_mapping_reads_isoform_info = {}
    for read_name in reads_isoform_info:
        if len(reads_isoform_info[read_name]) > 1:
            multi_mapping_reads_isoform_info[read_name] = reads_isoform_info[read_name]
    return isoform_len_dict,reads_isoform_info,multi_mapping_reads_isoform_info,expression_dict,unique_mapping_expression_dict
def get_read_len_dist(reads_isoform_info,isoform_len_dict):
    read_len_dict = {}
    read_len_dist_dict = {}
    for read in reads_isoform_info:
        [isoform_len,f_end,t_end] = next(iter(reads_isoform_info[read].values()))
        read_len = isoform_len -  f_end - t_end
        read_len_dict[read] = read_len
        if read_len not in read_len_dist_dict:
            read_len_dist_dict[read_len] = 0
        read_len_dist_dict[read_len] += 1
    read_len_dist = pd.Series(read_len_dist_dict).to_frame()
    read_len_dist.columns = ['PDF']
    read_len_dist = read_len_dist.sort_index(ascending=True)
    isoform_len_df = pd.Series(isoform_len_dict)
    return read_len_dict,read_len_dist,isoform_len_df