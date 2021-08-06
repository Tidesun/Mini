import re
from collections import defaultdict
import numpy as np
import pandas as pd
def sync_reference_name(ref_name):
    ref_name = ref_name.upper()
    match = re.search("(?<=CHR).*", ref_name)
    if match:
        ref_name = match.group(0)
    if ref_name == "M":
        ref_name = "MT"
    if '_' in ref_name:
        ref_name = ''
    return ref_name
def get_very_short_isoforms(output_path,filtered_gene_regions_read_length,LR_gene_regions_dict,isoform_length_dict):
    short_isoform_genes = set()
    short_isoforms = set()
    mapped_genes = set()
    filtered_out_read_lengths = {}
    for rname in filtered_gene_regions_read_length:
        for gname in filtered_gene_regions_read_length[rname]:
            read_lengths = []
            isoform_list = set()
            short_isoforms_set = set()
            filtered_out_read_lengths[gname] = []
            for region in LR_gene_regions_dict[rname][gname]:
                read_lengths += filtered_gene_regions_read_length[rname][gname][region]
                for isoform in LR_gene_regions_dict[rname][gname][region]:
                    isoform_list.add(isoform)
            if len(read_lengths) > 0:       
                for isoform in isoform_list:
                    if isoform_length_dict[isoform]<max(read_lengths):
                        short_isoforms_set.add(isoform)
                        short_isoforms.add(isoform)
                if len(short_isoforms_set) >= 1:
                    short_isoform_genes.add(gname)
                mapped_genes.add(gname)
            filtered_out_read_lengths[gname] += read_lengths
    with open('{}/filtered_out_read_lengths.tsv'.format(output_path),'w') as f:
        for gname in filtered_out_read_lengths:
            f.write('{}\n'.format(gname))
            if len(filtered_out_read_lengths[gname]) == 0:
                f.write('Length threshold:\n')
            else:
                f.write('Length threshold:{}\n'.format(max(filtered_out_read_lengths[gname])))
            for read_length in filtered_out_read_lengths[gname]:
                f.write('{}\n'.format(read_length))
    with open('{}/short_isoform.tsv'.format(output_path),'w') as f:
        for isoform in short_isoforms:
            f.write('{}\n'.format(isoform))
    with open('{}/short_isoform_genes.tsv'.format(output_path),'w') as f:
        for gene in short_isoform_genes:
            f.write('{}\n'.format(gene))
    with open('{}/all_mapped_genes.tsv'.format(output_path),'w') as f:
        for gene in mapped_genes:
            f.write('{}\n'.format(gene))
def get_filtered_out_long_read_M_dist(output_path,filtered_gene_regions_read_length,LR_gene_regions_dict):
    filtered_LR_gene_regions_dict = LR_gene_regions_dict.copy()
    delta_s_read_lengths = defaultdict(lambda:[])
    for rname in filtered_LR_gene_regions_dict:
        for gname in filtered_LR_gene_regions_dict[rname]:
            for region in filtered_LR_gene_regions_dict[rname][gname].copy():
                if len(filtered_gene_regions_read_length[rname][gname][region]) == 0:
                    del filtered_LR_gene_regions_dict[rname][gname][region]
    isoform_region_dict = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:set())))
    for rname in filtered_LR_gene_regions_dict:
        for gname in filtered_LR_gene_regions_dict[rname]:
            for region in filtered_LR_gene_regions_dict[rname][gname]:
                for isoform in filtered_LR_gene_regions_dict[rname][gname][region]:
                    isoform_region_dict[rname][gname][isoform].add(region)
    isoform_region_max_length = {}
    for rname in isoform_region_dict:
        for gname in isoform_region_dict[rname]:
            for isoform in isoform_region_dict[rname][gname]:
                max_length = 0
                max_region = ''
                for region in isoform_region_dict[rname][gname][isoform]:
                    region_length = region.count(':')
                    if max_length < region_length:
                        max_length = region_length
                        max_region = region
                isoform_region_max_length[isoform] = (max_length,max_region)
    lr_unique_dist = defaultdict(lambda:0)
    lr_multi_dist = defaultdict(lambda:0)
    for rname in filtered_gene_regions_read_length:
        for gname in filtered_gene_regions_read_length[rname]:
            for region in filtered_gene_regions_read_length[rname][gname]:
                if region in filtered_LR_gene_regions_dict[rname][gname]:
                    region_count = len(filtered_gene_regions_read_length[rname][gname][region])
                    num_exons = region.count(':')
                    if len(filtered_LR_gene_regions_dict[rname][gname][region]) == 1:
                        isoform = list(filtered_LR_gene_regions_dict[rname][gname][region])[0]
                        M = isoform_region_max_length[isoform][0] - num_exons
                        lr_unique_dist[M] += region_count
                        if M == 0:
                            for read_length in filtered_gene_regions_read_length[rname][gname][region]:
                                delta_s_read_lengths[gname].append((read_length,isoform,region,isoform_region_max_length[isoform][0]))
                    else:
                        M_list = []
                        for isoform in filtered_LR_gene_regions_dict[rname][gname][region]:
                            M = isoform_region_max_length[isoform][0] - num_exons
                            M_list.append(M)
                        lr_multi_dist[np.median(M_list)] += region_count
                    
    unique_dist_df = pd.DataFrame({'count':pd.Series(lr_unique_dist)}).reset_index()
    multi_dist_df = pd.DataFrame({'count':pd.Series(lr_multi_dist)}).reset_index()
    with open('{}/delta_s_filtered_out_read_lengths.tsv'.format(output_path),'w') as f:
        for gname in delta_s_read_lengths:
            if len(delta_s_read_lengths[gname]) > 0:
                f.write('{}\n'.format(gname))
                for read_length,isoform,region,isoform_region_max_length in delta_s_read_lengths[gname]:
                    f.write('{}\t{}\t{}\t{}\n'.format(read_length,isoform,region,isoform_region_max_length))
    return unique_dist_df,multi_dist_df
def get_long_read_M_dist(long_read_gene_regions_read_count,LR_gene_regions_dict):
    filtered_LR_gene_regions_dict = LR_gene_regions_dict.copy()
    for rname in filtered_LR_gene_regions_dict:
        for gname in filtered_LR_gene_regions_dict[rname]:
            for region in filtered_LR_gene_regions_dict[rname][gname].copy():
                if region not in long_read_gene_regions_read_count[rname][gname]:
                    del filtered_LR_gene_regions_dict[rname][gname][region]
    isoform_region_dict = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:set())))
    for rname in filtered_LR_gene_regions_dict:
        for gname in filtered_LR_gene_regions_dict[rname]:
            for region in filtered_LR_gene_regions_dict[rname][gname]:
                for isoform in filtered_LR_gene_regions_dict[rname][gname][region]:
                    isoform_region_dict[rname][gname][isoform].add(region)
    isoform_region_max_length = {}
    for rname in isoform_region_dict:
        for gname in isoform_region_dict[rname]:
            for isoform in isoform_region_dict[rname][gname]:
                max_length = 0
                max_region = ''
                for region in isoform_region_dict[rname][gname][isoform]:
                    region_length = region.count(':')
                    if max_length < region_length:
                        max_length = region_length
                        max_region = region
                isoform_region_max_length[isoform] = (max_length,max_region)
    lr_unique_dist = defaultdict(lambda:0)
    lr_multi_dist = defaultdict(lambda:0)
    for rname in long_read_gene_regions_read_count:
        for gname in long_read_gene_regions_read_count[rname]:
            for region in long_read_gene_regions_read_count[rname][gname]:
                if long_read_gene_regions_read_count[rname][gname][region] > 0:
                    num_exons = region.count(':')
                    if len(filtered_LR_gene_regions_dict[rname][gname][region]) == 1:
                        isoform = list(filtered_LR_gene_regions_dict[rname][gname][region])[0]
                        M = isoform_region_max_length[isoform][0] - num_exons
                        # if  M == 87:
                        #     print(isoform)
                        #     print(region)
                        #     print(isoform_region_dict[rname][gname][isoform])
                        lr_unique_dist[M] += long_read_gene_regions_read_count[rname][gname][region]
                    else:
                        M_list = []
                        for isoform in filtered_LR_gene_regions_dict[rname][gname][region]:
                            M = isoform_region_max_length[isoform][0] - num_exons
                            M_list.append(M)
                        lr_multi_dist[np.median(M_list)] += long_read_gene_regions_read_count[rname][gname][region]
    unique_dist_df = pd.DataFrame({'count':pd.Series(lr_unique_dist)}).reset_index()
    multi_dist_df = pd.DataFrame({'count':pd.Series(lr_multi_dist)}).reset_index()
    return unique_dist_df,multi_dist_df
    