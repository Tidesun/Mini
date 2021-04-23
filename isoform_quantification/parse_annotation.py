#!/usr/bin/python
from operator import itemgetter, attrgetter
import re
import sys
import time
import numpy as np
import concurrent.futures

"""
Note: 
"""

### Compute number of bases to junction gap from an exon index (exclusive)
### exon_list format:  [start+point, end_point, len]
##########
def comp_len_to_junc_gap_forward(start_index, exon_list, isoform_exons, isoform_points_dict):
    if (start_index >= len(exon_list)):
        print('Warning: Called comp_len_to_junc_gap_forward with start_idx of out of bound') 
        return 0
    index = start_index
    l = 0
    next_index = index +1
    while (next_index< len(exon_list)):
        region_name_temp = 'P' + str(isoform_points_dict[exon_list[next_index][0]]) + ':P' + str(isoform_points_dict[exon_list[next_index][1]])
        if (region_name_temp in isoform_exons):   # this exon is in the isoform
            if ((exon_list[next_index][0] - 1) == exon_list[index][0]):
                l += exon_list[next_index][2]
                index = next_index
            else:
                break
        next_index += 1
    return l
##########
def comp_len_to_junc_gap_backward(start_index, exon_list, isoform_exons, isoform_points_dict):
    if (start_index >= len(exon_list)):
        print('Warning: Called comp_len_to_junc_gap_backward with start_idx of out of bound') 
        return 0
    index = start_index
    l = 0
    next_index = index -1
    while (next_index >= 0):
        region_name_temp = 'P' + str(isoform_points_dict[exon_list[next_index][0]]) + ':P' + str(isoform_points_dict[exon_list[next_index][1]])
        if (region_name_temp in isoform_exons):   # this exon is in the isoform
            if ((exon_list[next_index][1] + 1) == exon_list[index][0]):
                l += exon_list[next_index][2]
                index = next_index
            else:
                break
        next_index -= 1
    return l
##########
def split_and_sort_exons(gene_exons_dict):
    new_gene_exons_dict = {}
    for chr_name in list(gene_exons_dict.keys()):
        new_gene_exons_dict[chr_name] = {}
        for gene_name in list(gene_exons_dict[chr_name].keys()):
            exon_sorted = sorted(gene_exons_dict[chr_name][gene_name], key=itemgetter(0, 1))  # Sort by start position then end position
            i = 1
            while (i < len(exon_sorted)):
                if ((exon_sorted[i][0] == exon_sorted[i - 1][0]) and 
                    (exon_sorted[i][1] == exon_sorted[i - 1][1])):
                    del exon_sorted[i]  # Delete the repeated exon
                elif (exon_sorted[i][0] <= exon_sorted[i - 1][1]):
                    temp_exons = sorted([exon_sorted[i][0], exon_sorted[i - 1][0], exon_sorted[i][1], exon_sorted[i - 1][1]])
                    del exon_sorted[i - 1]  # Delete the two overlapping exons
                    del exon_sorted[i - 1]
                    if (temp_exons[0] == temp_exons[1]):  # Based on two exons overlap type, re-generate 2 or 3 new ones
                        exon_sorted.insert(i - 1, [temp_exons[1], temp_exons[2]])
                        exon_sorted.insert(i , [temp_exons[2] + 1, temp_exons[3]])
                    elif (temp_exons[2] == temp_exons[3]):
                        exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        exon_sorted.insert(i , [temp_exons[1], temp_exons[2]])
                    else:
                        exon_sorted.insert(i - 1, [temp_exons[0], temp_exons[1] - 1])
                        exon_sorted.insert(i, [temp_exons[1], temp_exons[2]])
                        exon_sorted.insert(i + 1, [temp_exons[2] + 1, temp_exons[3]])
                    exon_sorted = sorted(exon_sorted, key=itemgetter(0, 1))  # re-sort the exon positions
                else:
                    i += 1  # After sorting re-evaluate the same index unless there is no change
            new_gene_exons_dict[chr_name][gene_name] = [[start_pos,end_pos,end_pos - start_pos + 1] for [start_pos,end_pos] in exon_sorted if start_pos != end_pos]
    return new_gene_exons_dict
def cal_region_length(region_name,point_coord_dict):
    assert '-' in region_name or ':' in region_name
    assert 'P' in list(point_coord_dict.keys())[0]
    if '-' not in region_name:
        exons = [region_name]
    else:
        exons = [region for region in region_name.split('-') if ':' in region]
    length = 0
    for exon in exons:
        points = exon.split(':')
        length += point_coord_dict[points[-1]] - point_coord_dict[points[0]] + 1
    return length
def generate_multi_exons(single_gene_gene_points_dict,single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict):
    def concat_region_name(r1,r2):
        if '-' in r2:
            return r1 + r2.replace(r2.split('-')[0],'')
        else:
            return r1 + r2.replace(r2.split(':')[0],'')
    def is_neighbour_region(r1,r2):
        if '-' in r2:
            return r1.split('-')[-1] == r2.split('-')[0]
        else:
            r1_points = r1.replace(':','-').split('-')
            r2_points = r2.replace(':','-').split('-')
            return r1_points[-1] == r2_points[0]
    point_coord_dict = {}
    for p in single_gene_gene_points_dict:
        point_coord_dict['P{}'.format(single_gene_gene_points_dict[p])] = int(p)
    isoform_region_dict = {}
    isoform_region_len_dict = {}
    for region_name in single_gene_gene_regions_dict:
        for isoform_name in single_gene_gene_regions_dict[region_name]:
            if isoform_name in isoform_region_dict:
                isoform_region_dict[isoform_name].add(region_name)
                isoform_region_len_dict[isoform_name][region_name] = cal_region_length(region_name,point_coord_dict)
            else:
                isoform_region_dict[isoform_name] = set([region_name])
                isoform_region_len_dict[isoform_name] = {region_name:cal_region_length(region_name,point_coord_dict)}

    multi_exon_region_dict = isoform_region_dict.copy()
    multi_exon_region_len_dict = isoform_region_len_dict.copy()
    for isoform_name in isoform_region_dict:
        new_region_names = []
        new_region_len_dict = {}
        for region_name in isoform_region_dict[isoform_name]:
            new_region_name = region_name
            new_region_len = isoform_region_len_dict[isoform_name][region_name]
            region_candidates = set(isoform_region_dict[isoform_name])
            region_candidates.remove(region_name)
            while len(region_candidates) !=0:
                neighbour_region = ''
                for region_candidate in region_candidates:
                    if is_neighbour_region(new_region_name,region_candidate):
                        neighbour_region = region_candidate
                        break
                if neighbour_region != '':
                    region_candidates.remove(neighbour_region)
                    new_region_name = concat_region_name(new_region_name,neighbour_region)
                    if (new_region_name not in single_gene_gene_regions_dict):
                        new_region_names.append(new_region_name)
                        new_region_len += isoform_region_len_dict[isoform_name][neighbour_region]
                        new_region_len_dict[new_region_name] = new_region_len
                else:
                    break
        multi_exon_region_dict[isoform_name] = list(multi_exon_region_dict[isoform_name]) + new_region_names
        multi_exon_region_len_dict[isoform_name].update(new_region_len_dict)
    new_single_gene_isoforms_regions_len_dict = {}
    new_single_gene_gene_regions_dict = {}
    for isoform_name in multi_exon_region_dict:
        for region_name in multi_exon_region_dict[isoform_name]:
            if (multi_exon_region_len_dict[isoform_name][region_name] == 1):
                continue
            if region_name not in new_single_gene_gene_regions_dict:
                new_single_gene_gene_regions_dict[region_name] = set([isoform_name])
                new_single_gene_isoforms_regions_len_dict[region_name] = {isoform_name:multi_exon_region_len_dict[isoform_name][region_name]}
            else:
                new_single_gene_gene_regions_dict[region_name].add(isoform_name)
                new_single_gene_isoforms_regions_len_dict[region_name][isoform_name] = multi_exon_region_len_dict[isoform_name][region_name]
    return new_single_gene_isoforms_regions_len_dict,new_single_gene_gene_regions_dict



##########
def generate_exon_indicator_for_isoform_single_gene(single_gene_raw_isoform_exons_dict,exon_list,single_gene_gene_points_dict,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict = {},{}
    for isoform_name in single_gene_raw_isoform_exons_dict:
        isoform_exons = []
        #Question here?
#                 start_pos = [i-1 for i in single_gene_raw_isoform_exons_dict[isoform_name]['start_pos']]
        start_pos = single_gene_raw_isoform_exons_dict[isoform_name]['start_pos']
        end_pos = single_gene_raw_isoform_exons_dict[isoform_name]['end_pos']
        num_exons = len(single_gene_raw_isoform_exons_dict[isoform_name]['start_pos'])
        # Check the exon regions
        j = 0
        for i in range(len(exon_list)):
            flag = (j < num_exons)
            while (flag):
                p0 = exon_list[i][0]
                p1 = exon_list[i][1]
                l = exon_list[i][2]
                if ((p0 >= start_pos[j]) and
                    (p1 <= end_pos[j])):
                    region_name = 'P' + str(single_gene_gene_points_dict[p0]) + ':' + 'P' + str(single_gene_gene_points_dict[p1])
                    if (l >= READ_LEN):  # A valid region for exon
                        temp_isoform_name = set()
                        temp_isoform_name.add(isoform_name)
                        if (region_name in single_gene_gene_regions_dict):
                            temp_isoform_name = temp_isoform_name.union(single_gene_gene_regions_dict[region_name])
                        else:
                            single_gene_isoforms_regions_len_dict[region_name] = dict()
                        single_gene_isoforms_regions_len_dict[region_name][isoform_name] = l - READ_LEN + 1
                        single_gene_gene_regions_dict[region_name] = temp_isoform_name
                    isoform_exons.append(region_name)
                    flag = False
                elif (p0 == start_pos[j]):
                    print('Out-of-order exon position for isoform ' + isoform_name)
                    return {},{}
                elif (p0 < start_pos[j]):
                    flag = False
                else:
                    j += 1
                    flag = (j < num_exons)  
        # Check the junction regions
        for i in range(len(exon_list) - 1):     # the last exon can not have any junction
            p0 = exon_list[i][0]
            p1 = exon_list[i][1] + 1 # the end_pos is not included
            l  = exon_list[i][2]
            region_name_temp = 'P' + str(single_gene_gene_points_dict[p0]) + ':P' + str(single_gene_gene_points_dict[p1-1])
            if (region_name_temp not in isoform_exons):   # this exon is not in the isoform
                continue
            # compute start and end point of a region that starts with p1 point
            start = max(p1 - READ_LEN + 1, p0)
            end = p1 - READ_JUNC_MIN_MAP_LEN + min(comp_len_to_junc_gap_forward(i, exon_list, isoform_exons, single_gene_gene_points_dict),
                                                READ_JUNC_MIN_MAP_LEN - 1)
            if (end < start):
                continue    # exon is not long enough to be mapped by a read
            current_point = start
            while (current_point <= end):
                region_name = ""
                if ((current_point == p0) and
                    (p0 != (p1-1))):   # Special case when the exon length is 1 we dont want to repeat the point number
                    region_name += ('P' + str(single_gene_gene_points_dict[p0]) + '-')
                region_name += ('P' + str(single_gene_gene_points_dict[p1 - 1]) + '-')
                remain_len = READ_LEN - (p1 - current_point)
                j = i
                while (j < (len(exon_list) -1)):
                    j += 1
                    p0_temp = exon_list[j][0]
                    p1_temp = exon_list[j][1] + 1 # the end_pos is not included
                    l_temp  = exon_list[j][2]
                    region_name_temp = 'P' + str(single_gene_gene_points_dict[p0_temp]) + ':P' + str(single_gene_gene_points_dict[p1_temp-1])
                    if (region_name_temp not in isoform_exons):   # this exon is not in the isoform
                        continue
                    if (l_temp >= remain_len):
                        if ((comp_len_to_junc_gap_backward(j, exon_list, isoform_exons, single_gene_gene_points_dict) + remain_len ) >= READ_JUNC_MIN_MAP_LEN):
                            region_name += ('P' + str(single_gene_gene_points_dict[p0_temp]))
                            if ((l_temp == remain_len) and
                                (p0_temp != (p1_temp -1))):   # Special case when the exon length is 1 we dont want to repeat the point number
                                region_name += ('-P' + str(single_gene_gene_points_dict[p1_temp - 1]))
                            temp_isoform_name = set()
                            temp_isoform_name.add(isoform_name)
                            #temp_isoform_name = {isoform_name}                            
                            if (region_name in single_gene_gene_regions_dict):
                                temp_isoform_name = temp_isoform_name.union(single_gene_gene_regions_dict[region_name])
                            else:
                                single_gene_isoforms_regions_len_dict[region_name] = dict()
                            single_gene_gene_regions_dict[region_name] = temp_isoform_name
                            if (isoform_name in single_gene_isoforms_regions_len_dict[region_name]):
                                single_gene_isoforms_regions_len_dict[region_name][isoform_name] += 1
                            else:
                                single_gene_isoforms_regions_len_dict[region_name][isoform_name] = 1
                        break   # Found possible region for this current_point
                    else:
                        remain_len -= l_temp
                        region_name += ('P' + str(single_gene_gene_points_dict[p0_temp]) + '-')
                        if (p0_temp != (p1_temp -1)):   # Special case when the exon length is 1 we dont want to repeat the point number
                            region_name += ('P' + str(single_gene_gene_points_dict[p1_temp - 1]) + '-')
                current_point += 1 
    new_single_gene_isoforms_regions_len_dict,new_single_gene_gene_regions_dict = {},{}
    junctions,exons = [],[]
    for region_name in single_gene_gene_regions_dict:
        if '-' in region_name:
            if ('P0' not in region_name and 'P{}'.format(len(single_gene_gene_points_dict)-1) not in region_name):
                junctions.append(region_name)
        else:
            exons.append(region_name)
    # junctions = []
    # for junc1 in temp_junctions:
    #     is_valid_junc = True
    #     for junc2 in temp_junctions:
    #         if ((junc1 != junc2) and (junc1 in junc2)):
    #             is_valid_junc = False
    #             break
    #     if (is_valid_junc):
    #         junctions.append(junc1)

    for region_name in junctions:
        points = region_name.split('-')
        new_region_name = 'P{}:{}-{}:P{}'.format(int(points[0][1:])-1,points[0],points[-1],int(points[-1][1:])+1)
        new_single_gene_isoforms_regions_len_dict[new_region_name] = single_gene_isoforms_regions_len_dict[region_name].copy()
        new_single_gene_gene_regions_dict[new_region_name] = single_gene_gene_regions_dict[region_name].copy()
    for region_name in exons:
        new_single_gene_isoforms_regions_len_dict[region_name] = single_gene_isoforms_regions_len_dict[region_name].copy()
        new_single_gene_gene_regions_dict[region_name] = single_gene_gene_regions_dict[region_name].copy()
    if isoform_name == 'ENST00000484235':
        print(new_single_gene_gene_regions_dict)
    return new_single_gene_isoforms_regions_len_dict,new_single_gene_gene_regions_dict
def generate_exon_indicator_for_isoform_single_chr(raw_isoform_exons_dict,gene_exons_dict,gene_points_dict,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = {},{},{}
    for gene_name in raw_isoform_exons_dict:
        single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict = generate_exon_indicator_for_isoform_single_gene(raw_isoform_exons_dict[gene_name],gene_exons_dict[gene_name],gene_points_dict[gene_name],READ_LEN,READ_JUNC_MIN_MAP_LEN)
        isoforms_regions_len_dict[gene_name],gene_regions_dict[gene_name] = generate_multi_exons(gene_points_dict[gene_name],single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict)
        genes_regions_len_dict[gene_name] = {}
        point_coord_dict = {}
        for p in gene_points_dict[gene_name]:
            point_coord_dict['P{}'.format(gene_points_dict[gene_name][p])] = int(p)
        for region in gene_regions_dict[gene_name]:
            genes_regions_len_dict[gene_name][region] = cal_region_length(region,point_coord_dict)
    return isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict
def generate_exon_indicator_for_isoform(gene_exons_dict,gene_points_dict,raw_isoform_exons_dict,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        isoforms_regions_len_dict[chr_name],gene_regions_dict[chr_name],genes_regions_len_dict[chr_name] = {},{},{}

    # list_of_all_genes_chrs = [(gene_name,chr_name) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
    # list_of_args = [(raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name]) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
    # chunksize, extra = divmod(len(list_of_all_genes_chrs), 12)
    # if extra:
    #     chunksize += 1
    # print(chunksize)
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, executor.map(lambda args:generate_exon_indicator_for_isoform_single_gene(*args,READ_LEN,READ_JUNC_MIN_MAP_LEN),list_of_args,chunksize=chunksize)):
    #         isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name] = result
    
    # Multi processing
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=threads)
    future_generate_exon_indicator = {executor.submit(generate_exon_indicator_for_isoform_single_chr, raw_isoform_exons_dict[chr_name],gene_exons_dict[chr_name],gene_points_dict[chr_name],READ_LEN,READ_JUNC_MIN_MAP_LEN): chr_name for chr_name in raw_isoform_exons_dict}
    for future in concurrent.futures.as_completed(future_generate_exon_indicator):
        chr_name = future_generate_exon_indicator[future]
        try:
            isoforms_regions_len_dict[chr_name],gene_regions_dict[chr_name],genes_regions_len_dict[chr_name] = future.result()
        except Exception as exc:
            print('%r generated an exception: %s' % (chr_name, exc))
    # No parralize
    # for chr_name in raw_isoform_exons_dict:
    #     for gene_name in raw_isoform_exons_dict[chr_name]:
    #         isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name] = generate_exon_indicator_for_isoform_single_gene(raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],READ_LEN,READ_JUNC_MIN_MAP_LEN)
    return isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict
# def sanity_check_isoform_regions_length(gene_exons_dict, gene_regions_dict, gene_isoforms_dict, 
#                                         isoforms_regions_len_dict):
#     genes_regions_len_dict = {}
#     # Sanity check the isoforms regions length
#     for i in list(gene_exons_dict.keys()):
#         genes_regions_len_dict[i] = {}
#         for j in list(gene_exons_dict[i].keys()):
#             genes_regions_len_dict[i][j] = {}
#             for k in gene_regions_dict[i][j]:
#                 if k not in isoforms_regions_len_dict[i][j]:
#                     continue
#                 region_len = 0
#                 for m in gene_isoforms_dict[i][j]:
#                     if m not in isoforms_regions_len_dict[i][j][k]:
#                         continue
#                     if (region_len == 0):
#                         if (isoforms_regions_len_dict[i][j][k][m] > 0):
#                             genes_regions_len_dict[i][j][k] = isoforms_regions_len_dict[i][j][k][m]
#                             region_len = genes_regions_len_dict[i][j][k]
#     return genes_regions_len_dict

#######################
def parse_annotation(ref_annotation_path,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    #gene_points_dict store the index of the point value in ascending order for each gene
    file_read = open(ref_annotation_path, 'r')
    gene_exons_dict,gene_points_dict,gene_isoforms_dict,gene_isoforms_length_dict,raw_isoform_exons_dict = {},{},{},{},{}
    num_exons = 0
    for line in file_read:
        if line.lstrip()[0] == "#":
            continue
        fields = line.split('\t')
        if (fields[2] != 'exon'):
            continue
        num_exons += 1
        chr_name = fields[0]
        gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
        isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
        # use 1 based index
        start_pos = int(fields[3])
        end_pos = int(fields[4])
        
        #initialize dict
        if chr_name not in gene_exons_dict:
            gene_exons_dict[chr_name],gene_points_dict[chr_name],gene_isoforms_dict[chr_name],gene_isoforms_length_dict[chr_name],raw_isoform_exons_dict[chr_name] = {},{},{},{},{}
        if gene_name not in gene_exons_dict[chr_name]:
            gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],gene_isoforms_dict[chr_name][gene_name],gene_isoforms_length_dict[chr_name][gene_name],raw_isoform_exons_dict[chr_name][gene_name]= [],{},[],{},{}
        gene_exons_dict[chr_name][gene_name].append([start_pos, end_pos])
        if isoform_name not in gene_isoforms_length_dict[chr_name][gene_name]:
            gene_isoforms_length_dict[chr_name][gene_name][isoform_name] = 0
            raw_isoform_exons_dict[chr_name][gene_name][isoform_name] = {'region_pos':[]}
            gene_isoforms_dict[chr_name][gene_name].append(isoform_name)
        # note here the base on both end included in our system
        gene_isoforms_length_dict[chr_name][gene_name][isoform_name] += end_pos - start_pos + 1
        raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos'].append([start_pos,end_pos])
    file_read.close()
    for chr_name in raw_isoform_exons_dict:
        for gene_name in raw_isoform_exons_dict[chr_name]:
            for isoform_name in raw_isoform_exons_dict[chr_name][gene_name]:
                region_pos = sorted(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos'],key=itemgetter(0, 1))
                raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'] = [start_pos for [start_pos,end_pos] in region_pos]
                raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'] = [end_pos for [start_pos,end_pos] in region_pos]
                del raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['region_pos']
    gene_exons_dict = split_and_sort_exons(gene_exons_dict)
    # index the point position
    for chr_name in gene_exons_dict:
        for gene_name in gene_exons_dict[chr_name]:
            point_index = 0
            for [start_pos,end_pos,_] in gene_exons_dict[chr_name][gene_name]:
                gene_points_dict[chr_name][gene_name][start_pos] = point_index
                gene_points_dict[chr_name][gene_name][end_pos] = point_index + 1
                point_index += 2
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = generate_exon_indicator_for_isoform(gene_exons_dict, gene_points_dict, raw_isoform_exons_dict,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    return [gene_exons_dict, gene_points_dict, gene_isoforms_dict,genes_regions_len_dict,
            isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict]
