#!/usr/bin/python
from operator import itemgetter, attrgetter
import re
import sys
import time
import numpy as np
import concurrent.futures
from util import sync_reference_name
import pickle
##########
def split_and_sort_exons(gene_exons_dict):
    new_gene_exons_dict = {}
    for chr_name in gene_exons_dict:
        new_gene_exons_dict[chr_name] = {}
        for gene_name in gene_exons_dict[chr_name]:
            exon_sorted = sorted(gene_exons_dict[chr_name][gene_name], key=itemgetter(0, 1))  # Sort by start position then end position
            exon_points = []
            for [begin,end] in exon_sorted:
                exon_points.append((begin,'+'))
                exon_points.append((end,'-'))
            exon_points = sorted(exon_points,key=lambda x:x[0])
            new_exon_sorted = [] 
            last_begin = None
            for offset,pm in exon_points:
                if pm == '+':
                    if last_begin is not None:
                        new_exon_sorted.append([last_begin,offset])
                    last_begin = offset
                elif pm == '-':
                    new_exon_sorted.append([last_begin,offset])
                    last_begin = offset
            if last_begin is not None:
                new_exon_sorted.append([last_begin,offset])
            # with open('/fs/ess/scratch/PCON0009/haoran/IDP/run_TransELS/TransELS/H1.pkl','wb') as f:
            #     pickle.dump(new_exon_sorted,f)
            new_exon_sorted = sorted(new_exon_sorted, key=itemgetter(0, 1))
            new_gene_exons_dict[chr_name][gene_name] = [[start_pos,end_pos,end_pos - start_pos + 1] for [start_pos,end_pos] in new_exon_sorted if start_pos != end_pos]
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
##########
def generate_exon_indicator_for_isoform_single_gene(args):
    single_gene_raw_isoform_exons_dict,exon_list,single_gene_gene_points_dict,READ_LEN,READ_JUNC_MIN_MAP_LEN = args
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
                    # if (l >= READ_LEN):  # A valid region for exon
                    temp_isoform_name = set()
                    temp_isoform_name.add(isoform_name)
                    if (region_name in single_gene_gene_regions_dict):
                        temp_isoform_name = temp_isoform_name.union(single_gene_gene_regions_dict[region_name])
                    single_gene_gene_regions_dict[region_name] = temp_isoform_name
                    isoform_exons.append(region_name)
                    flag = False
                elif (p0 <= start_pos[j]):
                    flag = False
                else:
                    j += 1
                    flag = (j < num_exons) 
        # Check the junction regions
        for i in range(len(exon_list) - 1):     # the last exon can not have any junction
            p0 = exon_list[i][0]
            p1 = exon_list[i][1]
            l  = exon_list[i][2]
            exon_region_name = 'P' + str(single_gene_gene_points_dict[p0]) + ':P' + str(single_gene_gene_points_dict[p1])
            if (exon_region_name not in isoform_exons):   # this exon is not in the isoform
                continue
            if (l<READ_JUNC_MIN_MAP_LEN):
                continue
            region_name = 'P{}:P{}'.format(str(single_gene_gene_points_dict[p0]),str(single_gene_gene_points_dict[p1]))
            for j in range(i+1,len(exon_list)):
                p0_temp = exon_list[j][0]
                p1_temp = exon_list[j][1]
                l_temp  = exon_list[j][2]
                exon_region_name = 'P' + str(single_gene_gene_points_dict[p0_temp]) + ':P' + str(single_gene_gene_points_dict[p1_temp])
                if (exon_region_name not in isoform_exons):   # this exon is not in the isoform
                    continue
                if (region_name.split(':')[-1] == exon_region_name.split(':')[0]):
                    region_name += ':P{}'.format(str(single_gene_gene_points_dict[p1_temp]))
                else:
                    region_name += '-P{}:P{}'.format(str(single_gene_gene_points_dict[p0_temp]),str(single_gene_gene_points_dict[p1_temp]))
                temp_isoform_set = {isoform_name}
                if (region_name in single_gene_gene_regions_dict):
                    temp_isoform_set = temp_isoform_set.union(single_gene_gene_regions_dict[region_name])
                single_gene_gene_regions_dict[region_name] = temp_isoform_set
    single_gene_regions_len_dict = {}
    point_coord_dict = {}
    for p in single_gene_gene_points_dict:
        point_coord_dict['P{}'.format(single_gene_gene_points_dict[p])] = int(p)
    # filter out region with short ending exon
    region_names = list(single_gene_gene_regions_dict.keys())
    for region in region_names:
        if '-' in region:
            last_exon = region.split('-')[-1]
            last_exon_len = cal_region_length(last_exon,point_coord_dict)
            if last_exon_len < READ_JUNC_MIN_MAP_LEN:
                del single_gene_gene_regions_dict[region]
    for region in single_gene_gene_regions_dict:
        region_len = cal_region_length(region,point_coord_dict)
        single_gene_regions_len_dict[region] = region_len
        # if region is an exon
        if (region.count(':') == 1 and '-' not in region):
            for isoform_name in single_gene_gene_regions_dict[region]:
                if isoform_name in single_gene_isoforms_regions_len_dict:
                    single_gene_isoforms_regions_len_dict[isoform_name] += region_len
                else:
                    single_gene_isoforms_regions_len_dict[isoform_name] = region_len

    return single_gene_isoforms_regions_len_dict,single_gene_gene_regions_dict,single_gene_regions_len_dict
def generate_exon_indicator_for_isoform(gene_exons_dict,gene_points_dict,raw_isoform_exons_dict,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN):
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = {},{},{}
    for chr_name in raw_isoform_exons_dict:
        isoforms_regions_len_dict[chr_name],gene_regions_dict[chr_name],genes_regions_len_dict[chr_name] = {},{},{}
        
    if threads == 1:
        for chr_name in raw_isoform_exons_dict:
            for gene_name in raw_isoform_exons_dict[chr_name]:
                isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = generate_exon_indicator_for_isoform_single_gene((raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],READ_LEN,READ_JUNC_MIN_MAP_LEN))
    else:
        list_of_all_genes_chrs = [(gene_name,chr_name) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
        list_of_args = [(raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],READ_LEN,READ_JUNC_MIN_MAP_LEN) for chr_name in raw_isoform_exons_dict for gene_name in raw_isoform_exons_dict[chr_name]]
        chunksize, extra = divmod(len(list_of_all_genes_chrs), threads)
        if extra:
            chunksize += 1
        print('Using {} threads'.format(threads))
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            for (gene_name,chr_name), result in zip(list_of_all_genes_chrs, executor.map(generate_exon_indicator_for_isoform_single_gene,list_of_args,chunksize=chunksize)):
                isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = result
    # No parralize
    # for chr_name in raw_isoform_exons_dict:
    #     for gene_name in raw_isoform_exons_dict[chr_name]:
    #         isoforms_regions_len_dict[chr_name][gene_name],gene_regions_dict[chr_name][gene_name],genes_regions_len_dict[chr_name][gene_name] = generate_exon_indicator_for_isoform_single_gene((raw_isoform_exons_dict[chr_name][gene_name],gene_exons_dict[chr_name][gene_name],gene_points_dict[chr_name][gene_name],READ_LEN,READ_JUNC_MIN_MAP_LEN))
    return isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict
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
        chr_name = fields[0]
        converted_chr_name = sync_reference_name(fields[0])
        if (converted_chr_name.isnumeric()):
            chr_name = converted_chr_name
        # if chr_name == '':
        #     continue
        num_exons += 1
        gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
        isoform_name = re.findall('transcript_id "([^"]*)"', fields[8])[0]
        # use 1 based index
        start_pos = int(fields[3])
        end_pos = int(fields[4])
        if start_pos >= end_pos:
            continue
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
    raw_gene_exons_dict = gene_exons_dict.copy()
    gene_exons_dict = split_and_sort_exons(gene_exons_dict)
    # index the point position
    for chr_name in gene_exons_dict:
        for gene_name in gene_exons_dict[chr_name]:
            point_index = 0
            for [start_pos,end_pos,_] in gene_exons_dict[chr_name][gene_name]:
                for pos in [start_pos,end_pos]:
                    if pos not in gene_points_dict[chr_name][gene_name]:
                        gene_points_dict[chr_name][gene_name][pos] = point_index
                        point_index += 1
    isoforms_regions_len_dict,gene_regions_dict,genes_regions_len_dict = generate_exon_indicator_for_isoform(gene_exons_dict, gene_points_dict, raw_isoform_exons_dict,threads,READ_LEN,READ_JUNC_MIN_MAP_LEN)
    return [gene_exons_dict, gene_points_dict, gene_isoforms_dict,genes_regions_len_dict,
            isoforms_regions_len_dict, gene_regions_dict, gene_isoforms_length_dict,raw_isoform_exons_dict,raw_gene_exons_dict]
