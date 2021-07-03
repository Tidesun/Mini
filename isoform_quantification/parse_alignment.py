#!/usr/bin/python
import sys
import re
from operator import itemgetter, attrgetter
import bisect
import traceback
from util import sync_reference_name
# from memory_profiler import profile
valid_cigar = set("0123456789MNID")
read_len_margin = 0

### Adds missing starting points when exon length is 1
##########
def update_missing_points(points_idx, points):
    
    points_temp = []
    last_idx = int(points_idx[-1][1:])
    idx = 0
    i = 0
    while (idx <= last_idx):
        points_temp.append(points[i])
        if ('P' + str(idx) in points_idx):
            i += 1
        idx += 1
    return points_temp


### Extract the SAM read information
##########
def parse_read_line(line, READ_LEN):
    
    fields = line.split('\t')
    read_name = fields[0]
    rname = fields[2]
    converted_chr_name = sync_reference_name(rname)
    if (converted_chr_name.isnumeric()):
        rname = converted_chr_name
    read_start_pos = int(fields[3])
    cigar_field = fields[5]
    read_len_list = []
    read_len = 0
    # if (len(set(cigar_field) - valid_cigar) > 0):
    #     cigar_field = '*'
    if (cigar_field != '*'):
        cigar_list = re.split(r'(M|N|I|D|S|=|X|H|P)', cigar_field)
        read_len_list = []
        seg_len = 0
        read_len = 0 
        M = 1
        for idx in range(len(cigar_list)//2):
            if (cigar_list[2 * idx + 1] in ['M','=','X']):
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                seg_len += int(cigar_list[2 * idx])
                read_len += int(cigar_list[2 * idx])
                M = 1
            elif (cigar_list[2 * idx + 1] == 'N'):
                if (M == 1):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                seg_len += int(cigar_list[2 * idx])
                M = 0
            elif (cigar_list[2 * idx + 1] == 'D'):  # Deletion from reference
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                else:
                    seg_len += int(cigar_list[2 * idx])
            elif (cigar_list[2 * idx + 1] in ['I','S']):  # Insertion in reference
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                read_len +=  int(cigar_list[2 * idx])
        read_len_list.append(seg_len)                
    else:
        read_len_list = []
        
#     if (abs(read_len - READ_LEN) > read_len_margin):
#         read_len_list = []
    read_len_list = [i for i in read_len_list if i!=0 ]
    return [read_name, read_start_pos, rname, read_len_list]
##########


### Compute the mapped read length
##########
def comp_read_len(read_len_list):
    
    M= 1
    read_len = 0
    for i in read_len_list:
        if (M == 1):
            read_len += i
        M = 1 - M
    if M == 1:
        raise Exception('Warning: Invalid CIGAR format: ' + str(read_len_list)) 
        
    return read_len
##########


### Check if read contains junction gap (Which cause in invalid mapping)
##########
def check_read_contains_gap(points, points_idx, 
                            read_start_pos, read_end_pos):
    
    if (((points_idx % 2) == 0) and
        (points_idx > 0)):
        if (points[points_idx] != (points[points_idx-1] + 1)):    # this is not an extended exon
            return points[points_idx] != read_start_pos
        else:
            return False
        
    if (((points_idx % 2) == 1) and
        ((points_idx + 1) < len(points) )):
        if (points[points_idx] != (points[points_idx+1] - 1)):    # this is not an extended exon
            return points[points_idx] != read_end_pos
        else:
            return False

    return False 
#########

### Check if start read position is in exon boundary (point[point_idx] >= point) 
##########
def check_start_pos_in_exon_boundary(pos, points_idx, points):
    
    if ((points_idx % 2) == 0):
        return points[points_idx] == pos     
    else:
        return (pos > points[points_idx - 1])
### Check if end read position is in exon boundary (point[point_idx] > point) 
##########
def check_end_pos_in_exon_boundary(pos, points_idx, points):
    
    if ((points_idx % 2) == 0):
        return points[points_idx -1] == pos     
    else:
        return ((pos == points[points_idx -1]) or (pos < points[points_idx]))   # Note: If reaching the last index, it should be equal to end read pos 
    
### Map read to gene regions
##########
def map_read_to_exon_region(read_start_pos, read_len_list, points):
    
    region_name = ''
    if (len(read_len_list) > 1):       # read definitely maps to multiple junction
        return   region_name


    read_end_pos = read_start_pos + read_len_list[0] - 1
    for i in range(len(points)//2):
        if ((read_start_pos >= points[2*i]) and 
            (read_end_pos <= points[2*i+1])):
            region_name = 'P' + str(2*i) + ':P' + str(2*i+1)
            return region_name

    return region_name 
def compute_overlapped_length(region_dict,seg_start,seg_end):
    seg_length = seg_end - seg_start + 1
    total_length = max(seg_end,region_dict['end']) - min(seg_start,region_dict['start']) + 1
    overlapped_length = region_dict['length'] + seg_length - total_length
    # overlapped_prop = overlapped_length/region_dict['length']
    return overlapped_length

##########
def map_read_to_region(read_start_pos,read_len_list,points_dict,gene_interval_tree,gene_region_dict,read_name):
    read_segments = []
    curr_pos = read_start_pos
    for i in range(len(read_len_list)):
        if i % 2 == 0:
            read_segments.append([curr_pos,curr_pos+read_len_list[i]-1])
        curr_pos += read_len_list[i]
    all_possible_exon_regions = []
    # best_region = ''
    # best_overlapped_length = 0
    for j in range(len(read_segments)):
        [seg_start,seg_end] = read_segments[j]
        exons = [[exon.begin,exon.end - 1] for exon in gene_interval_tree.overlap(seg_start,seg_end)]
        possible_exon_regions = []
        if len(exons) > 0:
            exons = sorted(exons,key=lambda exon:exon[0],reverse=True)
            for exon_start_idx in range(len(exons)):
                [exon_start,exon_end] = exons[exon_start_idx]
                exon_region_name = 'P{}:P{}'.format(points_dict[exon_start],points_dict[exon_end])
                possible_exon_regions.append({'start':exon_start,'end':exon_end,'length':exon_end - exon_start + 1,'region':exon_region_name})
                curr_region_start,curr_region_end = exon_start,exon_end
                for [exon_start,exon_end] in exons[exon_start_idx + 1:]:
                    if exon_end == curr_region_start:
                        exon_region_name = 'P{}:{}'.format(points_dict[exon_start],exon_region_name)
                        possible_exon_regions.append({'start':exon_start,'end':curr_region_end,'length':curr_region_end - exon_start + 1,'region':exon_region_name})
                        curr_region_start = exon_start
                    else:
                        break
        # max_overlapped_length = 0
        # max_overlapped_prop = 0
        # max_overlapped_region_dict = {}
        # for region_dict in possible_exon_regions:
        #     seg_length = seg_end - seg_start + 1
        #     total_length = max(seg_end,region_dict['end']) - min(seg_start,region_dict['start']) + 1
        #     overlapped_length = region_dict['length'] + seg_length - total_length
        #     overlapped_prop = overlapped_length/region_dict['length']
        #     if overlapped_length > max_overlapped_length or (overlapped_length == max_overlapped_length and overlapped_prop > max_overlapped_prop):
        #         max_overlapped_length = overlapped_length
        #         max_overlapped_prop = overlapped_prop
        #         max_overlapped_region_dict = region_dict
        # if 'region' in max_overlapped_region_dict:
        #     max_overlapped_region = max_overlapped_region_dict['region']
        # else:
        #     max_overlapped_region = ''
        
        # if (best_region == ''):
        #     best_region = max_overlapped_region
        #     best_overlapped_length = max_overlapped_length
        # else:
        #     if (best_region.split(':')[-1] != max_overlapped_region.split(':')[0]) and max_overlapped_region!= '':
        #         best_region = '{}-{}'.format(best_region,max_overlapped_region)
        #         best_overlapped_length += max_overlapped_length
        all_possible_exon_regions.append(possible_exon_regions)
    best_regions = []
    for read_segment,possible_exon_regions in zip(read_segments,all_possible_exon_regions):
        [seg_start,seg_end] = read_segment
        if len(best_regions) == 0:
            for region_dict in possible_exon_regions:
                new_connected_region = region_dict['region']
                new_overlapped_length =  compute_overlapped_length(region_dict,seg_start,seg_end)
                new_mapped_region_length = region_dict['length']
                best_regions.append((new_connected_region,new_overlapped_length,new_mapped_region_length))
            best_regions = sorted(best_regions,key=lambda x:(x[1],x[1]/x[2]),reverse=True)
            best_regions = best_regions[:2]
        else:
            new_best_regions = []
            for (connected_region,overlapped_length,mapped_region_length) in best_regions:
                for region_dict in possible_exon_regions:
                    if int(connected_region.split(':')[-1][1:]) < int(region_dict['region'].split(':')[0][1:]):
                        new_connected_region = '{}-{}'.format(connected_region,region_dict['region'])
                        new_overlapped_length = overlapped_length + compute_overlapped_length(region_dict,seg_start,seg_end)
                        new_mapped_region_length = mapped_region_length + region_dict['length']
                        new_best_regions.append((new_connected_region,new_overlapped_length,new_mapped_region_length))
            # best_regions = new_best_regions
            best_regions = sorted(new_best_regions,key=lambda x:(x[1],x[1]/x[2]),reverse=True)
            best_regions = best_regions[:2]
    for (connected_region,overlapped_length,mapped_region_length) in best_regions:
        if connected_region in gene_region_dict:
            return connected_region,overlapped_length
    # # if fail searching
    # for (connected_region,overlapped_length,mapped_region_length) in best_regions:
    #     if connected_region in region:
    #         return connected_region,overlapped_length
    
    return '',0


    
            


    
### Map read to junc regions
### Algorithm:
###     It checks the start and end point of each read segment is inside of an exon region.
###     And checks gene Pi points inside a read segmnet maps either to the start/end point if
###     they are next to a junction gap.
##########
def map_read_to_junct_region(read_start_pos, read_len_list, points,read_name):
    
    region_name = ''
    read_len_idx = 0
    read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
    
    points_idx = 0
    read_len_idx += 1

    while (True):
        while (points[points_idx] < read_start_pos ):   # find the start position exon index
            points_idx += 1
        while(points[points_idx] <= read_end_pos):
            if not (((points_idx + 1) < len(points)) and 
                    (points[points_idx] == points[points_idx+1])):   # Special case when the exon length is 1 we dont want to repeat the same point (should be the end idx)
                region_name += 'P' + str(points_idx) + '-'
            if (check_read_contains_gap(points, points_idx, read_start_pos, read_end_pos)):
                return ''    # Not a valid read for this genome
            points_idx += 1
            if (points_idx >= len(points)):
                break       

            
        if (read_len_idx == len(read_len_list)):
            if (region_name != ''):
                if not (check_end_pos_in_exon_boundary(read_end_pos, points_idx, points)):
                    return ''
                region_name = region_name[:-1]
            return region_name
        read_start_pos = read_end_pos + read_len_list[read_len_idx] + 1 
        read_len_idx += 1
        read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
        read_len_idx += 1
##########
### Map read to gene regions
##########
# @profile
def map_read(gene_points_dict,gene_interval_tree_dict,gene_regions_dict, 
             start_pos_list, start_gname_list, end_pos_list, end_gname_list,
             READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,parsed_line):
    # mapping = {}
    [read_name, read_start_pos, rname, read_len_list] = parsed_line
    # mapping = {'read_name':read_name,'read_start_pos':read_start_pos,'rname':rname,'read_len':read_len_list,'mapping_area':[],'read_mapped':False}
    mapping = {'read_name':read_name,'read_mapped':False,'mapping_area':[]}
    if (rname not in CHR_LIST):
        return mapping
    if ((len(read_len_list) % 2) == 0):  # CIGAR should start and end in non-gap tag
        return mapping
    read_end_pos = read_start_pos - 1
    for i in read_len_list:
        read_end_pos += i
    read_length = comp_read_len(read_len_list)
    mapping['read_length'] = read_length
    start_index = bisect.bisect_right(start_pos_list[rname], read_start_pos)
    end_index = bisect.bisect_left(end_pos_list[rname], read_end_pos)
    gene_candidates = (set(end_gname_list[rname][end_index:]) & set(start_gname_list[rname][:start_index])) 
    best_overlapped_length = 0
    best_regions = []
    best_genes = []
    for gname in gene_candidates:
        points_dict = gene_points_dict[rname][gname]
        gene_interval_tree = gene_interval_tree_dict[rname][gname]
        temp_region,temp_overlapped_length = map_read_to_region(read_start_pos,read_len_list,points_dict,gene_interval_tree,gene_regions_dict[rname][gname],read_name)
        if temp_region == '':
            continue
        if temp_overlapped_length > best_overlapped_length:
            best_regions = [temp_region]
            best_overlapped_length = temp_overlapped_length
            best_genes = [gname]
        elif temp_overlapped_length == best_overlapped_length:
            best_regions.append(temp_region)
            best_genes.append(gname)

        # region_name = map_read_to_exon_region(read_start_pos, read_len_list, points)
        # if (region_name == ''):
        #     region_name =  map_read_to_junct_region(read_start_pos, read_len_list, points,read_name)
    if (len(best_regions) !=0):
        mapping['read_mapped'] = True
        is_false_mapped = True
        for best_gene,best_region in zip(best_genes,best_regions):
            mapping['mapping_area'].append({'chr_name':rname,'gene_name':best_gene,'region_name':best_region})
    return mapping
                #print 'Gname %s, %s, region %s mapped' % (rname, gname, region_name)
                #print 'Read: ' + line
            #else:
                #print 'Warning: Gname %s, %s, does not contain region %s ' % (gname, rname, region_name)
                #print '         Read: ' + line
        #else:
            #print 'Warning: Gname %s, %s, was not mapped.' % (gname, rname)
            #print '         Read: ' + line
def map_long_read(gene_points_dict,gene_interval_tree_dict,gene_regions_dict, 
             start_pos_list, start_gname_list, end_pos_list, end_gname_list,
             READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST,parsed_line):
    # mapping = {}
    [read_name, read_start_pos, rname, read_len_list] = parsed_line
    # mapping = {'read_name':read_name,'read_start_pos':read_start_pos,'rname':rname,'read_len':read_len_list,'mapping_area':[],'read_mapped':False}
    mapping = {'read_name':read_name,'read_mapped':False,'mapping_area':[]}
    if (rname not in CHR_LIST):
        return mapping
    
    if ((len(read_len_list) % 2) == 0):  # CIGAR should start and end in non-gap tag
        return mapping
    read_end_pos = read_start_pos - 1
    for i in read_len_list:
        read_end_pos += i
    read_length = comp_read_len(read_len_list)
    mapping['read_length'] = read_length

    gene_candidates = []
    start_index = bisect.bisect_right(start_pos_list[rname], read_start_pos)
    end_index = bisect.bisect_left(end_pos_list[rname], read_end_pos)
    gene_candidates = (set(end_gname_list[rname][end_index:]) & set(start_gname_list[rname][:start_index])) 
    best_overlapped_length = 0
    best_regions = []
    best_genes = []
    temp_region = ''
    for gname in gene_candidates:
        points_dict = gene_points_dict[rname][gname]
        gene_interval_tree = gene_interval_tree_dict[rname][gname]
        temp_region,temp_overlapped_length = map_read_to_region(read_start_pos,read_len_list,points_dict,gene_interval_tree,gene_regions_dict[rname][gname],read_name)
        if temp_region == '':
            continue
        if temp_overlapped_length > best_overlapped_length:
            best_regions = [temp_region]
            best_overlapped_length = temp_overlapped_length
            best_genes = [gname]
        elif temp_overlapped_length == best_overlapped_length:
            best_regions.append(temp_region)
            best_genes.append(gname)

        # region_name = map_read_to_exon_region(read_start_pos, read_len_list, points)
        # if (region_name == ''):
        #     region_name =  map_read_to_junct_region(read_start_pos, read_len_list, points,read_name)
    if (len(best_regions) !=0):
        mapping['read_mapped'] = True
        for best_gene,best_region in zip(best_genes,best_regions):
            mapping['mapping_area'].append({'chr_name':rname,'gene_name':best_gene,'region_name':best_region})
    return mapping
##########        
# def map_long_read_to_region(read_start_pos, read_len_list, points):
#     region_name = ''
#     read_len_idx = 0
#     read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
    
#     points_idx = 0
#     read_len_idx += 1
#     while (True):
#         while (points[points_idx] < read_start_pos ):   # find the start position exon index
#             points_idx += 1
#         while(points[points_idx] <= read_end_pos):
#             if not (((points_idx + 1) < len(points)) and 
#                     (points[points_idx] == points[points_idx+1])):   # Special case when the exon length is 1 we dont want to repeat the same point (should be the end idx)
#                 region_name += 'P' + str(points_idx) + '-'
#             points_idx += 1
#             if (points_idx >= len(points)):
#                 break
#         if (read_len_idx == len(read_len_list)):
#             if (region_name != ''):
#                 region_name = region_name[:-1]
#             return region_name
#         read_start_pos = read_end_pos + read_len_list[read_len_idx] + 1 
#         read_len_idx += 1
#         read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
#         read_len_idx += 1
    
# def map_long_read(parsed_line, gene_regions_read_count, gene_regions_read_length,gene_points_dict, gene_interval_tree_dict,
#              start_pos_list, start_gname_list, end_pos_list, end_gname_list,
#              READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST):
#     mapping = {}
#     [read_name, read_start_pos, rname, read_len_list] = parsed_line
#     mapping = {'read_name':read_name,'read_start_pos':read_start_pos,'rname':rname,'read_len':read_len_list,'mapping_area':[],'read_mapped':False}
#     return mapping,gene_regions_read_count,gene_regions_read_length
#     if (rname not in CHR_LIST):
#         return mapping,gene_regions_read_count,gene_regions_read_length
#     try:
#         read_length = comp_read_len(read_len_list)
#     except:
#         print(read_name)
#         read_length = 150
#     read_end_pos = read_start_pos  - 1 + read_length
    
#     gene_candidates = []
#     start_index = bisect.bisect_right(end_pos_list[rname], read_start_pos)
#     end_index = bisect.bisect_left(start_pos_list[rname], read_end_pos)
#     gene_candidates = (set(end_gname_list[rname][:end_index]) & set(start_gname_list[rname][start_index:])) 
#     for gname in gene_candidates:
#         points = gene_regions_points_list[rname][gname]
#         try:
#             region_name = map_long_read_to_region(read_start_pos, read_len_list, points)
#         except:
#             region_name = ''
#         if (region_name != ''):
#             sub_region_name = region_name
#             # for region_candidate in gene_regions_read_count[rname][gname]:
#             #     if (region_candidate.replace(":","-") in region_name) & (len(region_candidate) > len(sub_region_name)):
#             #         sub_region_name = region_candidate
#             gene_regions_read_count[rname][gname][sub_region_name] += 1
#             # gene_regions_read_length[rname][gname][sub_region_name].append({'read_name':read_name,'read_length':read_length})
#             gene_regions_read_length[rname][gname][sub_region_name].append(read_length)
#             mapping['read_mapped'] = True
#             mapping['mapping_area'].append({'gene_name':gname,'region_name':sub_region_name})
#     return mapping,gene_regions_read_count,gene_regions_read_length