#!/usr/bin/python
import sys
import re
from operator import itemgetter, attrgetter
import bisect
import traceback
from util import sync_reference_name
import config
# from memory_profiler import profile
# valid_cigar = set("0123456789MNID")
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
def parse_read_line(line ):
    
    fields = line.split('\t')
    read_name = fields[0]
    rname = fields[2]
    converted_chr_name = sync_reference_name(rname)
    if (converted_chr_name.isnumeric()):
        rname = converted_chr_name
    read_start_pos = int(fields[3])
    cigar_field = fields[5]
    read_len_list = []
    mapped_read_len = 0
    # if (len(set(cigar_field) - valid_cigar) > 0):
    #     cigar_field = '*'
    if (cigar_field != '*'):
        cigar_list = re.split(r'(M|N|I|D|S|=|X|H|P)', cigar_field)
        read_len_list = []
        seg_len = 0
        M = 1
        for idx in range(len(cigar_list)//2):
            if (cigar_list[2 * idx + 1] in ['M','=','X']):
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                seg_len += int(cigar_list[2 * idx])
                mapped_read_len += int(cigar_list[2 * idx])
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
                seg_len += int(cigar_list[2 * idx])
                M = 1
            elif (cigar_list[2 * idx + 1] in ['I','S','H']):  # Insertion in reference
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                mapped_read_len +=  int(cigar_list[2 * idx])
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
def map_read_to_region(read_start_pos,read_len_list,points_dict,gene_interval_tree,gene_region_dict,read_name,READ_JUNC_MIN_MAP_LEN):
    tolerance = config.isoform_start_end_site_tolerance
    junc_tolerance = config.junction_site_tolerance
    read_segments = []
    curr_pos = read_start_pos
    for i in range(len(read_len_list)):
        if i % 2 == 0:
            read_segments.append([curr_pos,curr_pos+read_len_list[i]-1])
        curr_pos += read_len_list[i]
    all_possible_exon_regions = []
    for i in range(len(read_segments)):
        [seg_start,seg_end] = read_segments[i]
        if seg_end - seg_start < READ_JUNC_MIN_MAP_LEN:
            return '',0,1
    for j in range(len(read_segments)):
        [seg_start,seg_end] = read_segments[j]
        exons = [[exon.begin,exon.end - 1] for exon in gene_interval_tree.overlap(seg_start-junc_tolerance,seg_end+junc_tolerance)]
        exons = sorted(exons,key=lambda exon:exon[0],reverse=False)
        possible_exon_regions = []
        if len(exons) > 0:
            exons = sorted(exons,key=lambda exon:exon[0],reverse=True)
            for exon_start_idx in range(len(exons)):
                temp_exon_regions = []
                [exon_start,exon_end] = exons[exon_start_idx]
                if j != len(read_segments) - 1:
                    if not abs(exon_end - seg_end) <= junc_tolerance:
                        continue
                if j == len(read_segments) - 1:
                    if exon_end < seg_end and seg_end - exon_end > tolerance:
                        continue
                exon_region_name = 'P{}:P{}'.format(points_dict[exon_start],points_dict[exon_end])
                temp_exon_regions.append({'start':exon_start,'end':exon_end,'length':exon_end - exon_start + 1,'region':exon_region_name})
                curr_region_start,curr_region_end = exon_start,exon_end
                for [exon_start,exon_end] in exons[exon_start_idx + 1:]:
                    if exon_end == curr_region_start:
                        exon_region_name = 'P{}:{}'.format(points_dict[exon_start],exon_region_name)
                        curr_region_start = exon_start
                        temp_exon_regions.append({'start':exon_start,'end':curr_region_end,'length':curr_region_end - exon_start + 1,'region':exon_region_name})
                    elif curr_region_start - exon_end <= junc_tolerance:
                        exon_region_name = 'P{}:P{}-{}'.format(points_dict[exon_start],points_dict[exon_end],exon_region_name)
                        curr_region_start = exon_start
                        temp_exon_regions.append({'start':exon_start,'end':curr_region_end,'length':curr_region_end - exon_start + 1,'region':exon_region_name})
                    else:
                        break
                for exon_region in temp_exon_regions:
                    exon_start = exon_region['start']
                    if j!= 0:
                        if not (abs(exon_start - seg_start) <= junc_tolerance):
                            continue
                    if j == 0:
                        if exon_start > seg_start and exon_start - seg_start > tolerance:
                            continue
                    possible_exon_regions.append(exon_region)
        all_possible_exon_regions.append(possible_exon_regions)
    for possible_exon_regions in all_possible_exon_regions:
        if len(possible_exon_regions) == 0:
            return '',0,1
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
            best_regions = best_regions
        else:
            new_best_regions = []
            for (connected_region,overlapped_length,mapped_region_length) in best_regions:
                for region_dict in possible_exon_regions:
                    if int(connected_region.split(':')[-1][1:]) < int(region_dict['region'].split(':')[0][1:]):
                        new_connected_region = '{}-{}'.format(connected_region,region_dict['region'])
                        new_overlapped_length = overlapped_length + compute_overlapped_length(region_dict,seg_start,seg_end)
                        new_mapped_region_length = mapped_region_length + region_dict['length']
                        new_best_regions.append((new_connected_region,new_overlapped_length,new_mapped_region_length))
            best_regions = sorted(new_best_regions,key=lambda x:(x[1],x[1]/x[2]),reverse=True)
            best_regions = best_regions
    for (connected_region,overlapped_length,mapped_region_length) in best_regions:
        if connected_region in gene_region_dict:
            return connected_region,overlapped_length,mapped_region_length
    
    return '',0,1


    
            


    
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
             READ_JUNC_MIN_MAP_LEN, CHR_LIST,parsed_line):
    # mapping = {}
    [read_name, read_start_pos, rname, read_len_list] = parsed_line
    tolerance = 20
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
    mapping['read_pos'] = (read_start_pos,read_end_pos)
    if read_end_pos - read_start_pos > 2 * tolerance:
        start_index = bisect.bisect_right(start_pos_list[rname], read_start_pos+tolerance)
        end_index = bisect.bisect_left(end_pos_list[rname], read_end_pos-tolerance)
    else:
        start_index = bisect.bisect_right(start_pos_list[rname], read_start_pos+2)
        end_index = bisect.bisect_left(end_pos_list[rname], read_end_pos-2)
    gene_candidates = (set(end_gname_list[rname][end_index:]) & set(start_gname_list[rname][:start_index])) 
    best_overlapped_length = 0
    best_mapped_region_length = 1
    best_regions = []
    best_genes = []
    gene_candidates = sorted(list(gene_candidates))
    for gname in gene_candidates:
        points_dict = gene_points_dict[rname][gname]
        gene_interval_tree = gene_interval_tree_dict[rname][gname]
        temp_region,temp_overlapped_length,temp_mapped_region_length= map_read_to_region(read_start_pos,read_len_list,points_dict,gene_interval_tree,gene_regions_dict[rname][gname],read_name,READ_JUNC_MIN_MAP_LEN)

        if temp_region == '':
            continue
        if abs(temp_overlapped_length-best_overlapped_length) <= 2:
            if abs(temp_mapped_region_length - best_mapped_region_length) <= 2:
                best_regions.append(temp_region)
                best_overlapped_length = temp_overlapped_length
                best_mapped_region_length = temp_mapped_region_length
                best_genes.append(gname)
            elif temp_mapped_region_length < best_mapped_region_length:
                best_regions = [temp_region]
                best_overlapped_length = temp_overlapped_length
                best_mapped_region_length = temp_mapped_region_length
                best_genes = [gname]                
        else:
            if temp_overlapped_length > best_overlapped_length:
                best_regions = [temp_region]
                best_overlapped_length = temp_overlapped_length
                best_mapped_region_length = temp_mapped_region_length
                best_genes = [gname]
    if (len(best_regions) !=0):
        mapping['read_mapped'] = True
        is_false_mapped = True
        for best_gene,best_region in zip(best_genes,best_regions):
            mapping['mapping_area'].append({'chr_name':rname,'gene_name':best_gene,'region_name':best_region})
    return mapping