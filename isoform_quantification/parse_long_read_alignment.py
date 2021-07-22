#!/usr/bin/python
import sys
import re
from operator import itemgetter, attrgetter
import bisect
import traceback

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
##TODO: BAD　need to be refactored
##########
def parse_read_line(line, READ_LEN):
    
    fields = line.split('\t')
    read_name = fields[0]
    read_start_pos = int(fields[3])
    rname = fields[2]
    cigar_field = fields[5]
    read_len_list = []
    read_len = 0
    if (len(set(cigar_field) - valid_cigar) > 0):
        cigar_field = '*'
    if (cigar_field != '*'):
        cigar_list = re.split(r'(M|N|I|D)', cigar_field)
        read_len_list = []
        seg_len = 0
        read_len = 0 
        M = 1
        for idx in range(len(cigar_list)//2):
            if (cigar_list[2 * idx + 1] == 'M'):
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
                seg_len += int(cigar_list[2 * idx])
            elif (cigar_list[2 * idx + 1] == 'I'):  # Insertion in reference
                if (M == 0):  # Mode is changed
                    read_len_list.append(seg_len)
                    seg_len = 0
                read_len +=  int(cigar_list[2 * idx])  
        read_len_list.append(seg_len)                
    else:
        read_len_list = []
        
#     if (abs(read_len - READ_LEN) > read_len_margin):
#         read_len_list = []
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
        print('Warning: Invalid CIGAR format: ' + str(read_len_list)) 
        
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
##########
    
    
### Map read to junc regions
### Algorithm:
###     It checks the start and end point of each read segment is inside of an exon region.
###     And checks gene Pi points inside a read segmnet maps either to the start/end point if
###     they are next to a junction gap.
##########
def map_read_to_junct_region(read_start_pos, read_len_list, points):
    
    region_name = ''
    read_len_idx = 0
    read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
    
    points_idx = 0
    read_len_idx += 1

    while (True):
        while (points[points_idx] < read_start_pos ):   # find the start position exon index
            points_idx += 1
#         if not (check_start_pos_in_exon_boundary(read_start_pos, points_idx, points)):
#             return ''
        while(points[points_idx] <= read_end_pos):
            if not (((points_idx + 1) < len(points)) and 
                    (points[points_idx] == points[points_idx+1])):   # Special case when the exon length is 1 we dont want to repeat the same point (should be the end idx)
                region_name += 'P' + str(points_idx) + '-'
#             if (check_read_contains_gap(points, points_idx, read_start_pos, read_end_pos)):
#                 return ''    # Not a valid read for this genome
            points_idx += 1
            if (points_idx >= len(points)):
                break       

            
        if (read_len_idx == len(read_len_list)):
            if (region_name != ''):
#                 if not (check_end_pos_in_exon_boundary(read_end_pos, points_idx, points)):
#                     return ''
                region_name = region_name[:-1]
            return region_name
        read_start_pos = read_end_pos + read_len_list[read_len_idx] + 1 
        read_len_idx += 1
        read_end_pos = read_start_pos + read_len_list[read_len_idx] - 1
        read_len_idx += 1
                
##########


### Map read to gene regions
##########
def map_read(line, gene_regions_read_count, gene_regions_points_list, 
             start_pos_list, start_gname_list, end_pos_list, end_gname_list,
             READ_LEN, READ_JUNC_MIN_MAP_LEN, CHR_LIST):
    mapping = {}
    [read_name, read_start_pos, rname, read_len_list] = parse_read_line(line, READ_LEN)
    mapping = {'read_name':read_name,'read_start_pos':read_start_pos,'rname':rname,'read_len':read_len_list,'mapping_area':[],'read_mapped':False}
    
    if (rname not in CHR_LIST):
        return mapping
    
    if ((len(read_len_list) % 2) == 0):  # CIGAR should start and end in non-gap tag
        return mapping
    # Ensure minimum READ_JUNC_MIN_MAP_LEN on each end
    if ((read_len_list[0] < READ_JUNC_MIN_MAP_LEN) or
        (read_len_list[-1] < READ_JUNC_MIN_MAP_LEN)):  
        return mapping
    read_end_pos = read_start_pos  - 1
    for i in read_len_list:
        read_end_pos += i
    
    gene_candidates = []
    start_index = bisect.bisect_right(end_pos_list[rname], read_start_pos)
    end_index = bisect.bisect_left(start_pos_list[rname], read_end_pos)
    gene_candidates = (set(end_gname_list[rname][:end_index]) & set(start_gname_list[rname][start_index:])) 
    for gname in gene_candidates:
        points = gene_regions_points_list[rname][gname]
        region_name = map_read_to_exon_region(read_start_pos, read_len_list, points)
        if (region_name == ''):
            region_name =  map_read_to_junct_region(read_start_pos, read_len_list, points)
        if (region_name != ''):
            if (region_name in gene_regions_read_count[rname][gname]):
                gene_regions_read_count[rname][gname][region_name] += 1
                mapping['read_mapped'] = True
                mapping['mapping_area'].append({'gene_name':gname,'region_name':region_name})
                #print 'Gname %s, %s, region %s mapped' % (rname, gname, region_name)
                #print 'Read: ' + line
            #else:
                #print 'Warning: Gname %s, %s, does not contain region %s ' % (gname, rname, region_name)
                #print '         Read: ' + line
        #else:
            #print 'Warning: Gname %s, %s, was not mapped.' % (gname, rname)
            #print '         Read: ' + line
    return mapping
##########        
