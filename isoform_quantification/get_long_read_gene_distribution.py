from collections import defaultdict
from intervaltree import IntervalTree
from parse_alignment import parse_read_line,comp_read_len
import re
def get_gene_range(annot_path):
    gene_range_dict = defaultdict(lambda:IntervalTree())
    with open(annot_path,'r') as f:
        for line in f.readlines():
            if line[0] == '#':
                continue
            fields = line.split('\t')
            if (fields[2] != 'gene'):
                continue
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            gene_name = re.findall('gene_id "([^"]*)"', fields[8])[0]
            gene_range_dict[fields[0]].addi(start_pos, end_pos+1, gene_name)
    return gene_range_dict
def get_long_read_gene_distribution(annot_path,alignment_path):
    gene_range_dict = get_gene_range(annot_path)
    gene_read_len_dict = defaultdict(lambda:defaultdict(lambda:[]))
    with open(alignment_path,'r') as f:
        for line in f:
            if line[0] == '@':
                continue
            fields = line.split('\t')
            if (fields[2] == '*'):
                continue
            [read_name, read_start_pos, rname, read_len_list] = parse_read_line(line,0)
            read_len = len(fields[9])
            read_end_pos = read_start_pos + read_len
            gene_candidates = gene_range_dict[rname].overlap(read_start_pos,read_end_pos)
            max_overlapped_length = 0
            max_gene = ''
            for gene in gene_candidates:
                total_length = max(read_end_pos,gene.end-1) - min(read_start_pos,gene.begin) + 1
                overlapped_length = (gene.end-1) - gene.begin + 1  + read_len - total_length
                if overlapped_length > max_overlapped_length:
                    max_gene = gene.data
                    max_overlapped_length = overlapped_length
            if max_gene != '':
                gene_read_len_dict[rname][max_gene].append(max_overlapped_length)
    
    gene_read_min_len_dict = defaultdict(dict)
    for rname in gene_read_len_dict:
        for gname in gene_read_len_dict[rname]:
            if (len(gene_read_len_dict[rname][gname])==0):
                continue
            sorted_lens = sorted(gene_read_len_dict[rname][gname])
            gene_read_min_len_dict[rname][gname] = sorted_lens[len(sorted_lens)//10]
    return gene_read_min_len_dict