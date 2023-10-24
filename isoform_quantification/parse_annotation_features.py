import numpy as np
import scipy.stats
import pickle
def get_stats(arr):
    if len(arr) == 0:
        return np.array([np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan')])
    s = scipy.stats.describe(arr)
    return np.array([s.minmax[0],s.minmax[1],s.mean,s.variance,s.skewness,s.kurtosis])
def mergeIntervals(intervals):
    intervals.sort()
    stack = []
    # insert first interval into stack
    stack.append(intervals[0])
    for i in intervals[1:]:
        # Check for overlapping interval,
        # if interval overlap
        if stack[-1][0] <= i[0] <= stack[-1][-1]:
            stack[-1][-1] = max(stack[-1][-1], i[-1])
        else:
            stack.append(i)
    return stack
def cal_length(stacks):
    s = 0
    for stack in stacks:
        s += stack[1] - stack[0] + 1
    return s
def save_gene_structure_features(gene_exons_dict,raw_gene_exons_dict,gene_isoforms_length_dict,genes_regions_len_dict,output_path):
    features = []
    for rname in gene_exons_dict:
        for gname in gene_exons_dict[rname]:
            exons = raw_gene_exons_dict[rname][gname]
            isoforms = gene_isoforms_length_dict[rname][gname]
            regions = genes_regions_len_dict[rname][gname]
            gene_length = cal_length(mergeIntervals(exons))
            isoforms_lengths = list(isoforms.values())
            exon_lengths = [e[1]-e[0]+1 for e in exons]
            region_lengths = list(regions.values())
            features.append([gname,gene_length,exon_lengths,isoforms_lengths,region_lengths])
    # np.savez_compressed(f'gene_features.npz', features=np.array(features),gnames=np.array(genes))
    with open(f'{output_path}/temp/machine_learning/gene_features.pkl','wb') as f:
        pickle.dump(features,f)