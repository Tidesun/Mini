import pickle
import glob
import numpy as np
import scipy.stats
def get_stats(arr):
    if len(arr) == 0:
        return np.array([np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan'),np.float('nan')])
    s = scipy.stats.describe(arr)
    return np.array([s.minmax[0],s.minmax[1],s.mean,s.variance,s.skewness,s.kurtosis])
def extract_features(output_path):
    with open(f'{output_path}/temp/machine_learning/gene_community_id_dict.pkl','rb') as f:
        gene_community_id_dict = pickle.load(f)
    SR_feature_dict = {}
    fpaths = glob.glob(f'{output_path}/temp/machine_learning/SR_feature_dict_*')
    for fpath in fpaths:
        with open(fpath,'rb') as f:
            gene_SR_feature_dict = pickle.load(f)
        for gname in gene_SR_feature_dict:
            community_id = gene_community_id_dict[gname]
            if community_id in SR_feature_dict:
                SR_feature_dict[community_id]['unique_mapping'] += list(gene_SR_feature_dict[gname]['unique_mapping'])
                SR_feature_dict[community_id]['multi_mapping'] += list(gene_SR_feature_dict[gname]['multi_mapping'])
            else:
                SR_feature_dict[community_id] = {}
                SR_feature_dict[community_id]['unique_mapping'] = list(gene_SR_feature_dict[gname]['unique_mapping'] )
                SR_feature_dict[community_id]['multi_mapping'] = list(gene_SR_feature_dict[gname]['multi_mapping'] )
    LR_feature_dict = {}
    fpaths = glob.glob(f'{output_path}/temp/machine_learning/LR_feature_dict_*')
    for fpath in fpaths:
        with open(fpath,'rb') as f:
            gene_LR_feature_dict = pickle.load(f)
        for gname in gene_LR_feature_dict:
            community_id = gene_community_id_dict[gname]
            if community_id in LR_feature_dict:
                LR_feature_dict[community_id]['unique_mapping'] += gene_LR_feature_dict[gname]['unique_mapping']
                LR_feature_dict[community_id]['multi_mapping'] += gene_LR_feature_dict[gname]['multi_mapping']
                LR_feature_dict[community_id]['3_end_truncation'] += gene_LR_feature_dict[gname]['3_end_truncation']
                LR_feature_dict[community_id]['5_end_truncation'] += gene_LR_feature_dict[gname]['5_end_truncation']
            else:
                LR_feature_dict[community_id] = gene_LR_feature_dict[gname]
    
    features = []
    community_ids = []
    for community_id in set(SR_feature_dict.keys()).union(LR_feature_dict.keys()):
        if community_id not in SR_feature_dict:
            SR_feature_dict[community_id] = {'unique_mapping':[],'multi_mapping':[]}
        if community_id not in LR_feature_dict:
            LR_feature_dict[community_id] = {'unique_mapping':[],'multi_mapping':[],'5_end_truncation':[],'3_end_truncation':[]}
        feat = np.concatenate([
                    [len(SR_feature_dict[community_id]['unique_mapping']),len(SR_feature_dict[community_id]['multi_mapping']),len(SR_feature_dict[community_id]['unique_mapping']) + len(SR_feature_dict[community_id]['multi_mapping'])],
                    [len(LR_feature_dict[community_id]['unique_mapping']),len(LR_feature_dict[community_id]['multi_mapping']),len(LR_feature_dict[community_id]['unique_mapping']) + len(LR_feature_dict[community_id]['multi_mapping'])],
                    get_stats(SR_feature_dict[community_id]['unique_mapping']),
                    get_stats(SR_feature_dict[community_id]['multi_mapping']),
                    get_stats(LR_feature_dict[community_id]['unique_mapping']),
                    get_stats(LR_feature_dict[community_id]['multi_mapping']),
                    get_stats(LR_feature_dict[community_id]['5_end_truncation']),
                    get_stats(LR_feature_dict[community_id]['3_end_truncation'])])
        features.append(feat)
        community_ids.append(community_id)
    features = np.array(features)
    np.savez_compressed(f'{output_path}/temp/machine_learning/features.npz', features=features,community_ids=np.array(community_ids))
    with open(f'{output_path}/temp/machine_learning/gene_community_id_dict.pkl','wb') as f:
        pickle.dump(gene_community_id_dict,f)