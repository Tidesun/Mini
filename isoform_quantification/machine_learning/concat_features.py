import pandas as pd
import numpy as np
import pickle
import scipy.stats
def get_stats(arr):
    if len(arr) == 0:
        return np.array([float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')])
    s = scipy.stats.describe(arr)
    return np.array([s.minmax[0],s.minmax[1],s.mean,s.variance,s.skewness,s.kurtosis])
def concat_features(output_path):
    with open(f'{output_path}/temp/machine_learning/gene_community_id_dict.pkl','rb') as f:
        gene_community_id_dict = pickle.load(f)
    com_df = pd.Series(gene_community_id_dict)
    gene_structure_features_dict = {}
    with open(f'{output_path}/temp/machine_learning/gene_features.pkl','rb') as f:
        features = pickle.load(f)
        for feature in features:
            gene_structure_features_dict[feature[0]] = feature[1:]
    structure_features = {}
    for gene,feature in gene_structure_features_dict.items():
        community_id = com_df.loc[gene]
        [gene_length,exon_length,isoform_length,region_length] = feature
        if community_id not in structure_features:
            structure_features[community_id] = [[],[],[],[]]
        structure_features[community_id][0].append(gene_length)
        structure_features[community_id][1] += isoform_length
        structure_features[community_id][2] += exon_length
        structure_features[community_id][3] += region_length
    features = []
    community_ids = []
    for community_id in sorted(structure_features.keys()):
        feat = np.concatenate([
            [len(structure_features[community_id][0]),len(structure_features[community_id][1]),len(structure_features[community_id][2]),len(structure_features[community_id][3])],
            get_stats(structure_features[community_id][0]),
            get_stats(structure_features[community_id][1]),
            get_stats(structure_features[community_id][2]),
            get_stats(structure_features[community_id][3])])
        features.append(feat)
        community_ids.append(community_id)
    gene_structure_feature_df = pd.DataFrame(np.concatenate([np.array([community_ids]).T,features],axis=1))
    gene_structure_feature_df = gene_structure_feature_df.set_index(0)
    # data feature
    features = np.load(f'{output_path}/temp/machine_learning/features.npz')['features']
    community_ids =np.load(f'{output_path}/temp/machine_learning/features.npz')['community_ids']
    data_feature_df = pd.DataFrame(np.concatenate([np.array([community_ids]).T,features],axis=1))
    data_feature_df = data_feature_df.set_index(0)
    # k value feature
    with open(f'{output_path}/temp/machine_learning/SR_kvalue_dict.pkl','rb') as f:
        SR_kval_df = pickle.load(f)
    Kval_features = {}
    for gene in com_df.index:
        community_id = com_df.loc[gene]
        if community_id not in Kval_features:
            Kval_features[community_id] = [[],[]]
        if gene in SR_kval_df:
            Kval_features[community_id][0].append(SR_kval_df.loc[gene])
        else:
            Kval_features[community_id][0].append(np.float('nan'))
    features = []
    community_ids = []
    for community_id in sorted(Kval_features.keys()):
        feat = get_stats(Kval_features[community_id][0])
        features.append(feat)
        community_ids.append(community_id)
    Kvalue_feature_df = pd.DataFrame(np.concatenate([np.array([community_ids]).T,features],axis=1))
    Kvalue_feature_df = Kvalue_feature_df.set_index(0)
    all_features_df = gene_structure_feature_df.join(data_feature_df,lsuffix='_struc',rsuffix='_data').join(Kvalue_feature_df)
    # predict alpha
    # all_features = []
    # all_community_ids = []
    # all_alphas = []
    # for i in range(0,101,5):
    #     alpha = str(i/100)
    #     if alpha == "0.0" or alpha == '1.0':
    #         alpha = alpha.split('.')[0]
    #     # features
    #     df = all_features_df.copy()
    #     df['alpha'] = i/100
    #     all_features.append(df.values)
    #     all_community_ids.append(df.index)
    #     all_alphas.append(df['alpha'])
    all_features_arr = all_features_df.values
    all_community_ids = all_features_df.index
    return all_features_arr,all_community_ids