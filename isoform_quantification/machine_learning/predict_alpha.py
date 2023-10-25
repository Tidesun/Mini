from machine_learning.extract_features import extract_features
from machine_learning.concat_features import concat_features
import config
import pickle
import pandas as pd
import time
import datetime
def predict_alpha(output_path):
    print('Using pretrained model in {} to predict the alpha...'.format(config.pretrained_model_path))
    start_time = time.time()
    extract_features(output_path)
    all_features_arr,all_community_ids,all_alphas = concat_features(output_path)
    # load models
    scaler_path = config.pretrained_model_path+'/preprocessor.pkl'
    model_path = config.pretrained_model_path+'/clf.pkl'
    with open(scaler_path,'rb') as f:
        [scaler,imputer] = pickle.load(f)
    with open(model_path,'rb') as f:
        clf = pickle.load(f)
    X_arr = imputer.transform(scaler.transform(all_features_arr))
    y_pred = clf.predict(X_arr)
    error_alpha = pd.DataFrame({'id':all_community_ids,'error':y_pred,'alpha':all_alphas})
    true_strategy_df = error_alpha.loc[error_alpha.groupby('id').apply(lambda rows:rows.error.idxmin())]
    true_strategy_df.to_csv(config.alpha_df_path,sep='\t',index=False)
    end_time = time.time()
    print('Done in {} seconds at {}'.format(end_time-start_time,str(datetime.datetime.now())),flush=True)
    