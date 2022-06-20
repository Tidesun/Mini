import os 
import numpy as np
import torch
from torch import nn
from torch.nn.utils.rnn import pad_sequence,pack_sequence,pack_padded_sequence,pad_packed_sequence
torch.set_default_dtype(torch.float64)
# from qpth.qp import QPFunction,QPSolvers
import warnings
warnings.filterwarnings('ignore')
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
class AbNET(nn.Module):

    def __init__(self, embedding_dim, hidden_dim_1,hidden_dim_2,n_output,n_replicates):
        super(AbNET, self).__init__()
        self.n_replicates = n_replicates
        self.lstm = nn.LSTM(embedding_dim, hidden_dim_1,bidirectional=True,num_layers=2).to(device)
        self.prelu = nn.LeakyReLU().to(device)
        self.batchnorm1 = nn.BatchNorm1d(hidden_dim_1+5, affine=False).to(device)
        self.batchnorm2 = nn.BatchNorm1d(hidden_dim_2, affine=False).to(device)
        self.fc1 = nn.Linear(hidden_dim_1+5,hidden_dim_2).to(device)
        self.fc2 = nn.Linear(hidden_dim_2,n_output).to(device)
        self.sigmoid = nn.Sigmoid().to(device)
        self.softmax = nn.Softmax(dim=1).to(device)
        self.init_weight()
    def init_weight(self):
        for name, param in self.lstm.named_parameters():
            if 'bias' in name:
                nn.init.constant(param, 0.0)
            elif 'weight' in name:
                nn.init.kaiming_uniform_(param)
        for fc in [self.fc1,self.fc2]:
            for name, param in fc.named_parameters():
                if 'bias' in name:
                    nn.init.constant(param, 0.0)
                elif 'weight' in name:
                    nn.init.kaiming_uniform_(param)
    def forward(self,all_packed_x,all_fixed_x,all_matrics):

        x = self.lstm(all_packed_x.double())
        x = x[1][0][-1,:].detach()
        assert not torch.isnan(x).any()
        x = torch.cat([x,all_fixed_x.double()],dim=1)
        assert not torch.isnan(x).any()
        x = self.batchnorm1(x)
        assert not torch.isnan(x).any()
        x = self.prelu(x)
        assert not torch.isnan(x).any()
        x = self.fc1(x)
#         print(x)
        assert not torch.isnan(x).any()
        x = self.batchnorm2(x)
        assert not torch.isnan(x).any()
        x = self.prelu(x)
        assert not torch.isnan(x).any()
        x= self.fc2(x)
        assert not torch.isnan(x).any()
#         x = self.sigmoid(x)
        x = self.softmax(x)
        assert not torch.isnan(x).any()
        x = x.view(x.shape[0]//self.n_replicates,self.n_replicates,2)
        x = torch.transpose(x,1,2)
        packed = pack_sequence(all_matrics,enforce_sorted=False).to(device)
        seq_unpacked, lens_unpacked = pad_packed_sequence(packed, batch_first=True)
        output = ((seq_unpacked[:,:,0,:] * x[:,[1],:]) + (seq_unpacked[:,:,1,:] * x[:,[0],:]))
#         list_of_output = [output[i,:lens_unpacked[i].item(),:] for i in range(lens_unpacked.shape[0])]
        return output,x[:,0,:]
def get_generalized_condition_number(A):
    cpu_A = A.cpu()
    _,s,_ = torch.linalg.svd(cpu_A)
    try:
        res = s[0]/s[s>0][-1]
    except Exception as e:
        print(cpu_A)
        print(e)
        return torch.Tensor([1])[0]
    return res
def get_features(all_variable_features,all_fixed_features,all_matrics):
    packed = pack_sequence(all_variable_features,enforce_sorted=False).to(device)
    seq_unpacked, lens_unpacked = pad_packed_sequence(packed, batch_first=True)
    normalized = nn.BatchNorm1d(seq_unpacked.shape[1]).to(device)(seq_unpacked)
    normalized_packed = pack_padded_sequence(normalized,lens_unpacked,batch_first=True,enforce_sorted=False)
    return normalized_packed,torch.stack(all_fixed_features),all_matrics
def get_batch_data(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list,all_num_exons_list,all_isoform_length_list):
    all_variable_features = []
    all_fixed_features = []
    all_matrics = []
    # for each gene
    for sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM,num_exons,isoform_length in zip(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list,all_num_exons_list,all_isoform_length_list):
        num_exons_per_transcript = torch.unsqueeze(num_exons,1).to(device)
        isoform_length_per_transcript = torch.unsqueeze(isoform_length,1).to(device)
        # for each replicate of gene
        multi_replicate_tpm = []
        for sr_A,lr_A,sr_b,lr_b,sr_tpm,lr_tpm,i in zip(sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM,range(len(sr_A_list))):
            sr_tensor = torch.stack([torch.flatten(sr_A),sr_b.repeat_interleave(sr_A.shape[1])],1).to(device)
            lr_tensor = torch.stack([torch.flatten(lr_A),lr_b.repeat_interleave(lr_A.shape[1])],1).to(device)
            sr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(sr_A,dim=0),1).to(device)
            lr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(lr_A,dim=0),1).to(device)
            coverage = torch.stack([sr_b.sum(),lr_b.sum()]).to(device)
            condition_number = torch.stack([get_generalized_condition_number(sr_A),get_generalized_condition_number(lr_A)]).to(device)
            Ab_tensor = torch.cat([sr_tensor,lr_tensor],0).to(device)
            transcript_features_tensor_1 = torch.cat([sr_region_per_transcript_tensor,lr_region_per_transcript_tensor],1)
            transcript_features_tensor_2 = torch.cat([num_exons_per_transcript,isoform_length_per_transcript],1)
            tpm_tensor = torch.stack([sr_tpm,lr_tpm],axis=1).to(device)
            num_isoforms = torch.Tensor([sr_A.shape[1]]).to(device)
            feature_tensor_1 = torch.cat([Ab_tensor,transcript_features_tensor_1,transcript_features_tensor_2])
            feature_tensor_2 = torch.cat([condition_number,coverage,num_isoforms])        
            all_variable_features.append(feature_tensor_1)
            all_fixed_features.append(feature_tensor_2)
            multi_replicate_tpm.append(tpm_tensor)
        all_matrics.append(torch.stack(multi_replicate_tpm,-1))
    return get_features(all_variable_features,all_fixed_features,all_matrics)
def predict_params(sr_A,sr_b,lr_A,lr_b,isoform_lengths,isoform_num_exons,model):
    sr_A_tensor = torch.DoubleTensor(sr_A)
    lr_A_tensor = torch.DoubleTensor(lr_A)
    sr_b_tensor = torch.DoubleTensor(sr_b)
    lr_b_tensor = torch.DoubleTensor(lr_b)
    sr_TPM_tensor = torch.ones(sr_A.shape[1]).double()
    lr_TPM_tensor = torch.ones(lr_A.shape[1]).double()
    all_packed_x,all_fixed_x,all_matrics = get_batch_data([[sr_A_tensor]],[[lr_A_tensor]],[[sr_b_tensor]],[[lr_b_tensor]],[[sr_TPM_tensor]],[[lr_TPM_tensor]],[torch.FloatTensor(isoform_num_exons)],[torch.FloatTensor(isoform_lengths)])
    try:
        _,params = model.forward(all_packed_x,all_fixed_x,all_matrics)
        alpha = params.item()
    except Exception as e:
        print(e)
        alpha = 0.5
    return alpha
def load_model(model_path):
    rnn = AbNET(2,64,32,2,1).to(device)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    checkpoint = torch.load(model_path,map_location='cpu')
    # checkpoint = torch.load("{}/models/{}".format(dir_path,model_path),map_location='cpu')
    rnn.load_state_dict(checkpoint['model_state_dict'])
    rnn.eval()
    # rnn.train(True)
    return rnn