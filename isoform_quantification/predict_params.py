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

    def __init__(self, embedding_dim, hidden_dim_1,hidden_dim_2,n_output):
        super(AbNET, self).__init__()
        self.lstm = nn.LSTM(embedding_dim, hidden_dim_1,bidirectional=True,num_layers=1).to(device)
        self.prelu = nn.LeakyReLU().to(device)
        self.batchnorm1 = nn.BatchNorm1d(hidden_dim_1+5).to(device)
        self.batchnorm2 = nn.BatchNorm1d(hidden_dim_2).to(device)
        self.fc1 = nn.Linear(hidden_dim_1+5,hidden_dim_2).to(device)
        self.fc2 = nn.Linear(hidden_dim_2,n_output).to(device)
        self.sigmoid = nn.Sigmoid().to(device)
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True,solver=QPSolvers.CVXPY)
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True)
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
    def quantify(self,matrics,list_of_alpha):
        outputs = []
        for i,(sr_tpm,lr_tpm),alpha in zip(range(len(matrics)),matrics,list_of_alpha):
            output = (1 - alpha) * sr_tpm + alpha * lr_tpm
#             if torch.isnan(output).any():
#                 print(sr_tpm)
#                 print(lr_tpm)
#                 print(alpha)
#             assert torch.isnan(output).any() == False
            outputs.append(output)
        return outputs
    def forward(self,list_of_x_packed,all_matrics):
        list_of_x,list_of_outputs = [],[]
        for x_packed in list_of_x_packed:
            
            feature_vector_1,feature_vector_2 = x_packed
            x = self.lstm(feature_vector_1)
            x = x[1][0][-1,:].detach()
#             if torch.isnan(x).any():
#                 print(x)
            x = torch.cat([x,feature_vector_2],dim=1)
            x = self.batchnorm1(x)
            x = self.prelu(x)
            x = self.fc1(x)
            x = self.batchnorm2(x)
            x = self.prelu(x)
            x= self.fc2(x)
            x = self.sigmoid(x)
#             assert torch.isnan(x).any() == False
            list_of_x.append(x)
        output_alpha = torch.cat(list_of_x,axis=1)
        # for matrics,alpha in zip(all_matrics,list_of_x):
        #     outputs = self.quantify(matrics,alpha)
        #     list_of_outputs.append(outputs)
        list_of_outputs_tensor = []
        # for i in range(len(list_of_outputs[0])):
        #     outputs = []
        #     for j in range(len(list_of_outputs)):
        #         outputs.append(list_of_outputs[j][i])
        #     list_of_outputs_tensor.append(torch.stack(outputs).T)
        return list_of_outputs_tensor,output_alpha
def get_generalized_condition_number(A):
    cpu_A = A.cpu()
    _,s,_ = torch.linalg.svd(cpu_A)
    return s[0]/s[s>0][-1]

def get_batch_data(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list):
    all_batches_1 = [[] for i in range(len(all_sr_A_list[0]))]
    all_batches_2 = [[] for i in range(len(all_sr_A_list[0]))]
    all_matrics = [[] for i in range(len(all_sr_A_list[0]))]
    # for each gene
    for sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM in zip(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list):
        # for each replicate of gene
        for sr_A,lr_A,sr_b,lr_b,sr_tpm,lr_tpm,i in zip(sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM,range(len(sr_A_list))):
            sr_A_tensor = torch.unsqueeze(torch.count_nonzero(sr_A,dim=1),1).to(device)
            lr_A_tensor = torch.unsqueeze(torch.count_nonzero(lr_A,dim=1),1).to(device)
            sr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(sr_A,dim=0),1).to(device)
            lr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(lr_A,dim=0),1).to(device)
            sr_b_tensor = torch.unsqueeze(sr_b,1).to(device)
            if sr_b_tensor.sum() != 0:
                sr_b_tensor = sr_b_tensor/sr_b_tensor.sum()
            lr_b_tensor = torch.unsqueeze(lr_b,1).to(device)
            if lr_b_tensor.sum() != 0:
                lr_b_tensor = lr_b_tensor/lr_b_tensor.sum()
            coverage = torch.stack([sr_b.sum(),lr_b.sum()]).to(device)
            condition_number = torch.stack([get_generalized_condition_number(sr_A),get_generalized_condition_number(lr_A)]).to(device)
            Ab_tensor = torch.cat([torch.cat([sr_A_tensor,lr_A_tensor],0),torch.cat([sr_b_tensor,lr_b_tensor],0)],1)
            num_region_per_transcripts_tensor = torch.cat([sr_region_per_transcript_tensor,lr_region_per_transcript_tensor],1)
            tpm_tensor = torch.stack([sr_tpm,lr_tpm],axis=1).to(device)
            num_isoforms = torch.Tensor([sr_A_tensor.shape[1]]).to(device)
            feature_tensor_1 = torch.cat([Ab_tensor,num_region_per_transcripts_tensor])
            feature_tensor_2 = torch.cat([condition_number,coverage,num_isoforms])        
            all_batches_1[i].append(feature_tensor_1)
            all_batches_2[i].append(feature_tensor_2)
            all_matrics[i].append((sr_tpm.to(device),lr_tpm.to(device)))
#             all_matrics[i].append((sr_A,lr_A,sr_b,lr_b))
    all_packed = []
    for batch,feature_batch in zip(all_batches_1,all_batches_2):
        packed = pack_sequence(batch,enforce_sorted=False).to(device)
#         seq_unpacked, lens_unpacked = pad_packed_sequence(packed, batch_first=True)
#         normalized = nn.BatchNorm1d(seq_unpacked.shape[1]).to(device)(seq_unpacked)
#         normalized_packed = pack_padded_sequence(normalized,lens_unpacked,batch_first=True,enforce_sorted=False)
        all_packed.append((packed.double(),torch.stack(feature_batch).double()))
#         all_packed.append(packed)
    return all_packed,all_matrics
def predict_params(sr_A,sr_b,lr_A,lr_b,model):
    sr_A_tensor = torch.DoubleTensor(sr_A)
    lr_A_tensor = torch.DoubleTensor(lr_A)
    sr_b_tensor = torch.DoubleTensor(sr_b)
    lr_b_tensor = torch.DoubleTensor(lr_b)
    sr_TPM_tensor = torch.ones(sr_A.shape[1]).double()
    lr_TPM_tensor = torch.ones(lr_A.shape[1]).double()
    packed,matrics = get_batch_data([[sr_A_tensor]],[[lr_A_tensor]],[[sr_b_tensor]],[[lr_b_tensor]],[[sr_TPM_tensor]],[[lr_TPM_tensor]])
    try:
        _,params = model.forward(packed,matrics)
        alpha = params.item()
    except Exception as e:
        print(e)
        alpha = 0.5
    return alpha
def load_model(model_path):
    rnn = AbNET(2,64,32,1).to(device)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    checkpoint = torch.load("{}/models/{}".format(dir_path,model_path),map_location='cpu')
    rnn.load_state_dict(checkpoint['model_state_dict'])
    rnn.eval()
    return rnn