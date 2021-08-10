import os 
import numpy as np
import torch
from torch import nn
from torch.nn.utils.rnn import pack_sequence
from qpth.qp import QPFunction,QPSolvers
import warnings
warnings.filterwarnings('ignore')
class AbLSTM(nn.Module):

    def __init__(self, embedding_dim, hidden_dim,n_output):
        super(AbLSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.lstm = nn.LSTM(embedding_dim, hidden_dim,bidirectional=True,num_layers=1)
        self.fc = torch.nn.Linear(hidden_dim, n_output)
        self.relu = nn.ReLU()
        self.prelu = nn.PReLU()
        self.sigmoid = nn.Sigmoid()
        self.batchnorm_1 = nn.BatchNorm1d(hidden_dim)
        self.batchnorm_2 = nn.BatchNorm1d(n_output)
        self.qpnet = QPFunction(verbose=False, check_Q_spd=True,solver=QPSolvers.CVXPY)
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True)
        self.double()
        self.init_weight()
    def init_weight(self):
        for name, param in self.lstm.named_parameters():
            if 'bias' in name:
                nn.init.constant(param, 0.0)
            elif 'weight' in name:
                nn.init.kaiming_uniform_(param)
    def quantify(self,matrics,x):
        outputs = []
        for i,(SR_isoform_region_matrix,LR_isoform_region_matrix,SR_region_read_count_matrix,LR_region_read_count_matrix) in zip(range(len(matrics)),matrics):
            SR_isoform_region_matrix = SR_isoform_region_matrix.double().cpu()
            LR_isoform_region_matrix = LR_isoform_region_matrix.double().cpu()
            SR_region_read_count_matrix = SR_region_read_count_matrix.double().cpu()
            if SR_region_read_count_matrix.sum() != 0:
                SR_region_read_count_matrix = SR_region_read_count_matrix/SR_region_read_count_matrix.sum()
            LR_region_read_count_matrix = LR_region_read_count_matrix.double().cpu()
            if LR_region_read_count_matrix.sum()!=0:
                LR_region_read_count_matrix = LR_region_read_count_matrix/LR_region_read_count_matrix.sum()
            num_isoforms = SR_isoform_region_matrix.shape[1]
            P = 1e-6
            alpha = x[i]
            beta_1 = 1e-6
#             beta = x[i,1]
#             beta_1 = map_beta(beta)
            Q = 2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta_1 * torch.eye(num_isoforms).cpu()
#             Q = 2 * alpha * torch.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * (1-alpha) * torch.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta_1 * torch.eye(num_isoforms).cpu()
            Q = Q.double()
            G = -torch.eye(num_isoforms).double().cpu()
            h = torch.zeros(num_isoforms).double().cpu()
#             c = -2 * alpha * torch.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * (1-alpha) * torch.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T).double()
            c = -2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T).double()
            e = torch.DoubleTensor().cpu()
            output = self.qpnet(Q, c,G,h,e,e)
            if output.sum() != 0:
                output = output/output.sum()
            outputs.append(output)
        return outputs
    def forward(self,x_packed,matrics):
        x = self.lstm(x_packed.double())[1][0][-1,:].detach()
        x = self.prelu(x)
        x = self.batchnorm_1(x)
        x = self.fc(x)
        x = self.batchnorm_2(x)
        x = self.sigmoid(x)
        outputs = self.quantify(matrics,x)
        return outputs,x
def get_generalized_condition_number(A):
    cpu_A = A.cpu()
    cpu_A = (1/cpu_A)*cpu_A
    cpu_A = torch.nan_to_num(cpu_A,0)
    cond_num = torch.linalg.cond(cpu_A)
    if torch.isinf(cond_num):
        _,s,_ = torch.linalg.svd(cpu_A)
        return s[0]/s[s>0][-1]
    else:
        return cond_num

def get_batch_data(sr_A_list,lr_A_list,sr_b_list,lr_b_list):
    batch = []
    matrics = []
    for sr_A,lr_A,sr_b,lr_b in zip(sr_A_list,lr_A_list,sr_b_list,lr_b_list):
        sr_A_tensor = torch.unsqueeze(sr_A.sum(axis=1),1).cpu()
        lr_A_tensor = torch.unsqueeze(lr_A.sum(axis=1),1).cpu()
        sr_b_tensor = torch.unsqueeze(nn.functional.normalize(sr_b,dim=0),1).cpu()
        lr_b_tensor = torch.unsqueeze(nn.functional.normalize(lr_b,dim=0),1).cpu()
        coverage = torch.stack([sr_b.sum(),lr_b.sum()])
        condition_number = torch.stack([get_generalized_condition_number(sr_A),get_generalized_condition_number(lr_A)]).cpu()
        feature_tensor = torch.cat([torch.cat([sr_A_tensor,lr_A_tensor],0),torch.cat([sr_b_tensor,lr_b_tensor],0)],1)
        feature_tensor = torch.cat([torch.unsqueeze(condition_number,0),torch.unsqueeze(coverage,0),feature_tensor])
        batch.append(feature_tensor)
        matrics.append((sr_A,lr_A,sr_b,lr_b))
    packed = pack_sequence(batch,enforce_sorted=False).cpu()
    return packed,matrics
def predict_params(sr_A,sr_b,lr_A,lr_b,model):
    sr_A_tensor = torch.FloatTensor(sr_A)
    lr_A_tensor = torch.FloatTensor(lr_A)
    sr_b_tensor = torch.FloatTensor(sr_b)
    lr_b_tensor = torch.FloatTensor(lr_b)
    packed,matrics = get_batch_data([sr_A_tensor],[lr_A_tensor],[sr_b_tensor],[lr_b_tensor])
    try:
        _,params = model.forward(packed,matrics)
        alpha = params.item()
        beta = 1e-6
    except Exception as e:
        print(e)
        alpha = 0.5
        beta = 1e-6
    return alpha,beta
def load_model(model_path):
    rnn = AbLSTM(2,512, 1).cpu()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    checkpoint = torch.load("{}/models/{}".format(dir_path,model_path),map_location='cpu')
    rnn.load_state_dict(checkpoint['model_state_dict'])
    rnn.eval()
    return rnn