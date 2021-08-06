import os 
import numpy as np
import torch
from torch import nn
from torch.nn.utils.rnn import pack_sequence
from qpth.qp import QPFunction,QPSolvers
import warnings
warnings.filterwarnings('ignore')
def map_beta(beta):
    return 1e-9+(1e-5-1e-9) * beta
class AbLSTM(nn.Module):

    def __init__(self, embedding_dim, hidden_dim,n_output):
        super(AbLSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.lstm = nn.LSTM(embedding_dim, hidden_dim,bidirectional=True)
        self.fc = torch.nn.Linear(hidden_dim, n_output)
        self.relu = nn.ReLU()
        self.prelu = nn.PReLU()
        self.sigmoid = nn.Sigmoid()
        self.batchnorm_1 = nn.BatchNorm1d(hidden_dim)
        self.batchnorm_2 = nn.BatchNorm1d(n_output)
        self.qpnet = QPFunction(verbose=False, check_Q_spd=True,solver=QPSolvers.CVXPY)
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True)
        self.double()
    def forward(self,x_packed,matrics):
        x = self.lstm(x_packed.double())[1][0][-1,:].detach()
        x = self.prelu(x)
        x = self.batchnorm_1(x)
        x = self.fc(x)
        x = self.batchnorm_2(x)
        x = self.sigmoid(x)
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
            alpha = x[i,0]
            beta = x[i,1]
            Q = 2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * map_beta(beta) * torch.eye(num_isoforms).cpu()
            Q = Q.double()
            c = -2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T).double()
            lb = np.zeros(num_isoforms)
            e = torch.DoubleTensor().cpu()
            G = - torch.cat((SR_isoform_region_matrix[SR_region_read_count_matrix>0,:], LR_isoform_region_matrix[LR_region_read_count_matrix>0,:]), dim=0).double()
            h = - np.ones(G.shape[0])/(1/P)
            G = G.double().cpu()
            h = torch.DoubleTensor(h).cpu()
            output = self.qpnet(Q, c,G,h,e,e)
            if output.sum() != 0:
                output = output/output.sum()
            outputs.append(output)
        return outputs,x
def get_batch_data(sr_A_list,lr_A_list,sr_b_list,lr_b_list):
    batch = []
    matrics = []
    for sr_A,lr_A,sr_b,lr_b in zip(sr_A_list,lr_A_list,sr_b_list,lr_b_list):
        sr_A_tensor = torch.unsqueeze(sr_A.sum(axis=1),1).cpu()
        lr_A_tensor = torch.unsqueeze(lr_A.sum(axis=1),1).cpu()
#         sr_A_tensor = sr_A.cpu()
#         lr_A_tensor = lr_A.cpu()
        sr_b_tensor = torch.unsqueeze(sr_b,1).cpu()
        lr_b_tensor = torch.unsqueeze(lr_b,1).cpu()
        feature_tensor = torch.cat([torch.cat([sr_A_tensor,lr_A_tensor],0),torch.cat([sr_b_tensor,lr_b_tensor],0)],1)
        batch.append(feature_tensor)
        matrics.append((sr_A,lr_A,sr_b,lr_b))
    packed = pack_sequence(batch,enforce_sorted=False).cpu()
    return packed,matrics
def predict_params(sr_A,sr_b,lr_A,lr_b,model):
    sr_A_tensor = torch.IntTensor(sr_A)
    lr_A_tensor = torch.IntTensor(lr_A)
    sr_b_tensor = nn.functional.normalize(torch.FloatTensor(sr_b),dim=0)
    lr_b_tensor = nn.functional.normalize(torch.FloatTensor(lr_b),dim=0)
    packed,matrics = get_batch_data([sr_A_tensor],[lr_A_tensor],[sr_b_tensor],[lr_b_tensor])
    try:
        _,params = model.forward(packed,matrics)
        alpha = params[0,0].item()
        beta = map_beta(params[0,1]).item()
    except Exception as e:
        print(e)
        alpha = 0.5
        beta = 1e-6
    return alpha,beta
def load_model():
    rnn = AbLSTM(2,16, 2)
    dir_path = os.path.dirname(os.path.realpath(__file__))
    checkpoint = torch.load("{}/models/model.pt".format(dir_path),map_location='cpu')
    rnn.load_state_dict(checkpoint['model_state_dict'])
    rnn.eval()
    return rnn
    