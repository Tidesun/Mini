import os 
import numpy as np
import torch
from torch import nn
from torch.nn.utils.rnn import pad_sequence,pack_sequence,pack_padded_sequence,pad_packed_sequence
# from qpth.qp import QPFunction,QPSolvers
import warnings
warnings.filterwarnings('ignore')
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
class GetLastLayerOut(nn.Module):
    def __init__(self):
        super(GetLastLayerOut, self).__init__()
    def forward(self, x):
        return x[1][0][-1,:].detach()
class AbLSTM(nn.Module):

    def __init__(self, embedding_dim, hidden_dim,n_output):
        super(AbLSTM, self).__init__()
        self.hidden_dim = hidden_dim
        self.model = nn.Sequential(
        nn.LSTM(embedding_dim, hidden_dim,bidirectional=True,num_layers=3),
        GetLastLayerOut(),
        nn.PReLU(),
        nn.BatchNorm1d(hidden_dim),
        nn.Linear(hidden_dim, n_output),
        nn.BatchNorm1d(n_output),
#         nn.Softmax()
        
        nn.Sigmoid() \
        ).to(device).double()
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True,solver=QPSolvers.CVXPY)
#         self.qpnet = QPFunction(verbose=False, check_Q_spd=True)
#         self.init_weight()
    def init_weight(self):
        for name, param in self.lstm.named_parameters():
            if 'bias' in name:
                nn.init.constant(param, 0.0)
            elif 'weight' in name:
                nn.init.kaiming_uniform_(param)
#     def quantify(self,matrics,x):
#         outputs = []
#         for i,(SR_isoform_region_matrix,LR_isoform_region_matrix,SR_region_read_count_matrix,LR_region_read_count_matrix) in zip(range(len(matrics)),matrics):
#             SR_isoform_region_matrix = SR_isoform_region_matrix.double().cuda()
#             LR_isoform_region_matrix = LR_isoform_region_matrix.double().cuda()
#             SR_region_read_count_matrix = SR_region_read_count_matrix.double().cuda()
#             if SR_region_read_count_matrix.sum() != 0:
#                 SR_region_read_count_matrix = SR_region_read_count_matrix/SR_region_read_count_matrix.sum()
#             LR_region_read_count_matrix = LR_region_read_count_matrix.double().cuda()
#             if LR_region_read_count_matrix.sum()!=0:
#                 LR_region_read_count_matrix = LR_region_read_count_matrix/LR_region_read_count_matrix.sum()
#             num_isoforms = SR_isoform_region_matrix.shape[1]
#             P = 1e-6
#             alpha = x[i]
#             beta_1 = 1e-6
# #             beta = x[i,1]
# #             beta_1 = map_beta(beta)
#             Q = 2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta_1 * torch.eye(num_isoforms).cuda()
# #             Q = 2 * alpha * torch.matmul(SR_isoform_region_matrix.T,SR_isoform_region_matrix) + 2 * (1-alpha) * torch.matmul(LR_isoform_region_matrix.T,LR_isoform_region_matrix) + 2 * beta_1 * torch.eye(num_isoforms).cuda()
#             Q = Q.double()
#             G = -torch.eye(num_isoforms).double().cuda()
#             h = torch.zeros(num_isoforms).double().cuda()
# #             c = -2 * alpha * torch.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * (1-alpha) * torch.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T).double()
#             c = -2 * (1.0 - alpha) * torch.matmul(SR_isoform_region_matrix.T,SR_region_read_count_matrix.T) - 2 * alpha * torch.matmul(LR_isoform_region_matrix.T,LR_region_read_count_matrix.T).double()
#             e = torch.DoubleTensor().cuda()
#             output = self.qpnet(Q, c,G,h,e,e)
#             if output.sum() != 0:
#                 output = output/output.sum()
#             outputs.append(output)
#         return outputs
    def quantify(self,matrics,list_of_alpha):
        outputs = []
        for i,(sr_tpm,lr_tpm),alpha in zip(range(len(matrics)),matrics,list_of_alpha):
#             if alpha > 1-0.01:
#                 output = alpha * lr_tpm
#             elif alpha < 0.01:
#                 output = (1 - alpha) * sr_tpm
#             else:
            output = (1 - alpha) * sr_tpm + alpha * lr_tpm
            assert torch.isnan(output).any() == False
#             if output.sum() != 0:
#                 output = output/output.sum()
            outputs.append(output)
        return outputs
    def forward(self,list_of_x_packed,all_matrics):
        list_of_x,list_of_outputs = [],[]
        for x_packed in list_of_x_packed:
            x = self.model(x_packed.double())
            assert torch.isnan(x).any() == False
            list_of_x.append(x)
        output_alpha = torch.cat(list_of_x,axis=1)
        for matrics,alpha in zip(all_matrics,list_of_x):
            outputs = self.quantify(matrics,alpha)
            list_of_outputs.append(outputs)
        list_of_outputs_tensor = []
        for i in range(len(list_of_outputs[0])):
            outputs = []
            for j in range(len(list_of_outputs)):
                outputs.append(list_of_outputs[j][i])
            list_of_outputs_tensor.append(torch.stack(outputs).T)
        return list_of_outputs_tensor,output_alpha
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

def get_batch_data(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list):
    all_batches = [[] for i in range(len(all_sr_A_list[0]))]
    all_matrics = [[] for i in range(len(all_sr_A_list[0]))]
    # for each gene
    for sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM in zip(all_sr_A_list,all_lr_A_list,all_sr_b_list,all_lr_b_list,all_sr_tpm_list,all_lr_tpm_list):
        # for each replicate of gene
        for sr_A,lr_A,sr_b,lr_b,sr_tpm,lr_tpm,i in zip(sr_A_list,lr_A_list,sr_b_list,lr_b_list,sr_TPM,lr_TPM,range(len(sr_A_list))):
            sr_A_tensor = torch.unsqueeze(torch.count_nonzero(sr_A,dim=1),1).to(device)
            lr_A_tensor = torch.unsqueeze(torch.count_nonzero(lr_A,dim=1),1).to(device)
            sr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(sr_A,dim=0),1).to(device)
            lr_region_per_transcript_tensor = torch.unsqueeze(torch.count_nonzero(lr_A,dim=0),1).to(device)
            sr_b_tensor = torch.unsqueeze(nn.functional.normalize(sr_b,dim=0),1).to(device)
            lr_b_tensor = torch.unsqueeze(nn.functional.normalize(lr_b,dim=0),1).to(device)
            coverage = torch.stack([sr_b.sum(),lr_b.sum()]).to(device)
            condition_number = torch.stack([get_generalized_condition_number(sr_A),get_generalized_condition_number(lr_A)]).to(device)
            feature_tensor = torch.cat([torch.cat([sr_A_tensor,lr_A_tensor],0),torch.cat([sr_b_tensor,lr_b_tensor],0)],1)
            num_region_per_transcripts_tensor = torch.cat([sr_region_per_transcript_tensor,lr_region_per_transcript_tensor],1)
            feature_tensor = torch.cat([torch.unsqueeze(condition_number,0),torch.unsqueeze(coverage,0),num_region_per_transcripts_tensor,feature_tensor])
            all_batches[i].append(feature_tensor)
            all_matrics[i].append((sr_tpm.to(device),lr_tpm.to(device)))
#             all_matrics[i].append((sr_A,lr_A,sr_b,lr_b))
    all_packed = []
    for batch in all_batches:
        packed = pack_sequence(batch,enforce_sorted=False).to(device)
        seq_unpacked, lens_unpacked = pad_packed_sequence(packed, batch_first=True)
        normalized = nn.BatchNorm1d(seq_unpacked.shape[1]).to(device)(seq_unpacked)
        normalized_packed = pack_padded_sequence(normalized,lens_unpacked,batch_first=True,enforce_sorted=False)
        all_packed.append(normalized_packed)
    return all_packed,all_matrics
def predict_params(sr_A,sr_b,lr_A,lr_b,model):
    sr_A_tensor = torch.FloatTensor(sr_A)
    lr_A_tensor = torch.FloatTensor(lr_A)
    sr_b_tensor = torch.FloatTensor(sr_b)
    lr_b_tensor = torch.FloatTensor(lr_b)
    dump_sr_tpm = torch.ones(sr_A_tensor.shape[1])
    dump_lr_tpm = torch.ones(lr_A_tensor.shape[1])
    packed,matrics = get_batch_data([[sr_A_tensor]],[[lr_A_tensor]],[[sr_b_tensor]],[[lr_b_tensor]],[[dump_sr_tpm]],[[dump_lr_tpm]])
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
    rnn = AbLSTM(2,256, 1).cpu()
    dir_path = os.path.dirname(os.path.realpath(__file__))
    checkpoint = torch.load("{}/models/{}".format(dir_path,model_path),map_location='cpu')
    rnn.load_state_dict(checkpoint['model_state_dict'])
    rnn.eval()
    return rnn