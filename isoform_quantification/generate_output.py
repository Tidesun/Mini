from pathlib import Path
import numpy as np
import dill as pickle
import io
def generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,info_dict_list,same_structure_isoform_dict,removed_gene_isoform_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    [raw_gene_num_exon_dict,gene_num_exon_dict,gene_num_isoform_dict,raw_isoform_num_exon_dict,isoform_length_dict,num_isoforms_dict] = info_dict_list
    out_dict = short_read_gene_matrix_dict.copy()
    bio = io.BytesIO()
    for chr in out_dict:
        for gene in out_dict[chr]:
            bio.write(str.encode('{}\n'.format(gene)))
            np.savetxt(bio, out_dict[chr][gene]['isoform_region_matrix'],fmt='%.d',delimiter=',')
    
    mystr = bio.getvalue().decode('latin1')
    with open(output_path+'/sr_A.out','w') as f:
        f.write(mystr)
    bio = io.BytesIO()
    for chr in long_read_gene_matrix_dict:
        for gene in long_read_gene_matrix_dict[chr]:
            bio.write(str.encode('{}\n'.format(gene)))
            np.savetxt(bio, long_read_gene_matrix_dict[chr][gene]['isoform_region_matrix'],fmt='%.d',delimiter=',')
    
    mystr = bio.getvalue().decode('latin1')
    with open(output_path+'/lr_A.out','w') as f:
        f.write(mystr)
    list_of_all_genes_chrs = []
    for chr_name in long_read_gene_matrix_dict:
        if chr_name in short_read_gene_matrix_dict:
            for gene_name in long_read_gene_matrix_dict[chr_name]:
                if gene_name in short_read_gene_matrix_dict[chr_name]:
                    list_of_all_genes_chrs.append((gene_name,chr_name))
    
    with open(output_path+"/kvalues_gene.out",'w') as f:
        f.write('Gene\tChr\tNum_isoforms\tNum_exons\tNum_split_exons\tSR_singular_value_product\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tSR_A_dim\tLR_singular_value_product\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\tLR_A_dim\n')
        for (gene_name,chr_name) in list_of_all_genes_chrs:
        # for chr_name in short_read_gene_matrix_dict:
        #     for gene_name in short_read_gene_matrix_dict[chr_name]:
            num_isoforms,num_exons,num_split_exons = gene_num_isoform_dict[chr_name][gene_name],raw_gene_num_exon_dict[chr_name][gene_name],gene_num_exon_dict[chr_name][gene_name]
            SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_singular_value_product = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
            LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_singular_value_product = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
            SR_A_dim = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix'].shape
            LR_A_dim = long_read_gene_matrix_dict[chr_name][gene_name]['isoform_region_matrix'].shape
            f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name,chr_name,num_isoforms,num_exons,num_split_exons,SR_singular_value_product,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_A_dim,LR_singular_value_product,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_A_dim))
        for chr_name in removed_gene_isoform_dict:
            for gene_name in removed_gene_isoform_dict[chr_name]:
                info_dict = removed_gene_isoform_dict[chr_name][gene_name]['info']
                num_isoforms,num_exons,num_split_exons = info_dict['num_isoforms'],info_dict['num_exons'],info_dict['num_split_exons']
                SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_singular_value_product = 'NA','NA','NA','NA'
                LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_singular_value_product = 'NA','NA','NA','NA'
                SR_A_dim = 'NA'
                LR_A_dim = 'NA'
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name,chr_name,num_isoforms,num_exons,num_split_exons,SR_singular_value_product,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_A_dim,LR_singular_value_product,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_A_dim))
    with open(output_path+"/kvalues_isoform.out",'w') as f:
        f.write('Isoform\tGene\tChr\tNum_exons\tIsoform_length\tNum_isoforms\tSR_singular_value_product\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_singular_value_product\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\n')
        for (gene_name,chr_name) in list_of_all_genes_chrs:
        # for chr_name in short_read_gene_matrix_dict:
        #     for gene_name in short_read_gene_matrix_dict[chr_name]:
            for isoform_name in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']:
                num_exons,isoform_length,num_isoforms = raw_isoform_num_exon_dict[isoform_name],isoform_length_dict[isoform_name],num_isoforms_dict[isoform_name]
                SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_singular_value_product = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_singular_value_product = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,num_exons,isoform_length,num_isoforms,SR_singular_value_product,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_singular_value_product,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
        for chr_name in removed_gene_isoform_dict:
            for gene_name in removed_gene_isoform_dict[chr_name]:
                isoform_info_dict = removed_gene_isoform_dict[chr_name][gene_name]['isoform_info']
                for isoform_name in isoform_info_dict:
                    num_exons,isoform_length,num_isoforms = isoform_info_dict[isoform_name]['num_exons'],isoform_info_dict[isoform_name]['isoform_length'],removed_gene_isoform_dict[chr_name][gene_name]['info']['num_isoforms']
                    SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,SR_singular_value_product = 'NA','NA','NA','NA'
                    LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number,LR_singular_value_product = 'NA','NA','NA','NA'
                    f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,num_exons,isoform_length,num_isoforms,SR_singular_value_product,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_singular_value_product,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
def generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,list_of_all_genes_chrs,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict,same_structure_isoform_dict,removed_gene_isoform_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    # with open(output_path+'lr.pkl','wb') as f:
    #     pickle.dump(long_read_gene_matrix_dict,f)
    with open(output_path+"/expression_gene.out",'w') as f_gene:
        with open(output_path+"/expression_isoform.out",'w') as f_isoform:
            f_gene.write('Gene\tChr\tTPM\tSR_expected_counts\tLR_expected_counts\n')
            f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tSR_expected_counts\tLR_expected_counts\tAlpha\n')
            # f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\n')
            for gene_name,chr_name in list_of_all_genes_chrs:
                tpm_sum = 0
                sr_expected_counts_sum = 0
                lr_expected_counts_sum = 0
                for isoform_name in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']:
                    start_pos = min(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                    end_pos = max(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'])
                    isoform_len = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                    isoform_index = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics'][isoform_name]
                    tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'][isoform_index]
                    sr_expected_counts = gene_isoform_tpm_expression_dict[chr_name][gene_name]['SR_expected_counts'][isoform_index]
                    lr_expected_counts = gene_isoform_tpm_expression_dict[chr_name][gene_name]['LR_expected_counts'][isoform_index]
                    alpha = gene_isoform_tpm_expression_dict[chr_name][gene_name]['alpha']
                    tpm_sum += tpm
                    sr_expected_counts_sum += sr_expected_counts
                    lr_expected_counts_sum += lr_expected_counts
                    if chr_name in same_structure_isoform_dict:
                        if gene_name in same_structure_isoform_dict[chr_name]:
                            if isoform_name in same_structure_isoform_dict[chr_name][gene_name]:
                                num_same_structure_isoforms = len(same_structure_isoform_dict[chr_name][gene_name][isoform_name])
                                tpm = tpm/(num_same_structure_isoforms+1)
                                sr_expected_counts = sr_expected_counts/(num_same_structure_isoforms+1)
                                lr_expected_counts = lr_expected_counts/(num_same_structure_isoforms+1)
                                for same_structure_isoform in same_structure_isoform_dict[chr_name][gene_name][isoform_name]:
                                    f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(same_structure_isoform,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,sr_expected_counts,lr_expected_counts,alpha))
                    f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,sr_expected_counts,lr_expected_counts,alpha))
                    # SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                    # LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']

                    # f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
                f_gene.write('{}\t{}\t{}\t{}\t{}\n'.format(gene_name,chr_name,tpm_sum,sr_expected_counts_sum,lr_expected_counts_sum))
                
            # censored gene output  
            for chr_name in removed_gene_isoform_dict:
                for gene_name in removed_gene_isoform_dict[chr_name]:
                    f_gene.write('{}\t{}\t0\t0\t0\n'.format(gene_name,chr_name))
                    for isoform_name in removed_gene_isoform_dict[chr_name][gene_name]['isoforms_length_dict']:
                        start_pos = min(removed_gene_isoform_dict[chr_name][gene_name]['raw_isoform_exons_dict'][isoform_name]['start_pos'])
                        end_pos = max(removed_gene_isoform_dict[chr_name][gene_name]['raw_isoform_exons_dict'][isoform_name]['end_pos'])
                        isoform_len = removed_gene_isoform_dict[chr_name][gene_name]['isoforms_length_dict'][isoform_name]
                        tpm = 0
                        sr_expected_counts = 0
                        lr_expected_counts = 0
                        f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm,sr_expected_counts,lr_expected_counts))
            
    # with open(output_path+"/expression_isoform_result1.out",'w') as f:
    #     f.write('Isoform\tGene\tChr\tTPM\tTPM_SR\tTPM_LR\n')
    #     for chr_name in short_read_gene_matrix_dict:
    #             for gene_name in short_read_gene_matrix_dict[chr_name]:
    #                 for isoform_name in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']:
    #                     isoform_index = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics'][isoform_name]
    #                     tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_by_gene']
    #                     perfect_tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['perfect_tpm_by_gene']
    #                     tpm_sr = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_sr']
    #                     tpm_lr = gene_isoform_tpm_expression_dict[chr_name][gene_name][isoform_index]['tpm_lr']