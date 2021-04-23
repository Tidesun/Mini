from pathlib import Path
def generate_TrEESR_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    with open(output_path+"/kvalues.out",'w') as f:
        f.write('Gene\tChr\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\n')
        for chr_name in short_read_gene_matrix_dict:
            for gene_name in short_read_gene_matrix_dict[chr_name]:
                SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                f.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(gene_name,chr_name,SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
def generate_TransELS_output(output_path,short_read_gene_matrix_dict,long_read_gene_matrix_dict,gene_isoform_tpm_expression_dict,raw_isoform_exons_dict,gene_isoforms_length_dict):
    Path(output_path).mkdir(parents=True, exist_ok=True)
    with open(output_path+"/expression_gene.out",'w') as f_gene:
        with open(output_path+"/expression_isoform.out",'w') as f_isoform:
            f_gene.write('Gene\tChr\tTPM\n')
            f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\n')
            # f_isoform.write('Isoform\tGene\tChr\tStart\tEnd\tIsoform_length\tTPM\tSR_k_value\tSR_regular_condition_number\tSR_generalized_condition_number\tLR_k_value\tLR_regular_condition_number\tLR_generalized_condition_number\n')
            for chr_name in short_read_gene_matrix_dict:
                for gene_name in short_read_gene_matrix_dict[chr_name]:
                    tpm_sum = 0
                    for isoform_name in short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']:
                        start_pos = min(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['start_pos'])
                        end_pos = max(raw_isoform_exons_dict[chr_name][gene_name][isoform_name]['end_pos'])
                        isoform_len = gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
                        isoform_index = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics'][isoform_name]
                        tpm = gene_isoform_tpm_expression_dict[chr_name][gene_name]['tpm'][isoform_index]
                        tpm_sum += tpm

                        SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number = short_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                        LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number = long_read_gene_matrix_dict[chr_name][gene_name]['condition_number']
                        f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(isoform_name,gene_name,chr_name,start_pos,end_pos,isoform_len,tpm))

                        # f_isoform.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(SR_kvalue,SR_regular_condition_number,SR_generalized_condition_number,LR_kvalue,LR_regular_condition_number,LR_generalized_condition_number))
                    f_gene.write('%s\t%s\t%f\n'%(gene_name,chr_name,tpm_sum))
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