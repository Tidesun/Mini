from libraries.GTFBasics import GTFFile
import subprocess
import os
import numpy as np
def get_SR_expression_arr(expr_dict,short_read_gene_matrix_dict,gene_isoforms_length_dict):
    gene_isoform_expression_dict = {}
    for chr_name in short_read_gene_matrix_dict:
        gene_isoform_expression_dict[chr_name] = {}
        for gene_name in short_read_gene_matrix_dict[chr_name]:
            isoform_names_indics = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']
            SR_expression_arr = np.zeros((len(isoform_names_indics)))
            SR_expected_counts_arr = np.zeros((len(isoform_names_indics)))
            
            for isoform_name in isoform_names_indics:
                if isoform_name in expr_dict:
                    isoform_index = isoform_names_indics[isoform_name]
                    SR_expression_arr[isoform_index] = expr_dict[isoform_name]
                    SR_expected_counts_arr[isoform_index] = gene_isoforms_length_dict[chr_name][gene_name][isoform_name] * expr_dict[isoform_name]
            gene_isoform_expression_dict[chr_name][gene_name] = {}
            gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts'] = SR_expected_counts_arr
            gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] = SR_expression_arr
    return gene_isoform_expression_dict
def SR_external_quantification(short_read_gene_matrix_dict,gene_isoforms_length_dict,SR_quantification_option,SR_fastq_list,SR_read_len,ref_annotation,ref_genome,output_dir,threads):
    external_bin_path = os.path.dirname(os.path.realpath(__file__))+'/external_bin'
    if SR_quantification_option == 'Kallisto':
        gtf_file = GTFFile(ref_annotation,ref_genome)
        ref_transcriptome = f'{output_dir}/ref_transcriptome.fa'
        gtf_file.write_fa(ref_transcriptome)
        subprocess.run([f'{external_bin_path}/kallisto', "index",ref_transcriptome,'-i',f'{output_dir}/kallisto.index'])
        if len(SR_fastq_list) == 1:
            subprocess.run([f'{external_bin_path}/kallisto', "quant",'--single','-l',str(SR_read_len),'-s','0.01','-i',f'{output_dir}/kallisto.index','-o',f'{output_dir}/kallisto','-t',str(threads),\
                '--gtf',ref_annotation,SR_fastq_list[0]])
        else:
            subprocess.run([f'{external_bin_path}/kallisto', "quant",'-i',f'{output_dir}/kallisto.index','-o',f'{output_dir}/kallisto','-t',str(threads),\
                '--gtf',ref_annotation,SR_fastq_list[0],SR_fastq_list[1]])
        expr_dict = {}
        with open(f'{output_dir}/kallisto/abundance.tsv') as f:
            f.readline()
            for line in f:
                fields = line.split('\t')
                transcript_id,tpm = fields[0],fields[4]
                expr_dict[transcript_id] = float(tpm)
        gene_isoform_expression_dict = get_SR_expression_arr(expr_dict,short_read_gene_matrix_dict,gene_isoforms_length_dict)
        return gene_isoform_expression_dict