from libraries.GTFBasics import GTFFile
import subprocess
import os
import numpy as np
from pathlib import Path
def get_SR_expression_arr(expr_dict,short_read_gene_matrix_dict,gene_isoforms_length_dict):
    gene_isoform_expression_dict = {}
    expr_sum = 0
    for chr_name in short_read_gene_matrix_dict:
        gene_isoform_expression_dict[chr_name] = {}
        for gene_name in short_read_gene_matrix_dict[chr_name]:
            isoform_names_indics = short_read_gene_matrix_dict[chr_name][gene_name]['isoform_names_indics']
            SR_expression_arr = np.zeros((len(isoform_names_indics)))
            SR_expected_counts_arr = np.zeros((len(isoform_names_indics)))
            
            for isoform_name in isoform_names_indics:
                if isoform_name in expr_dict:
                    isoform_index = isoform_names_indics[isoform_name]
                    SR_expected_counts_arr[isoform_index] = expr_dict[isoform_name]
                    SR_expression_arr[isoform_index] = expr_dict[isoform_name] / gene_isoforms_length_dict[chr_name][gene_name][isoform_name]
            gene_isoform_expression_dict[chr_name][gene_name] = {}
            gene_isoform_expression_dict[chr_name][gene_name]['SR_expected_counts'] = SR_expected_counts_arr
            gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] = SR_expression_arr
            expr_sum += SR_expression_arr.sum()
    if expr_sum != 0:
        for chr_name in gene_isoform_expression_dict:
            for gene_name in gene_isoform_expression_dict[chr_name]:
                gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression'] = (gene_isoform_expression_dict[chr_name][gene_name]['SR_isoform_expression']/expr_sum)*1e6
    return gene_isoform_expression_dict
def quantify_by_kallisto(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,SR_read_len,threads):
    gtf_file = GTFFile(ref_annotation,ref_genome)
    ref_transcriptome = f'{temp_dir}/ref_transcriptome.fa'
    gtf_file.write_fa(ref_transcriptome)
    subprocess.run([f'{external_bin_path}/kallisto', "index",ref_transcriptome,'-i',f'{temp_dir}/kallisto.index'])
    if len(SR_fastq_list) == 1:
        subprocess.run([f'{external_bin_path}/kallisto', "quant",'--single','-l',str(SR_read_len),'-s','0.01','-i',f'{temp_dir}/kallisto.index','-o',f'{temp_dir}/kallisto','-t',str(threads),\
            '--gtf',ref_annotation,SR_fastq_list[0]])
    else:
        subprocess.run([f'{external_bin_path}/kallisto', "quant",'-i',f'{temp_dir}/kallisto.index','-o',f'{temp_dir}/kallisto','-t',str(threads),\
            '--gtf',ref_annotation,SR_fastq_list[0],SR_fastq_list[1]])
    expr_dict = {}
    with open(f'{temp_dir}/kallisto/abundance.tsv') as f:
        f.readline()
        for line in f:
            fields = line.split('\t')
            transcript_id,read_counts = fields[0],fields[3]
            expr_dict[transcript_id] = float(read_counts)
    return expr_dict
def quantify_by_salmon(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,threads):
    gtf_file = GTFFile(ref_annotation,ref_genome)
    ref_transcriptome = f'{temp_dir}/ref_transcriptome.fa'
    gtf_file.write_fa(ref_transcriptome)
    subprocess.run([f'{external_bin_path}/salmon-latest_linux_x86_64/bin/salmon', "index",'-p',str(threads),'--gencode',"-t",ref_transcriptome,'-i',f'{temp_dir}/salmon.index'])
    if len(SR_fastq_list) == 1:
        subprocess.run([f'{external_bin_path}/salmon-latest_linux_x86_64/bin/salmon', "quant",'-i',f'{temp_dir}/salmon.index','-o',f'{temp_dir}/salmon','-p',str(threads),'-l','A',\
            '-r',SR_fastq_list[0]])
    else:
        subprocess.run([f'{external_bin_path}/salmon-latest_linux_x86_64/bin/salmon', "quant",'-i',f'{temp_dir}/salmon.index','-o',f'{temp_dir}/salmon','-p',str(threads),'-l','A',\
            '-1',SR_fastq_list[0],'-2',SR_fastq_list[1]])
    expr_dict = {}
    with open(f'{temp_dir}/salmon/quant.sf') as f:
        f.readline()
        for line in f:
            fields = line.split('\t')
            transcript_id,read_counts = fields[0],fields[4]
            expr_dict[transcript_id] = float(read_counts)
    return expr_dict
def quantify_by_rsem(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,threads):
    # gtf_file = GTFFile(ref_annotation,ref_genome)
    # ref_transcriptome = f'{temp_dir}/ref_transcriptome.fa'
    # gtf_file.write_fa(ref_transcriptome)
    subprocess.run([f'{external_bin_path}/RSEM-1.3.3/bin/rsem-prepare-reference', "--gtf",ref_annotation,'-p',str(threads),'--hisat2-hca','--hisat2-path',f'{external_bin_path}/hisat2/',ref_genome,'rsem_index'],cwd=temp_dir)
    if 'fastq' in SR_fastq_list[0] or 'fq' in SR_fastq_list[0]:
        if len(SR_fastq_list) == 1:
            subprocess.run([f'{external_bin_path}/RSEM-1.3.3/bin/rsem-calculate-expression','-p',str(threads),'--hisat2-hca','--hisat2-path',f'{external_bin_path}/hisat2/',\
                SR_fastq_list[0],'rsem_index','rsem'],cwd=temp_dir)
        else:
            subprocess.run([f'{external_bin_path}/RSEM-1.3.3/bin/rsem-calculate-expression','-p',str(threads),'--hisat2-hca','--hisat2-path',f'{external_bin_path}/hisat2/',\
                '--paired-end',SR_fastq_list[0],SR_fastq_list[1],'rsem_index','rsem'],cwd=temp_dir)
    else:
        if len(SR_fastq_list) == 1:
            subprocess.run([f'{external_bin_path}/RSEM-1.3.3/bin/rsem-calculate-expression','-p',str(threads),'--no-qualities','--hisat2-hca','--hisat2-path',f'{external_bin_path}/hisat2/',\
                SR_fastq_list[0],'rsem_index','rsem'],cwd=temp_dir)
        else:
            subprocess.run([f'{external_bin_path}/RSEM-1.3.3/bin/rsem-calculate-expression','-p',str(threads),'--no-qualities','--hisat2-hca','--hisat2-path',f'{external_bin_path}/hisat2/',\
                '--paired-end',SR_fastq_list[0],SR_fastq_list[1],'rsem_index','rsem'],cwd=temp_dir)
    expr_dict = {}
    with open(f'{temp_dir}/rsem.isoforms.results') as f:
        f.readline()
        for line in f:
            fields = line.split('\t')
            transcript_id,read_counts = fields[0],fields[4]
            expr_dict[transcript_id] = float(read_counts)
    return expr_dict
    
def SR_external_quantification(short_read_gene_matrix_dict,gene_isoforms_length_dict,SR_quantification_option,SR_fastq_list,SR_read_len,ref_annotation,ref_genome,output_dir,threads):
    external_bin_path = os.path.dirname(os.path.realpath(__file__))+'/external_bin'
    temp_dir = f'{output_dir}/temp/'
    Path(temp_dir).mkdir(exist_ok=True,parents=True)
    if SR_quantification_option == 'Kallisto':
        expr_dict = quantify_by_kallisto(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,SR_read_len,threads)
    elif SR_quantification_option == 'Salmon':
        expr_dict = quantify_by_salmon(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,threads)
    elif SR_quantification_option == 'RSEM':
        expr_dict = quantify_by_rsem(ref_annotation,ref_genome,temp_dir,external_bin_path,SR_fastq_list,threads)
    gene_isoform_expression_dict = get_SR_expression_arr(expr_dict,short_read_gene_matrix_dict,gene_isoforms_length_dict)
    return gene_isoform_expression_dict