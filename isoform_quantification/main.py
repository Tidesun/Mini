import argparse
from TrEESR import TrEESR
from TransELS import TransELS
from EM import EM,EM_SR,EM_hybrid
import config
import os
# import os
# os.system("taskset -p 0xfffff %d" % os.getpid())
# affinity_mask = os.sched_getaffinity(0)
# os.sched_setaffinity(0, affinity_mask)
def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_TrEESR = subparsers.add_parser('cal_K_value', aliases=['TrEESR'],help='Calculate K values')
    # parser_TransELS = subparsers.add_parser('quantify', aliases=['TransELS'],help='Isoform quantification')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for calculation of K value')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=False)
    optional_TrEESR.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False)
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_TrEESR.add_argument('--sr_region_selection',type=str, default='real_data',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    optional_TrEESR.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:True][True,False]")
    optional_TrEESR.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_TrEESR.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_TrEESR.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_TrEESR.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_TrEESR.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    optional_TrEESR.add_argument('--keep_sr_exon_region',type=str, default='True',help="Keep exon region for SR if using real data to filter region [default:True][True,False]")
    optional_TrEESR.add_argument('--use_weight_matrix',type=str, default='True',help="Whether use weight matrix[default:True][True,False]")
    optional_TrEESR.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    optional_TrEESR.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    # weight_path = os.path.dirname(os.path.realpath(__file__))+'/weights/nanosim_weight_dict.pkl'
    # assert os.path.exists(weight_path)
    # optional_TrEESR.add_argument('--region_weight_path',type=str, default=weight_path,help="Mili LR region weight path")

    # requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for isoform quantification')
    # requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    # requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    # requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    # optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    # optional_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    # optional_TransELS.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    # optional_TransELS.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    # optional_TransELS.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    # optional_TransELS.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    # optional_TransELS.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, Kallisto,Salmon, RSEM] [default:Kallisto]",default='Kallisto')
    # # optional_TransELS.add_argument('--kallisto_index', type=str, help="Kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/KallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    # optional_TransELS.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    # optional_TransELS.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    # optional_TransELS.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    # optional_TransELS.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    # optional_TransELS.add_argument('--training',type=str,default='False', help="Generate training dict")
    # optional_TransELS.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    # optional_TransELS.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    # optional_TransELS.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    # optional_TransELS.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    # optional_TransELS.add_argument('--use_weight_matrix',type=str, default='True',help="Whether use weight matrix[default:True][True,False]")
    # optional_TransELS.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    # optional_TransELS.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    # optional_TransELS.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    # optional_TransELS.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    # optional_TransELS.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    # optional_TransELS.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    # optional_TransELS.add_argument('--sr_region_selection',type=str, default='real_data',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    # optional_TransELS.add_argument('--keep_sr_exon_region',type=str, default='True',help="Keep exon region for SR if using real data to filter region [default:True][True,False]")
    # optional_TransELS.add_argument('--region_weight_path',type=str, default=weight_path,help="Mili LR region weight path")
   
    parser_EM = subparsers.add_parser('quantify', aliases=['EM'],help='Isoform quantification by EM algorithm')
    requiredNamed_EM = parser_EM.add_argument_group('required named arguments for isoform quantification')
    requiredNamed_EM.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_EM.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_EM = parser_EM.add_argument_group('optional arguments')
    optional_EM.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=False,default=None)
    optional_EM.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    optional_EM.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    optional_EM.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    optional_EM.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    optional_EM.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    optional_EM.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, Kallisto,Salmon, RSEM] [default:Kallisto]",default='Kallisto')
    # optional_EM.add_argument('--kallisto_index', type=str, help="Kallisto index",default='/fs/project/PCON0009/Yunhao/Project/Mili/Annotation/KallistoIndex/gencode.v39.transcripts.clean.dedup.m')
    optional_EM.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    optional_EM.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    optional_EM.add_argument('--filtering',type=str,default='False', help="Whether the very short long reads will be filtered[default:False][True,False]")
    optional_EM.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    optional_EM.add_argument('--training',type=str,default='False', help="Generate training dict")
    optional_EM.add_argument('--DL_model',type=str,default=None,help='DL model to use')
    optional_EM.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    optional_EM.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_EM.add_argument('--READ_JUNC_MIN_MAP_LEN',type=int, default=1,help="minimum mapped read length to consider a junction")
    optional_EM.add_argument('--use_weight_matrix',type=str, default='True',help="Whether use weight matrix[default:True][True,False]")
    optional_EM.add_argument('--normalize_lr_A',type=str, default='True',help="Whether normalize lr A [default:True] [True,False]")
    optional_EM.add_argument('--same_struc_isoform_handling',type=str, default='merge',help="How to handle isoforms with same structures within a gene[default:merge][merge,keep]")
    optional_EM.add_argument('--add_full_length_region',type=str, default='all',help="Whether add full length region[default:all] [all,nonfullrank,none]")
    optional_EM.add_argument('--multi_exon_region_weight',type=str, default='regular',help="The weight in matrix A for multi_exon_region[default:regular][regular,minus_inner_region]")
    optional_EM.add_argument('--output_matrix_info',type=str, default='False',help="Whether output matrix info [default:False] [True,False]")
    optional_EM.add_argument('--normalize_sr_A',type=str, default='True',help="Whether normalize sr A [default:False] [True,False]")
    optional_EM.add_argument('--sr_region_selection',type=str, default='real_data',help="SR region selection methods [default:real_data][read_length,num_exons,real_data]")
    optional_EM.add_argument('--keep_sr_exon_region',type=str, default='True',help="Keep exon region for SR if using real data to filter region [default:True][True,False]")
    # optional_EM.add_argument('--region_weight_path',type=str, default=weight_path,help="Mili LR region weight path")
    optional_EM.add_argument('--EM_choice',type=str, default='LIQA_modified',help="EM_choice[LIQA,LIQA_modified,LR]")
    optional_EM.add_argument('--iter_theta',type=str, default='False',help="Whether use updated theta to re-calculate conditional prob [True,False]")
    optional_EM.add_argument('--kde_path',type=str, default='/fs/project/PCON0009/Au-scratch2/haoran/_projects/long_reads_rna_seq_simulator/models/kde_H1-hESC_dRNA',help="KDE model path")
    optional_EM.add_argument('--eff_len_option',type=str, default='Kallisto',help="Calculation of effective length option [Kallisto,RSEM]")
    optional_EM.add_argument('--EM_SR_num_iters',type=int, default=200,help="Number of EM SR iterations")
    optional_EM.add_argument('--alpha_df_path',type=str, default=None,help="Alpha df path")
    optional_EM.add_argument('--inital_theta',type=str, default='LR',help="inital_theta [LR,SR]")

    args = parser.parse_args()
    if args.filtering == 'True':
        args.filtering = True
    else:
        args.filtering = False
    config.same_struc_isoform_handling = args.same_struc_isoform_handling
    config.READ_JUNC_MIN_MAP_LEN = args.READ_JUNC_MIN_MAP_LEN
    config.multi_exon_region_weight = args.multi_exon_region_weight
    config.sr_region_selection = args.sr_region_selection
    # config.region_weight_path = args.region_weight_path
    if args.output_matrix_info == 'True':
        config.output_matrix_info = True
    else:
        config.output_matrix_info = False
    if args.keep_sr_exon_region == 'True':
        config.keep_sr_exon_region = True
    else:
        config.keep_sr_exon_region = False
    if args.normalize_sr_A == 'True':
        config.normalize_sr_A = True
    else:
        config.normalize_sr_A = False
    if args.normalize_lr_A == 'True':
        config.normalize_lr_A = True
    else:
        config.normalize_lr_A = False
    if args.use_weight_matrix == 'True':
        config.use_weight_matrix = True
    else:
        config.use_weight_matrix = False
    config.add_full_length_region = args.add_full_length_region
    # config.kallisto_index = args.kallisto_index
    print('\n'.join(f'{k}={v}' for k, v in vars(args).items()))
    if args.subparser_name in ['cal_K_value','TrEESR']:
        print('Calculate K values')
        TrEESR(args.gtf_annotation_path,args.output_path,args.short_read_sam_path,args.long_read_sam_path,args.sr_region_selection,args.filtering,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    # elif args.subparser_name in ['quantify','TransELS']:
    #     if args.training == 'True':
    #         args.training = True
    #     else:
    #         args.training = False
    #     print('Isoform quantification',flush=True)
    #     if (args.short_read_sam_path is None) or (args.alpha == 1.0):
    #         args.alpha = 1.0
    #         args.SR_quantification_option = 'Mili'
    #     if args.alpha != 1.0:
    #         if args.short_read_sam_path is None:
    #             raise Exception('You need to provide a short read alignemnt file if the alpha is not 1!')
    #         if args.SR_quantification_option != 'Mili':
    #             if (args.short_read_fastq is None) and (args.short_read_mate1_fastq is None or args.short_read_mate2_fastq is None):
    #                 raise Exception('You need to provide the single end or paried end SR fastq if using other short read quantification options!')
    #             if (args.reference_genome is None):
    #                 raise Exception('You need to provide the reference genome if using other short read quantification options!')
    #     if (args.alpha == 'adaptive'):
    #         alpha = 'adaptive'
    #     else:
    #         try:
    #             alpha = float(args.alpha)
    #         except:
    #             raise Exception('Alpha given is not numeric')
    #     if (args.beta == 'adaptive'):
    #         beta = 'adaptive'
    #     else:
    #         try:
    #             beta = float(args.beta)
    #         except:
    #             raise Exception('Beta given is not numeric')
    #     if args.SR_quantification_option not in ['Mili','Kallisto','Salmon','RSEM']:
    #         raise Exception('SR_quantification_option is not valid.Options: [Mili, Kallisto,Salmon, RSEM]')
    #     if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
    #         args.multi_mapping_filtering = 'no_filtering'
    #     SR_fastq_list = []
    #     if args.short_read_fastq is not None:
    #         SR_fastq_list = [args.short_read_fastq]
    #     elif args.short_read_mate1_fastq is not None:
    #         SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
    #     if args.DL_model is None:
    #         args.DL_model = args.SR_quantification_option + '.pt'
    #     TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN)
    elif args.subparser_name in ['quantify','EM']:
        config.kde_path = args.kde_path
        if args.training == 'True':
            args.training = True
        else:
            args.training = False
        print('Isoform quantification by Mini',flush=True)
        if (args.short_read_sam_path is None) or (args.alpha == 1.0):
            args.alpha = 1.0
            args.SR_quantification_option = 'Mini'
        # if args.alpha != 1.0:
        #     if args.short_read_sam_path is None:
        #         raise Exception('You need to provide a short read alignemnt file if the alpha is not 1!')
        #     if args.SR_quantification_option != 'Mili':
        #         if (args.short_read_fastq is None) and (args.short_read_mate1_fastq is None or args.short_read_mate2_fastq is None):
        #             raise Exception('You need to provide the single end or paried end SR fastq if using other short read quantification options!')
        #         if (args.reference_genome is None):
        #             raise Exception('You need to provide the reference genome if using other short read quantification options!')
        if (args.alpha == 'adaptive'):
            alpha = 'adaptive'
        else:
            try:
                alpha = float(args.alpha)
            except:
                raise Exception('Alpha given is not numeric')
        if (args.beta == 'adaptive'):
            beta = 'adaptive'
        else:
            try:
                beta = float(args.beta)
            except:
                raise Exception('Beta given is not numeric')
        # if args.SR_quantification_option not in ['Mili','Kallisto','Salmon','RSEM']:
        #     raise Exception('SR_quantification_option is not valid.Options: [Mili, Kallisto,Salmon, RSEM]')
        if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
            args.multi_mapping_filtering = 'no_filtering'
        SR_fastq_list = []
        if args.short_read_fastq is not None:
            SR_fastq_list = [args.short_read_fastq]
        elif args.short_read_mate1_fastq is not None:
            SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
        if args.DL_model is None:
            args.DL_model = args.SR_quantification_option + '.pt'
        config.alpha = args.alpha
        config.alpha_df_path = args.alpha_df_path
        config.inital_theta = args.inital_theta
        if args.EM_choice == 'SR':
            config.eff_len_option = args.eff_len_option
            config.EM_SR_num_iters = args.EM_SR_num_iters
            args.long_read_sam_path = None
            args.alpha = 0
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        elif args.EM_choice == 'hybrid':
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)
        else:
            if args.EM_choice == 'LR':
                args.EM_choice = 'LIQA_modified'
            args.short_read_sam_path = None
            args.alpha = 1
            EM_hybrid(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads,READ_JUNC_MIN_MAP_LEN=args.READ_JUNC_MIN_MAP_LEN,EM_choice=args.EM_choice,iter_theta=args.iter_theta)    
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()
