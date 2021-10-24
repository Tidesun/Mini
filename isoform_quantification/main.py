import argparse
from TrEESR import TrEESR
from TransELS import TransELS
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
    parser_TransELS = subparsers.add_parser('quantify', aliases=['TransELS'],help='Isoform quantification')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for calculation of K value')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    optional_TrEESR.add_argument('--sr_region_selection',type=str, default='read_length',help="SR region selection methods [default:read_length][read_length,num_exons]")
    optional_TrEESR.add_argument('--filtering',type=bool,default=True, help="Whether the very short long reads will be filtered[default:True][True,False]")

    requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for isoform quantification')
    requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    optional_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",default=None)
    optional_TransELS.add_argument('-srfastq','--short_read_fastq', type=str, help="The path of short read fastq file",default=None)
    optional_TransELS.add_argument('-sr_m1','--short_read_mate1_fastq', type=str, help="The path of short read mate 1 fastq file",default=None)
    optional_TransELS.add_argument('-sr_m2','--short_read_mate2_fastq', type=str, help="The path of short read mate 2 fastq file",default=None)

    optional_TransELS.add_argument('-ref_genome','--reference_genome', type=str, help="The path of reference genome file",default=None)
    optional_TransELS.add_argument('--SR_quantification_option', type=str, help="SR quantification option[Options: Mili, Kallisto,Salmon, RSEM] [default:Mili]",default='Mili')
    optional_TransELS.add_argument('--alpha',type=str,default='adaptive', help="Alpha[default:adaptive]: SR and LR balance parameter")
    optional_TransELS.add_argument('--beta',type=str, default='1e-6',help="Beta[default:1e-6]: L2 regularization parameter")
    optional_TransELS.add_argument('--filtering',type=bool,default=False, help="Whether the very short long reads will be filtered[default:False][True,False]")
    optional_TransELS.add_argument('--multi_mapping_filtering',type=str,default='best', help="How to filter multi-mapping reads[default:best][unique_only,best]")
    optional_TransELS.add_argument('--training',type=bool,default=False, help="Generate training dict")
    optional_TransELS.add_argument('--DL_model',type=str,default='rsem_model_20.pt',help='DL model to use')
    optional_TransELS.add_argument('--assign_unique_mapping_option',type=str,default='linear_model',help='How to assign unique mapping reads [Options:linear_model,manual_assign] [default:linear_model]')
    optional_TransELS.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    args = parser.parse_args()
    if args.subparser_name == ['cal_K_value','TrEESR']:
        print('Calculate K values')
        TrEESR(args.gtf_annotation_path,args.output_path,args.long_read_sam_path,args.sr_region_selection,args.filtering,args.threads)
    elif args.subparser_name in ['quantify','TransELS']:
        print('Isoform quantification',flush=True)
        if (args.short_read_sam_path is None):
            args.alpha = 1.0
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
        if (args.multi_mapping_filtering is None) or (not args.multi_mapping_filtering in ['unique_only','best']):
            args.multi_mapping_filtering = 'no_filtering'
        SR_fastq_list = []
        if args.short_read_fastq is not None:
            SR_fastq_list = [args.short_read_fastq]
        elif args.short_read_mate1_fastq is not None:
            SR_fastq_list = [args.short_read_mate1_fastq,args.short_read_mate2_fastq]
        print(args)
        TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,alpha,beta,1e-6,args.filtering,args.multi_mapping_filtering,args.SR_quantification_option,SR_fastq_list,args.reference_genome,args.training,args.DL_model,args.assign_unique_mapping_option,args.threads)
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()
