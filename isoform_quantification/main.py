import argparse
from TrEESR import TrEESR
from TransELS import TransELS

def parse_arguments():
    """
    Parse the arguments
    """
    parser = argparse.ArgumentParser(description="Isoform quantification tools",add_help=True)
    subparsers = parser.add_subparsers(help='sub-command help',dest="subparser_name")
    parser_TrEESR = subparsers.add_parser('TrEESR', help='TrEESR')
    parser_TransELS = subparsers.add_parser('TransELS', help='TransELS')
    
    requiredNamed_TrEESR = parser_TrEESR.add_argument_group('required named arguments for TrEESR')
    requiredNamed_TrEESR.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TrEESR.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    optional_TrEESR = parser_TrEESR.add_argument_group('optional arguments')
    optional_TrEESR.add_argument('-t','--threads',type=int, default=1,help="Number of threads")

    requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for TrEESR')
    requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=True)
    requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    optional_TransELS.add_argument('--b_cal_method',type=str,default='original', help="Region expression calculation method ['original','coverage','div_read_length']")
    optional_TransELS.add_argument('--alpha',type=str,default='adaptive', help="Alpha")
    optional_TransELS.add_argument('--beta',type=str, default='adaptive',help="Beta")
    optional_TransELS.add_argument('--P',type=float, default=1e-6,help="P")
    optional_TransELS.add_argument('-t','--threads',type=int, default=1,help="Number of threads")
    
    args = parser.parse_args()
    if args.subparser_name == 'TrEESR':
        print('Using TrEESR')
        TrEESR(args.gtf_annotation_path,args.output_path,args.threads)
    elif args.subparser_name == 'TransELS':
        print('Using TransELS')
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
        TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,args.b_cal_method,alpha,beta,args.P,args.threads)
    else:
        parser.print_help()
if __name__ == "__main__":
    parse_arguments()