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
    
    requiredNamed_TransELS = parser_TransELS.add_argument_group('required named arguments for TrEESR')
    requiredNamed_TransELS.add_argument('-gtf','--gtf_annotation_path', type=str, help="The path of annotation file",required=True)
    requiredNamed_TransELS.add_argument('-srsam','--short_read_sam_path', type=str, help="The path of short read sam file",required=True)
    requiredNamed_TransELS.add_argument('-lrsam','--long_read_sam_path', type=str, help="The path of long read sam file",required=True)
    requiredNamed_TransELS.add_argument('-o','--output_path', type=str, help="The path of output directory",required=True)
    
    optional_TransELS = parser_TransELS.add_argument_group('optional arguments')
    optional_TransELS.add_argument('--b_cal_method',type=str,default='coverage', help="Region expression calculation method ['coverage','div_read_length']")
    optional_TransELS.add_argument('--alpha',type=str,default='adaptive', help="Alpha")
    optional_TransELS.add_argument('--beta',type=str, default='adaptive',help="Beta")
    optional_TransELS.add_argument('--P',type=float, default=1e-6,help="P")
#     optional.add_argument(
#     '-h',
#     '--help',
#     action='help',
#     default=argparse.SUPPRESS,
#     help='show this help message and exit')
    
    args = parser.parse_args()
    if args.subparser_name == 'TrEESR':
        print('Using TrEESR')
        print(args)
        TrEESR(args.gtf_annotation_path,args.output_path)
    else:
        print('Using TransELS')
        print(args)
        TransELS(args.gtf_annotation_path,args.short_read_sam_path,args.long_read_sam_path,args.output_path,args.b_cal_method,args.alpha,args.beta,args.P)
if __name__ == "__main__":
    parse_arguments()
