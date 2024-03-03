import numpy as np
import pysam
def get_eff_len_kallisto(isoform_len,fraglen,fragsd):
    if (isoform_len < fraglen - 5*fragsd):
        return (1)
    sum_numerator = 0
    sum_denominator = 0
    for j in range(1,isoform_len+1):
        sum_numerator = sum_numerator + j * 1/(fragsd*np.sqrt(2*np.pi)) * np.e**(-1/2*((j-fraglen)/fragsd)**2)
        sum_denominator = sum_denominator + 1/(fragsd*np.sqrt(2*np.pi)) * np.e**(-1/2*((j-fraglen)/fragsd)**2)
    eff_len = isoform_len - sum_numerator/sum_denominator +1
    if (eff_len < 1):
        eff_len = 1
    return eff_len
def get_eff_len_rsem(isoform_len,fraglen,fragsd):
    if (isoform_len < fraglen - 5*fragsd):
        return (1)
    sum_numerator = 0
    sum_denominator = 0
    for j in range(1,isoform_len+1):
        sum_numerator = sum_numerator + (isoform_len-j+1) * 1/(fragsd*np.sqrt(2*np.pi)) * np.e**(-1/2*((j-fraglen)/fragsd)**2)
        sum_denominator = sum_denominator + 1/(fragsd*np.sqrt(2*np.pi)) * np.e**(-1/2*((j-fraglen)/fragsd)**2)
    eff_len = sum_numerator/sum_denominator
    if (eff_len < 1):
        eff_len = 1
    return eff_len
def get_eff_len_dict(SR_sam,mean_f_len,std_f_len,option='kallisto'):
    eff_len_dict = {}
    with pysam.AlignmentFile(SR_sam, "r") as f:
        for isoform,isoform_len in zip(f.references,f.lengths):
            if '|' in isoform:
                isoform = isoform.split('|')[0]
            if option == 'kallisto':
                eff_len_dict[isoform] = get_eff_len_kallisto(isoform_len,mean_f_len,std_f_len)
            elif option == 'RSEM':
                eff_len_dict[isoform] = get_eff_len_rsem(isoform_len,mean_f_len,std_f_len)
    return eff_len_dict