import numpy as np
import pysam
def get_eff_len_kallisto(isoform_len, fraglen, fragsd):
    if isoform_len < fraglen - 5*fragsd:
        return 1
    sqrt_2pi = np.sqrt(2*np.pi)
    exp_term = np.exp(-1/2*((np.arange(1, isoform_len+1)-fraglen)/fragsd)**2)
    numerator = np.arange(1, isoform_len+1) * exp_term / (fragsd*sqrt_2pi)
    denominator = exp_term / (fragsd*sqrt_2pi)
    eff_len = isoform_len - np.sum(numerator) / np.sum(denominator) + 1
    return max(eff_len, 1)
def get_eff_len_rsem(isoform_len, fraglen, fragsd):
    if isoform_len < fraglen - 5*fragsd:
        return 1
    sqrt_2pi = np.sqrt(2*np.pi)
    exp_term = np.exp(-1/2*((np.arange(1, isoform_len+1)-fraglen)/fragsd)**2)
    numerator = (isoform_len - np.arange(1, isoform_len+1) + 1) * exp_term / (fragsd*sqrt_2pi)
    denominator = exp_term / (fragsd*sqrt_2pi)
    eff_len = np.sum(numerator) / np.sum(denominator)
    return max(eff_len, 1)
def get_eff_len_dict(SR_sam,mean_f_len,std_f_len,option='Kallisto'):
    eff_len_dict = {}
    with pysam.AlignmentFile(SR_sam, "r") as f:
        for isoform,isoform_len in zip(f.references,f.lengths):
            if '|' in isoform:
                isoform = isoform.split('|')[0]
            if option == 'Kallisto':
                eff_len_dict[isoform] = get_eff_len_kallisto(isoform_len,mean_f_len,std_f_len)
            elif option == 'RSEM':
                eff_len_dict[isoform] = get_eff_len_rsem(isoform_len,mean_f_len,std_f_len)
    return eff_len_dict