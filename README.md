# Mini
**M͟i͟**xed Bayesian **n**etwork for **i̲**soform quantification (**Mini**)

**Mini** achieves highly-accurate gene isoform quantification using (1) hybrid sequencing (long-read + short-read) data or (2) long-read-alone data.

**Mini** features: 
1. Novel **K-value** metric: Be able to identify gene isoforms with erroneous quantification by other short-reads-alone methods.
2. **Mixed Bayesian network**: Combine the strengths of long reads and short reads to achieve high accuracy in the estimation of isoform abundances, especially for lowly expressed or erronenous isoforms identified through K-value.
## Installation
```
git clone https://github.com/Augroup/Mini.git
cd Mini
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
## Isoform quantification by Mini
### Quantify using long reads data
```
python isoform_quantification/main.py quantify \
-gtf GTF_ANNOTATION_PATH \
-lrsam LONG_READS_SAM_PATH \
-t 1 \
-o OUTPUT_PATH

arguments:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of isoform annotation file in GTF format
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file mapping to reference genome.
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory
```

### Quantify using short and long reads data in hybrid mode
```
python isoform_quantification/main.py quantify \
-gtf GTF_ANNOTATION_PATH \
-lrsam LONG_READS_SAM_PATH \
-srsam SHORT_READS_SAM_PATH \
--EM_choice hybrid \
-t 1 \
-o OUTPUT_PATH

arguments:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of isoform annotation file in GTF format
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file mapping to reference genome.
  -srsam SHORT_READ_SAM_PATH, --short_read_sam_path SHORT_READ_SAM_PATH
                        The path of short read sam file mapping to reference transcriptome.
  -t THREADS, --threads THREADS
                        Number of threads. Default is 1.
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory
```
### Advanced parameters for hybrid quantification
```
optional arguments
  --iter_theta ITER_THETA
                        Whether use updated theta to re-calculate conditional probablity [True,False].
                        See the online methods for detailed information. Default is False to 
                        speed up the calculation.
  --eff_len_option EFF_LEN_OPTION
                        How to calculate the effective length [Kallisto,RSEM]. Choose Kallisto 
                        or RSEM to calculate the effective length in the same way as the 
                        corresponding method. Default is Kallisto.
  --EM_SR_num_iters EM_SR_NUM_ITERS
                        Number of maximum iterations for EM algorithm. Default is 200.
  --inital_theta INITAL_THETA
                        Inital_theta [LR,SR]. Set the initial theta based on the isoform 
                        expression given long reads(LR) or short reads (SR). Default is LR.
```

## Calculate K-value
```
usage: main.py cal_K_value [-h] -gtf GTF_ANNOTATION_PATH -o OUTPUT_PATH
                           [-lrsam LONG_READ_SAM_PATH] [-t THREADS]
                           [--sr_region_selection SR_REGION_SELECTION]
                           [--filtering FILTERING]

optional arguments:
  -h, --help            show this help message and exit

required named arguments for calculation of K value:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

optional arguments:
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file
  -t THREADS, --threads THREADS
                        Number of threads
  --sr_region_selection SR_REGION_SELECTION
                        SR region selection methods
                        [default:read_length][read_length,num_exons]
  --filtering FILTERING
                        Whether the very short long reads will be
                        filtered[default:True][True,False]
```
For gene that only consists of a very short exon, and that exon was filtered out by the strategy defined by `--sr_region_selection`, the k value will not be calculated and a `NA` value will be given.
