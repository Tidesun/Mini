# miniQuant
**M͟i͟**xed Bayesian **n̲**etwork for **i̲**soform quantification (**miniQuant**) provides a highly-accurate bioinformatics tool for transcript abundance estimation.

**miniQuant** features: 
1. Novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study
2. **Mixed Bayesian network**: a novel mixed Bayesian network model for transcript abundance estimation that can be applied to three different data scenarios: long-read-alone, short-read-alone and hybrid (i.e., long reads plus short reads) integrating the strengths of both long reads and short reads.
## Installation
```
git clone https://github.com/Augroup/miniQuant.git
cd miniQuant
wget -qO- https://miniquant.s3.us-east-2.amazonaws.com/pretrained_models.tar.gz | tar xvz
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
## Isoform quantification by miniQuant
### Quantify using long reads data
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
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
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf GTF_ANNOTATION_PATH \
-lrsam LONG_READS_SAM_PATH \
-srsam SHORT_READS_SAM_PATH \
--pretrained_model_path PRETRAINED_MODEL_PATH \
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
  --pretrained_model_path PRETRAINED_MODEL_PATH
                        The pretrained model path to identify the alpha. default: cDNA-ONT. \n
                        Can be one of the options [cDNA-ONT,dRNA-ONT,cDNA-PacBio] or file path of pretrained model.
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
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py cal_K_value 
-gtf GTF_ANNOTATION_PATH -o OUTPUT_PATH
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
