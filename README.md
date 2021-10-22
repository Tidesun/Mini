# Mili
**M͟i͟**xed **l̲**inear modelling for **i̲**soform quantification (**Mili**)

**Mili** achieves highly-accurate gene isoform quantification using (1) hybrid sequencing (long-read + short-read) data or (2) long-read-alone data.

**Mili** features: 
1. Novel **K-value** metric: Be able to identify gene isoforms with erroneous quantification by other short-reads-alone methods.
2. **Mixed linear model**: Combine the strengths of long reads and short reads to achieve high accuracy in the estimation of isoform abundances, especially for lowly expressed or erronenous isoforms identified through K-value.
## Installation
```
git clone https://github.com/Augroup/Mili.git
cd Mili
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
## Isoform quantification by Mili
```
usage: main.py quantify [-h] -gtf GTF_ANNOTATION_PATH -lrsam
                        LONG_READ_SAM_PATH -o OUTPUT_PATH
                        [-srsam SHORT_READ_SAM_PATH]
                        [-srfastq SHORT_READ_FASTQ]
                        [-sr_m1 SHORT_READ_MATE1_FASTQ]
                        [-sr_m2 SHORT_READ_MATE2_FASTQ]
                        [-ref_genome REFERENCE_GENOME]
                        [--SR_quantification_option SR_QUANTIFICATION_OPTION]
                        [--alpha ALPHA] [--beta BETA] [--filtering FILTERING]
                        [--multi_mapping_filtering MULTI_MAPPING_FILTERING]
                        [--training TRAINING] [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit

required named arguments for TransELS:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

optional arguments:
  -srsam SHORT_READ_SAM_PATH, --short_read_sam_path SHORT_READ_SAM_PATH
                        The path of short read sam file
  -srfastq SHORT_READ_FASTQ, --short_read_fastq SHORT_READ_FASTQ
                        The path of short read fastq file [Required for quantification by other tools.Or use paired end reads]
  -sr_m1 SHORT_READ_MATE1_FASTQ, --short_read_mate1_fastq SHORT_READ_MATE1_FASTQ
                        The path of short read mate 1 fastq file [Required for quantification by other tools.Or use single end reads]
  -sr_m2 SHORT_READ_MATE2_FASTQ, --short_read_mate2_fastq SHORT_READ_MATE2_FASTQ
                        The path of short read mate 2 fastq file [Required for quantification by other tools.Or use single end reads]
  -ref_genome REFERENCE_GENOME, --reference_genome REFERENCE_GENOME
                        The path of reference genome file [Required for quantification by other tools]
  --SR_quantification_option SR_QUANTIFICATION_OPTION
                        SR quantification option[Options: Mili, Kallisto] [default:Mili]
  --alpha ALPHA         Alpha[default:adaptive]: SR and LR balance parameter
  --beta BETA           Beta[default:1e-6]: L2 regularization parameter
  --filtering FILTERING
                        Whether the very short long reads will be
                        filtered[default:False][True,False]
  --multi_mapping_filtering MULTI_MAPPING_FILTERING
                        How to filter multi-mapping
                        reads[default:best][unique_only,best]
  --training TRAINING   Generate training dict
  -t THREADS, --threads THREADS
                        Number of threads
```
The alpha and beta parameters can be set as float. <br>
The `alpha` should be set between 0 and 1, where 0 indicates using only the short reads and 1 indicates using only the long reads for quantification. <br>
The `beta` can be set as a small float between 1e-9 to 1e-2. <br>
By leaving `alpha` and `beta` default(adaptive), the alpha and beta will be obtained by deep learning model to get the optimal performance. <br>
The `multi_mapping_filtering` can be set to `unique_only`, which will use only uniquely-mapping reads; or `best`, which will use the best MAPQ alignment from the multi-mapping reads; or do not set, which will not perform any filtering on multi_mapping reads. Default is `best`<br>
Set `SR_quantification_option` to use other tools for short reads quantification. For other tools, you need also provide `reference_genome` and `srfastq` for single end reads / `sr_m1` and `sr_m2` for paired end reads.

## Calculate K-value
```
usage: main.py cal_K_value [-h] -gtf GTF_ANNOTATION_PATH -lrsam LONG_READ_SAM_PATH
                      -o OUTPUT_PATH [-t THREADS]
                      [--sr_region_selection SR_REGION_SELECTION]
                      [--filtering FILTERING]

optional arguments:
  -h, --help            show this help message and exit

required named arguments for TrEESR:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

optional arguments:
  -t THREADS, --threads THREADS
                        Number of threads
  --sr_region_selection SR_REGION_SELECTION
                        SR region selection methods
                        [default:read_length][read_length,num_exons]
  --filtering FILTERING
                        Whether the very short long reads will be
                        filtered[default:True][True,False]
```
