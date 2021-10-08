# isoform_quantification
## Installation
```
git clone https://github.com/Tidesun/isoform_quantification.git
cd isoform_quantification
python -m venv base
source base/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```
## TrEESR
```
usage: main.py TrEESR [-h] -gtf GTF_ANNOTATION_PATH -lrsam LONG_READ_SAM_PATH
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
## TransELS
```
usage: main.py TransELS [-h] -gtf GTF_ANNOTATION_PATH -srsam
                        SHORT_READ_SAM_PATH -lrsam LONG_READ_SAM_PATH -o
                        OUTPUT_PATH [--b_cal_method B_CAL_METHOD]
                        [--alpha ALPHA] [--beta BETA] [--P P] [-t THREADS]

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
  -srsam SHORT_READ_SAM_PATH, --short_read_sam_path SHORT_READ_SAM_PATH
                        The path of short read sam file
  --alpha ALPHA         Alpha[default:adaptive]: SR and LR balance parameter
  --beta BETA           Beta[default:adaptive]: L2 regularization parameter
  --filtering FILTERING
                        Whether the very short long reads will be filtered[default:True][True,False]
  --multi_mapping_filtering MULTI_MAPPING_FILTERING
                        How to filter multi-mapping
                        reads[default:best][unique_only,best]
  -t THREADS, --threads THREADS
                        Number of threads
```
The alpha and beta parameters can be set as float. <br>
The `alpha` should be set between 0 and 1, where 0 indicates using only the short reads and 1 indicates using only the long reads for quantification. <br>
The `beta` can be set as a small float between 1e-9 to 1e-2. <br>
By leaving `alpha` and `beta` default(adaptive), the alpha and beta will be obtained by deep learning model to get the optimal performance. <br>
The `multi_mapping_filtering` can be set to `unique_only`, which will use only uniquely-mapping reads; or `best`, which will use the best MAPQ alignment from the multi-mapping reads; or do not set, which will not perform any filtering on multi_mapping reads. Default is `best`