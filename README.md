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
usage: main.py TrEESR [-h] -gtf GTF_ANNOTATION_PATH -o OUTPUT_PATH
                      [-t THREADS]

optional arguments:
  -h, --help            show this help message and exit

required named arguments for TrEESR:
  -gtf GTF_ANNOTATION_PATH, --gtf_annotation_path GTF_ANNOTATION_PATH
                        The path of annotation file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

optional arguments:
  -t THREADS, --threads THREADS
                        Number of threads
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
  -srsam SHORT_READ_SAM_PATH, --short_read_sam_path SHORT_READ_SAM_PATH
                        The path of short read sam file
  -lrsam LONG_READ_SAM_PATH, --long_read_sam_path LONG_READ_SAM_PATH
                        The path of long read sam file
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        The path of output directory

optional arguments:
  --b_cal_method B_CAL_METHOD
                        Region expression calculation method
                        ['original','coverage','div_read_length']
  --alpha ALPHA         Alpha
  --beta BETA           Beta
  --P P                 P
  -t THREADS, --threads THREADS
                        Number of threads
```
