# miniQuant
**M͟i͟**xed Bayesian **n̲**etwork for **i̲**soform quantification (**miniQuant**) provides a highly-accurate bioinformatics tool for transcript abundance estimation.

**miniQuant** features: 
1. Novel **K-value** metric: a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study
2. **Mixed Bayesian network**: a novel mixed Bayesian network model for transcript abundance estimation that can be applied to different data scenarios: long-read-alone and hybrid (i.e., long reads plus short reads) integrating the strengths of both long reads and short reads.
## Requirements
### Dependency
```
Python>=3.6
Linux operating system
```
The software has been tested with following software version
```
Python==3.9.7
minimap2==2.24
bowtie2==2.4.1
```
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

## Data Preparation
<b>Required:</b>
* long reads alignment data mapped to reference genome in SAM format, example data can be found in `miniQuant/example/LR.sam`
* gene isoform annotation in GTF format, example data can be found in `miniQuant/example/annotation.gtf`
<br>
<b>Optional:</b>
* short reads alignment data mapped to reference transcriptome in SAM format, example data can be found in `miniQuant/example/SR.sam`
<br>
<b>Sequence alignment recommendation:</b>

#### use `minimap2` to map long reads data (e.g. `ENCFF714YOZ.fastq.gz`) to reference genome (e.g. `GRCh38.primary_assembly.genome.fa`)
##### For dRNA-ONT data
```
minimap2 -a --MD -t 10 -N 0 -u f -x splice -o LR.sam 
GRCh38.primary_assembly.genome.fa ENCFF714YOZ.fastq.gz
```
##### For cDNA-ONT or cDNA-PacBio data
```
minimap2 -a --MD -t 10 -N 0 -x splice -o LR.sam 
GRCh38.primary_assembly.genome.fa ENCFF714YOZ.fastq.gz
```
#### use `Bowtie2` to map short reads data (e.g. paired end reads: `ENCFF892WVN.fastq.gz` and `ENCFF481BLH.fastq.gz`) to reference transcriptome (e.g. `gencode.v39.transcripts.fa`)
```
bowtie2-build -f 
gencode.v39.transcripts.fa bowtie2_index

bowtie2 -f --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 10 -k 200 \
-x bowtie2_index -1 ENCFF892WVN.fastq.gz -2 ENCFF481BLH.fastq.gz > SR.sam

```
## Isoform quantification by miniQuant
*miniQuant* provides two options for isoform quantification: 
1. quantify by long reads data alone.

2. quantify using short and long reads data in hybrid mode.

### 1. If quantify using long reads data alone

#### Example: quantify using long reads data (`miniQuant/example/LR.sam`) with annotation (e.g. `miniQuant/example/annotation.gtf`), results in `miniQuant_res` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf miniQuant/example/annotation.gtf \
-lrsam miniQuant/example/LR.sam \
-t 1 \
-o miniQuant_res

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
#### Results explanation 
Isoform quantification abundance <br>
`miniQuant_res/Isoform_abundance.out`
```
Isoform	Gene	TPM
ENST00000373020.9	ENSG00000000003.15	710234.9711212328
ENST00000494424.1	ENSG00000000003.15	0.06848555891537092
ENST00000496771.5	ENSG00000000003.15	103773.58490566035
ENST00000612152.4	ENSG00000000003.15	3.2608726820945185e-20
ENST00000614008.4	ENSG00000000003.15	181274.39435547238
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `TPM`: isoform TPM

### 2. If quantify using short and long reads data in hybrid mode

#### Example: quantify using short reads (e.g. `miniQuant/example/SR.sam`) and long reads data (e.g. `miniQuant/example/SR.sam`) by annotation (e.g. `miniQuant/example/annotation.gtf`), results in `miniQuant_res_hybrid` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf miniQuant/example/annotation.gtf \
-lrsam miniQuant/example/LR.sam \
-srsam miniQuant/example/SR.sam \
--pretrained_model_path dRNA-ONT \
--EM_choice hybrid \
-t 1 \
-o miniQuant_res_hybrid

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
#### Advanced parameters for hybrid quantification
```
optional arguments
  --eff_len_option EFF_LEN_OPTION
                        How to calculate the effective length [kallisto,RSEM]. Choose kallisto 
                        or RSEM to calculate the effective length in the same way as the 
                        corresponding method. Default is kallisto.
  --EM_SR_num_iters EM_SR_NUM_ITERS
                        Number of maximum iterations for EM algorithm. Default is 200.
```
#### Results explanation 
Isoform quantification abundance <br>
`miniQuant_res_hybrid/Isoform_abundance.out`
```
Isoform	Gene	Effective length	TPM
ENST00000373020.9	ENSG00000000003.15	3535.9141630901286	728571.217176296
ENST00000494424.1	ENSG00000000003.15	587.9141630901288	0.08438441577205032
ENST00000496771.5	ENSG00000000003.15	792.9141630901288	97566.37817660523
ENST00000612152.4	ENSG00000000003.15	3563.9141630901286	0.008307999564622307
ENST00000614008.4	ENSG00000000003.15	667.9141630901288	173862.31195468327
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `Effective length`: isoform effective length
* `TPM`: isoform TPM
## Calculate K-value by miniQuant
**K-value** is a key feature of the sequence share pattern that causes particularly high abundance estimation error, allowing us to identify a problematic set of gene isoforms with erroneous quantification that researchers should take extra attention in the study. K-value can be calculated given a gene isoforms annotation in GTF format
#### Example: calculate K-value given annotation in GTF format (e.g. `miniQuant/example/annotation.gtf`)
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py cal_K_value \
-gtf miniQuant/example/annotation.gtf \
-t 1 \
-o miniQuant_kvalue

optional arguments:
  -t THREADS, --threads THREADS
                        Number of threads
  --sr_region_selection SR_REGION_SELECTION
                        SR region selection methods
                        [default:read_length][read_length,num_exons]
```
#### Results explanation 
Kvalue<br>
`miniQuant_kvalue/kvalues.out`
```
Gene	Chr	Num_isoforms	Kvalue
ENSG00000000003.15	chrX	5	14.263027941780145
```
* `Gene`: gene ID
* `Chr`: chromsosome ID
* `Num_isoforms`: number of isoforms in the gene
* `Kvalue`: K value <br>

*For gene that consists only short isoforms (i.e. all isoforms with length <150 bp), K-value will not be calculated and a `NA` value will be given.
