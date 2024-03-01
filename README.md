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
## Data Preparation
<b>REQUIRED:</b>
* long reads data in FASTQ format, example data can be downloaded from [ENCODE](https://www.encodeproject.org/files/ENCFF714YOZ/@@download/ENCFF714YOZ.fastq.gz)
* reference genome in FASTA format, can be downloaded from [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.annotation.gtf.gz) for human
* gene isoform annotation in GTF format, can be downloaded from [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz) for human
<br>

<b>OPTIONAL:</b>
* short reads data in FASTQ format, example data can be downloaded from ENCODE [pair1](https://www.encodeproject.org/files/ENCFF892WVN/@@download/ENCFF892WVN.fastq.gz) and [pair2](https://www.encodeproject.org/files/ENCFF481BLH/@@download/ENCFF481BLH.fastq.gz)
* reference transcriptome in FASTA format, can be downloaded from [GENCODE](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.transcripts.fa.gz) for human
## Isoform quantification by miniQuant

### 1. If quantify using long reads data alone

#### Step 0: use `minimap2` to map long reads data (e.g. `ENCFF714YOZ.fastq.gz`) to reference genome (e.g. `GRCh38.primary_assembly.genome.fa`)
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
#### Step 1: quantify using long reads data (`LR.sam`) by GENCODE annotation (e.g. `gencode.v39.annotation.gtf`), results in `miniQuant_res` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf gencode.v39.annotation.gtf \
-lrsam LR.sam \
-t 10 \
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
ENST00000002501.11	ENSG00000003249.15	28.213415252281333
ENST00000002829.8	ENSG00000001617.12	5.806563371373737
ENST00000003912.7	ENSG00000001461.17	3.066295789453613e-07
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `TPM`: isoform TPM

### 2. If quantify using short and long reads data in hybrid mode

#### Step 0: use `Bowtie2` to map short reads data (e.g. paired end reads: `ENCFF892WVN.fastq.gz` and `ENCFF481BLH.fastq.gz`) to reference transcriptome (e.g. `gencode.v39.transcripts.fa`)
```
bowtie2-build -f 
gencode.v39.transcripts.fa bowtie2_index

bowtie2 -f --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 10 -k 200 \
-x bowtie2_index -1 ENCFF892WVN.fastq.gz -2 ENCFF481BLH.fastq.gz > SR.sam

```
#### Step 1: quantify using short reads (e.g. `SR.sam`) and long reads data (e.g. `LR.sam`) by GENCODE annotation (e.g. `gencode.v39.annotation.gtf`), results in `miniQuant_res` folder
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py quantify \
-gtf gencode.v39.annotation.gtf \
-lrsam LR.sam \
-srsam SR.sam \
--pretrained_model_path PRETRAINED_MODEL_PATH \
--EM_choice hybrid \
-t 10 \
-o miniQuant_res

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
                        How to calculate the effective length [Kallisto,RSEM]. Choose Kallisto 
                        or RSEM to calculate the effective length in the same way as the 
                        corresponding method. Default is Kallisto.
  --EM_SR_num_iters EM_SR_NUM_ITERS
                        Number of maximum iterations for EM algorithm. Default is 200.
```
#### Results explanation 
Isoform quantification abundance <br>
`miniQuant_res/Isoform_abundance.out`
```
Isoform	Gene	Effective length	TPM
ENST00000003084.11	ENSG00000001626.17  5834.082567393011	0.02413423715032308
ENST00000005257.7	ENSG00000006451.8	2551.0825673930112  54.70210536985984	
ENST00000005386.8	ENSG00000005175.10	4103.082567393011	7.728764629572447
```
* `Isoform`: isoform ID
* `Gene`: gene ID
* `Effective length`: isoform effective length
* `TPM`: isoform TPM
## Calculate K-value by GENCODE annotation
```
source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py cal_K_value \
-gtf gencode.v39.annotation.gtf \
-t 10 \
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

```
* `Gene`: gene ID
* `Chr`: chr ID
* `Num_isoforms`: number of isoforms in the gene
* `Kvalue`: K value <br>

*For gene that only consists of a very short exon, and that exon was filtered out by the strategy defined by `--sr_region_selection`, the k value will not be calculated and a `NA` value will be given.
