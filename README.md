# hellobc
This is a toolkit for upstream analysis of scRNA-seq data developed by Youchun Xu's group, Tsinghua University.

The toolkit contains the following functions:
* Parsing user-defined barcode patterns
* Quantifying RNA
* Cell calling
* Generating count matrix

# Preperation
Before using hibc, make sure you have installed STAR, featureCounts.
```
# install STAR
# install featureCounts
```

# Installation
We recommend installing hibc in a conda environment.
```
conda create -n hellobc python=3.8
conda activate hellobc
conda install bioconda::seqkit
conda install bioconda::samtools
pip install hellobc
```

# Get started
In order to run the whole pipeline of upstream analysis, you can create a new .py file with the following code:
```
from hibc.pipline.runpip import hellobc_pipeline

workst = "/mnt/sda/xxt/workst/sevci/20240517"                     # output directory
sampid = "ZT-136-ZT-136-1"                                        # sample name
fq1 = "/mnt/sda/xxt/data/seqdata/biozt/ZT-136-ZT-136-1/ZT-136-ZT-136-1_R1.fastq.gz"
fq2 = "/mnt/sda/xxt/data/seqdata/biozt/ZT-136-ZT-136-1/ZT-136-ZT-136-1_R2.fastq.gz"
bcpat = '36B1[1:12; 19:30; 37:48] 8U1[49:56]'                     # barcode pattern of biocapital
# bcpat = '16B1[1:16] 12U1[17:28]'                                # barcode pattern of 10x chem_v2

# path to STAR, featureCounts. If not setted, the STAR and featureCounts in the environmental variables will be used.
star_path = '/mnt/sda/xxt/soft/STAR-2.6.1a/bin/Linux_x86_64_static/STAR'
fc_path = '/mnt/sda/xxt/soft/subread-2.0.6-Linux-x86_64/bin/featureCounts'

# human genome
gtf_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
star_index_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-GRCh38-2020-A/star"

# mouse genome
# mouse_gtf_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-mm10-2020-A/genes/genes.gtf"
# mouse_star_index_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-mm10-2020-A/star"

xu_pipeline(workst, sampid, fq1, fq2, bcpat, gtf_dir, star_index_dir, init_method='ordmag', min_umis=500, max_adj_pvalue=0.01, star_path=star_path, fc_path=fc_path, corenum=-1)
```
# More details
## Input format of barcode pattern
The composition of the barcode format string:
* Letters: 'B' means cell barcode, 'U' means umi.
* The number in front of the letter: total length (number of bases) of barcode or umi.
* The number after the letter: can only choose between 1 or 2. 1 means barcode/umi is in read1, 2 means read2.
* Contents in [ ]: position of each barcode/umi segment in the sequence, separated by ';'.
* For example:
```
# If you have 3 discrete barcode segments (each segment is 12nt) and 1 segment of 8nt umi, both barcode and umi are in read1:
bcpat = '36B1[1:12; 19:30; 37:48] 8U1[49:56]'
# If you have 1 segment of 16nt barcode and 1 segment of 12nt umi, both barcode and umi are in read1:
bcpat = '16B1[1:16] 12U1[17:28]'
```

## Cell calling parameters setting
There are two parameters that can significantly affect the final number of cells identified:
* min_umis (default value: 500): The larger of the value, the fewer cells will be called.
* max_adj_pvalue (default value: 0.01): The smaller of the value, the fewer cells will be called.
* These 2 parameters can be setted when calling the 'hellobc_pipeline' function:
```
# For example, here we set the min_umis=100, max_adj_pvalue=0.01
hellobc_pipeline(workst, sampid, fq1, fq2, bcpat, gtf_dir, star_index_dir, min_umis=100, max_adj_pvalue=0.01)
```