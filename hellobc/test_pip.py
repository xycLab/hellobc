from hellobc.pipline.runpip import xu_pipeline

workst = "/mnt/sda/xxt/workst/sevci/20240517"
sampid = "ZT-135-ZT-135-1"
# fq1 = "/mnt/sda/xxt/data/seqdata/biozt/221307CSB_HCCO3/221307CSB_HCCO3_R1.fastq.gz"
# fq2 = "/mnt/sda/xxt/data/seqdata/biozt/221307CSB_HCCO3/221307CSB_HCCO3_R2.fastq.gz"
fq1 = "/mnt/sda/xxt/data/seqdata/biozt/ZT-135-ZT-135-1/ZT-135-ZT-135-1_R1.fastq.gz"
fq2 = "/mnt/sda/xxt/data/seqdata/biozt/ZT-135-ZT-135-1/ZT-135-ZT-135-1_R2.fastq.gz"
bcpat = '36B1[1:12; 19:30; 37:48] 8U1[49:56]'                                                       # biozt
# bcpat = '16B1[1:16] 12U1[17:28]'                                                                    # 10x chem_v2

gtf_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
star_index_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-GRCh38-2020-A/star"

mouse_gtf_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-mm10-2020-A/genes/genes.gtf"
mouse_star_index_dir = "/mnt/sda/xxt/data/seqdata/ref/refdata-gex-mm10-2020-A/star"

star_path = '/mnt/sda/xxt/soft/STAR-2.6.1a/bin/Linux_x86_64_static/STAR'
fc_path = '/mnt/sda/xxt/soft/subread-2.0.6-Linux-x86_64/bin/featureCounts'

xu_pipeline(workst, sampid, fq1, fq2, bcpat, gtf_dir, star_index_dir, init_method='ordmag', min_umis=500, max_adj_pvalue=0.01, star_path=star_path, fc_path=fc_path, corenum=-1)
