from hellobc.demulti.parsebc import run_parseFastq
from hellobc.utils.sysop import run_mapping
from hellobc.getcount.bcgenecount import getCountMatrix

MIN_UMIS = 500
MAX_ADJ_PVALUE = 0.01                                   # p-value threshold

def xu_pipeline(workst:str, sampid:str, fq1:str, fq2:str, bcpattern:str, gtf:str, staridx:str, 
                min_umis=MIN_UMIS, max_adj_pvalue=MAX_ADJ_PVALUE, init_method:str='ordmag', star_path:str='STAR', fc_path:str='featureCounts', corenum:int=-1, remaintmp:bool=False):
    run_parseFastq(workst=workst, sampid=sampid, fq1=fq1, fq2=fq2, bcpat=bcpattern, corenum=corenum, remain_tmp=remaintmp)
    run_mapping(workst=workst, sampid=sampid, corenum=corenum, in_fq2=f"{workst}/{sampid}/01_parsed/bcparsed_R2.fastq", in_fq1='',
                  gtf_dir=gtf, star_index_dir=staridx, star_path=star_path, fc_path=fc_path)
    getCountMatrix(workst=workst, sampid=sampid, gtf_path=gtf, init_method=init_method, remain_tmp=remaintmp, min_umis=min_umis, max_adj_pvalue=max_adj_pvalue,
                   corenum=corenum, self_anno=True, get_raw=True, countbam=True, get_filtered=True)

