import os
import multiprocessing

def mkdir_p(dir:str) -> None:
    if not os.path.exists(dir):
        os.makedirs(dir)

def get_process_num(corenum:int) -> int:
    cpu_num = multiprocessing.cpu_count()
    if corenum == -1:
        process_num = cpu_num - 2
    elif corenum > cpu_num:
        raise ValueError("parameter corenum must less than or equal to cpu cores.")
    elif corenum < 1:
        raise ValueError("at least one cpu core should be used.")
    else:
        process_num = corenum

    return process_num

# ------------------------------ command line ----------------------------- #
def run_STAR(threads:int, star_index_dir:str, in_fq1:str, in_fq2:str, star_path:str='STAR') -> str:

    command = f"""
    {star_path} \\
    --runThreadN {threads} \\
    --genomeDir {star_index_dir} \\
    --readFilesIn {in_fq1} {in_fq2} \\
    --outFilterMultimapNmax 1 \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix STAR. \\
    --outSAMunmapped Within"""

    os.system(command=command)
    return command


def run_featureCounts(threads:int, gtf_dir:str, name_or_id:str="gene_id", fc_path:str='featureCounts') -> str:
    
    command = f"""
    {fc_path} -a {gtf_dir} \\
    -g {name_or_id} \\
    -t gene \\
    -o gene_assigned_name \\
    -R BAM \\
    STAR.Aligned.sortedByCoord.out.bam \\
    -T {threads}"""

    os.system(command=command)
    return command


def run_samtools(threads:int) -> list:

    command1 = f"""
    samtools sort \\
    -@ {threads} \\
    -o assigned_sorted.bam \\
    -O BAM STAR.Aligned.sortedByCoord.out.bam.featureCounts.bam"""
    os.system(command=command1)
    command2 = f"""
    samtools index \\
    -@ {threads} \\
    assigned_sorted.bam"""
    os.system(command=command2)

    command3 = f"""

    """

    return [command1, command2]


def run_mapping(workst:str, sampid:str, corenum:int, in_fq1:str, in_fq2:str,
                gtf_dir:int, star_index_dir:str, name_or_id:str="gene_id", star_path:str='STAR', fc_path:str='featureCounts',
                STAR:bool=True, featureCounts:bool=True, samtools:bool=True) -> None:
    
    threads = get_process_num(corenum)
    fc_threads = threads
    if threads > 64:
        fc_threads = 64
    elif threads < 1:
        fc_threads = 1
    mkdir_p(dir=f"{workst}/{sampid}/02_mapping")
    os.chdir(path=f"{workst}/{sampid}/02_mapping")
    # os.system(command="source activate upstream")
    if STAR:
        run_STAR(threads=threads, star_index_dir=star_index_dir, in_fq1=in_fq1, in_fq2=in_fq2, star_path=star_path)
    if featureCounts:
        run_featureCounts(threads=fc_threads, gtf_dir=gtf_dir, name_or_id=name_or_id, fc_path=fc_path)
    if samtools:
        run_samtools(threads=threads)
    os.chdir(path=f"{workst}/{sampid}")


def run_pipeline():

    # parsebc.py: extract fq file
    # run_mappint()
    # getcount.bcgenecount.py  ->  raw_mtx                   
    # getcount.dropcalling.py  ->  filtered_mtx
    # getcount.bcrankplot.py   ->  Barcode Rank Plot
    pass
