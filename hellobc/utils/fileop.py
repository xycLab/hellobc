import numpy as np
import os
import pyranges as pr
from . import sysop

# ------------------------ basic file operations ---------------------- #
def infoPasedReads(nreads:int, unit:int) -> None:
    if nreads%unit == 0:
        print(f"parsed {nreads} reads.")


def addOneLine(fname:str, theline:str) -> None:
    if theline[-1] != '\n':
        theline += '\n'
    with open(fname, "a") as f:
        f.write(theline)
        

def splitList2N(rawlist:list, splitnum:int) -> list:
    arraytmp = np.array(rawlist)
    return np.array_split(arraytmp, splitnum)


def splitDictbyKey(rawdict:dict, subkeys:list):
    subdict = {key: rawdict[key] for key in subkeys}
    return subdict


def mergePlaneFiles(ftype:str, outpath:str='.', outprefix:str='', nfiles:int=0, infolder:str='.', inprefix:str='merged',  namelist:list=list()) -> None:
    
    os_command_str = 'cat'
    if outprefix == '':
        outprefix = inprefix
    
    if len(namelist) == 0:
        for i in range(nfiles):
            fname = f"{infolder}/{inprefix}_{i}.{ftype}"       
            os_command_str += f" {fname}"
                
        os_command_str += f" > {outpath}/{outprefix}.{ftype}"

        os.system(command=os_command_str)
    else:
        os_command_str += ' '
        os_command_str += ' '.join(namelist)
    
        os_command_str += f" > {outpath}/{outprefix}.{ftype}"
        # print(os_command_str)

        os.system(command=os_command_str)


def splitFastqFile(rawfq:str, nfiles:int, outpath:str, threads:int) -> None:
    os_command_str = 'seqkit split2'
    #os_command_str += " --force"
    os_command_str += f" --by-part {nfiles}"
    os_command_str += f" --threads {threads}"
    os_command_str += f" {rawfq}"
    os_command_str += f" --out-dir {outpath}"
    os.system(command=os_command_str)


def mergeFastqFile(fqlist:list, inprefix:str, outpath:str, ftype:str) -> None:
    
    os_command_str = 'cat'
    os_command_str += ' '
    os_command_str += ' '.join(fqlist)
    
    os_command_str += f" > {outpath}/{inprefix}.{ftype}"

    os.system(command=os_command_str)


def zeroFill3(num):
 
    return str(num).zfill(3)


# --------------------------- gtf file parse ---------------------------- #
def generate_gtf_anno(workst:str, sampid:str, f_gtf:str, anno_prefix:str='') -> str:

    # file path                                            
    out_file = f"{workst}/{sampid}/02_mapping/{anno_prefix}.self.anno.txt"
    # sysop.mkdir_p(dir=f"{workst}/refanno")
    my_gtf = pr.read_gtf(f_gtf)                                         
    df_tmp = my_gtf.df[['gene_name', 'gene_id']]
    df_tmp = df_tmp.drop_duplicates()
    gene_name = df_tmp['gene_name']
    gene = df_tmp['gene_id']

    with open(out_file, "w") as f:
        f.write("gene\tgene_name\n")
        for g, g_name in zip(gene, gene_name):
            f.write(f"{g}\t{g_name}\n")

    print("Parse Gtf Done!")
    return anno_prefix
