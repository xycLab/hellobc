import pysam
import pandas as pd
import numpy as np
import multiprocessing
import shutil
import matplotlib.pyplot as plt
from hellobc.fileprocess.sam_process import bamLine
from .dropcalling import emptyDrops_drop_calling
from .dropcalling import countMatrix
from .bcrkplt import bc_curve_plt
from hellobc.utils.fileop import generate_gtf_anno
import hellobc.utils as ut

MIN_UMIS = 500
MAX_ADJ_PVALUE = 0.01                                   # p-value threshold

nreads = 0

def generateCountTsv(tsvfile:str, countdict:dict) -> str:
    
    with open(tsvfile, "w") as f:
        
        f.write("gene\tcell\tcount\n")
        dict1 = sorted(countdict.items(), key=lambda x: x[0])
        
        for gene, CBdict in dict1:
            dict2 = sorted(CBdict.items(), key=lambda x: x[0])
            for cb, umiset in dict2:
                f.write(f"{gene}\t{cb}\t{len(umiset)}\n")
    return tsvfile


def makeMtxPath(mtxpath:str, mtxtype:str, pnum:int) -> str:
    
    if mtxtype == 'raw':
        mtx_subpath = f"{mtxpath}/raw_feature_bc_matrix"
    elif mtxtype == 'filtered':
        mtx_subpath = f"{mtxpath}/filtered_feature_bc_matrix"
    else:
        raise ValueError("parameter mtxtype must be one of raw or filtered.")
    
    if pnum > 1:
        mtx_subpath += "/tmp"
    
    ut.sysop.mkdir_p(dir=mtx_subpath)

    return mtx_subpath


def writeFeatureTsv(mtx_type:str, mtx_subpath:str, annodf:pd.DataFrame, gene_list:list, gene_tag_type:str, gene_tag_type_ano:str,
                    init_bcs=np.array([]), eval_bcs=np.array([])) -> None:

    tsvpath = f"{mtx_subpath}/features.tsv"
    
    if mtx_type == 'filtered':
        with open(tsvpath, "w") as f:
            for gene_lab in gene_list:
                another_gene_lab = annodf.query(f'{gene_tag_type}=="{gene_lab}"')[gene_tag_type_ano].array[0]
                if gene_tag_type == 'gene':
                    f.write(f"{gene_lab}\t{another_gene_lab}\tGene Expression\n")
                else:
                    f.write(f"{another_gene_lab}\t{gene_lab}\tGene Expression\n")
    else:
        with open(tsvpath, "w") as f:
            for gene_lab in gene_list:
                another_gene_lab = annodf.query(f'{gene_tag_type}=="{gene_lab}"')[gene_tag_type_ano].array[0]
                if gene_tag_type == 'gene':
                    f.write(f"{gene_lab}\t{another_gene_lab}\tGene Expression\n")
                else:
                    f.write(f"{another_gene_lab}\t{gene_lab}\tGene Expression\n")


def writeBarcodeTsv(mtx_type:str, mtx_subpath:str, CB_list:list, lanenum:int, init_bcs=np.array([]), eval_bcs=np.array([])) -> list:
        
    tsvpath = f"{mtx_subpath}/barcodes.tsv"
    
    if mtx_type == 'filtered':
        fil_CB_list = list()
        with open(tsvpath, "w") as f:
            lane_num = lanenum
            for cb in CB_list:
                cb_lane = cb + f'-{lane_num}'
                if cb_lane in init_bcs or cb_lane in eval_bcs:
                    fil_CB_list.append(cb)
                    f.write(f"{cb}-{lane_num}\n")
        return fil_CB_list
    else:
        with open(tsvpath, "w") as f:
            lane_num = lanenum
            for cb in CB_list:
                f.write(f"{cb}-{lane_num}\n")
        return CB_list


def chunk_writeMatrixMtx(mtx_type:str, mtx_subpath:str, gene_list:list, CB_list:list, countlist:list, pid:int,
                   lanenum:int=1, init_bcs=np.array([]), eval_bcs=np.array([])):
    if mtx_type == 'filtered':
        mtxpath = f"{mtx_subpath}/tmp/filtered_matrix_{pid}.mtx"
        gene_dict = dict(zip(gene_list, list(range(1, len(gene_list)+1))))
        CB_dict = dict(zip(CB_list, list(range(1, len(CB_list)+1))))
        with open(mtxpath, "w") as f:        
            list1 = sorted(countlist)
            for gene, CBdict in list1:
                dict2 = sorted(CBdict.items(), key=lambda x: x[0])
                for cb, umiset in dict2:
                    cb_lane = cb + f'-{lanenum}'
                    if cb_lane in init_bcs or cb_lane in eval_bcs:
                        f.write(f"{gene_dict[gene]}\t{CB_dict[cb]}\t{len(umiset)}\n")
    return mtxpath
    

def writeMatrixMtx(mtx_type:str, mtx_subpath:str, gene_list:list, CB_list:list, countdict:dict, nlines:int,
                   lanenum:int=1, fil_pnum:int=1, fil_tmp:bool=False, init_bcs=np.array([]), eval_bcs=np.array([])) -> int:

    mtxpath = f"{mtx_subpath}/matrix.mtx"
    gene_dict = dict(zip(gene_list, list(range(1, len(gene_list)+1))))
    CB_dict = dict(zip(CB_list, list(range(1, len(CB_list)+1))))

    if mtx_type == 'filtered':
        ut.sysop.mkdir_p(dir=f"{mtx_subpath}/tmp")
        head_fname = f"{mtx_subpath}/tmp/filtered_head.mtx"
        with open(head_fname, "w") as f:
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write("%\n")
            f.write(f"{len(gene_list)}\t{len(init_bcs)+len(eval_bcs)}\t{nlines}\n") 
        
        count_list = list(countdict.items())
        countlist_list = list()
        presults = list()
        pool = multiprocessing.Pool(processes=fil_pnum)
        chunk_size = round(len(countdict)/fil_pnum)
        print(f"chunk_size = {chunk_size}")
        for fil_pid in range(fil_pnum):
            p_sidx = fil_pid * chunk_size
            if fil_pid == fil_pnum - 1:
                p_eidx = -1
            else:
                p_eidx = (fil_pid + 1) * chunk_size
            print(f"chunk of No.{fil_pid}: {p_sidx} : {p_eidx-1}")
            countlist_list.append(count_list[p_sidx:p_eidx])
            
            presults.append(pool.apply_async(func=chunk_writeMatrixMtx, args=(mtx_type, mtx_subpath, gene_list, CB_list, countlist_list[fil_pid], fil_pid, lanenum, init_bcs, eval_bcs)))
        
        pool.close()
        pool.join()

        pres_get_list = list()
        pres_get_list.append(head_fname)
        for i_res in presults:
            pres_get_list.append(i_res.get())

        ut.fileop.mergePlaneFiles(namelist=pres_get_list, outpath=mtx_subpath, outprefix="matrix", ftype="mtx")
        if not fil_tmp:
            shutil.rmtree(path=f"{mtx_subpath}/tmp")
            pass
        
    else:
        with open(mtxpath, "w") as f:        
            dict1 = sorted(countdict.items(), key=lambda x: x[0])            
            f.write("%%MatrixMarket matrix coordinate integer general\n")
            f.write("%\n")
            #f.write("\n")
            f.write(f"{len(gene_list)}\t{len(CB_list)}\t{nlines}\n")
            for gene, CBdict in dict1:
                dict2 = sorted(CBdict.items(), key=lambda x: x[0])
                for cb, umiset in dict2:
                    f.write(f"{gene_dict[gene]}\t{CB_dict[cb]}\t{len(umiset)}\n")
    return nlines


def count_fil_nlines(countlist:list, init_bcs, eval_bcs, lanenum, fil_pid:int) -> int:

    nlines_fil = 0
    print(f"p No.{fil_pid} started......")
    for tup in countlist:
        for cb in tup[1].keys():
            cb_lane = cb + f'-{lanenum}'
            if cb_lane in init_bcs or cb_lane in eval_bcs:
                nlines_fil += 1
    return nlines_fil


def generate_legacy_matrix(workst:str, sampid:str, annofile:str, gene_list:list, CB_list:list, countdict:dict, nlines:int, lanenum:int=1, min_umis=MIN_UMIS, max_adj_pvalue=MAX_ADJ_PVALUE,
                           init_method:str='ordmag', fil_pnum:int=1, fil_tmp:bool=False, get_raw:bool=True, get_filtered:bool=True, bcplt:bool=True):

    mtxpath = f"{workst}/{sampid}/03_count"
    print("Start generating filtered mtx...")

    annodf = pd.read_csv(annofile, header=0, sep='\t')
    ut.sysop.mkdir_p(f"{mtxpath}/raw_feature_bc_matrix")
    ut.sysop.mkdir_p(f"{mtxpath}/filtered_feature_bc_matrix")
    raw_mtx_subpath = makeMtxPath(mtxpath=mtxpath, mtxtype='raw', pnum=1)
    filtered_mtx_subpath = makeMtxPath(mtxpath=mtxpath, mtxtype='filtered', pnum=1)

    if annodf['gene'].isin([gene_list[0]]).any():
        gene_tag_type = 'gene'
        gene_tag_type_ano = 'gene_name'
    else:
        gene_tag_type = 'gene_name'
        gene_tag_type_ano = 'gene'

    # generate raw mtx
    if get_raw:
        writeFeatureTsv(mtx_type='raw', mtx_subpath=raw_mtx_subpath, annodf=annodf, gene_list=gene_list, gene_tag_type=gene_tag_type,
                        gene_tag_type_ano=gene_tag_type_ano)
        writeBarcodeTsv(mtx_type='raw', mtx_subpath=raw_mtx_subpath, CB_list=CB_list, lanenum=lanenum)
        writeMatrixMtx(mtx_type='raw', mtx_subpath=raw_mtx_subpath, gene_list=gene_list, CB_list=CB_list, countdict=countdict, nlines=nlines, lanenum=lanenum)
        print(f"in raw mtx, nlines = {nlines}")
    
        print("Done!")
    raw_mtx_path = f"{mtxpath}/raw_feature_bc_matrix"

    # -------------------------- drop calling ---------------------------- #    
    raw_mtx = countMatrix.readLegMatrix(f_mtx = raw_mtx_path)

    init, eval = emptyDrops_drop_calling(raw_mtx=raw_mtx, init_method=init_method, min_umis=min_umis, max_adj_pvalue=max_adj_pvalue)
    ut.sysop.mkdir_p(dir=f"{mtxpath}/plt")
    np.save(f"{mtxpath}/plt/init_idx.npy", init)
    np.save(f"{mtxpath}/plt/eval_idx.npy", eval)
    print(f"---------------------------------------\nCalled cells: {len(init)+len(eval)}\n---------------------------------------")
    init_bcs = raw_mtx.bcs[init]
    eval_bcs = raw_mtx.bcs[eval]

    # generate filtered mtx
    if get_filtered:
        count_list = list(countdict.items())
        countlist_list = list()
        presults = list()
        pool = multiprocessing.Pool(processes=fil_pnum)
        chunk_size = round(len(countdict)/fil_pnum)
        print(f"chunk_size = {chunk_size}")
        for fil_pid in range(fil_pnum):
            p_sidx = fil_pid * chunk_size
            if fil_pid == fil_pnum - 1:
                p_eidx = -1
            else:
                p_eidx = (fil_pid + 1) * chunk_size
            print(f"chunk of No.{fil_pid}: {p_sidx} : {p_eidx}")
            countlist_list.append(count_list[p_sidx:p_eidx])
            
            presults.append(pool.apply_async(func=count_fil_nlines, args=(countlist_list[fil_pid], init_bcs, eval_bcs, lanenum, fil_pid)))
        
        pool.close()
        pool.join()

        pres_get_list = list()
        for i_res in presults:
            pres_get_list.append(i_res.get())
        nlines_fil = sum(pres_get_list)

        print("Writing filtered fearutes")
        writeFeatureTsv(mtx_type='filtered', mtx_subpath=filtered_mtx_subpath, annodf=annodf, gene_list=gene_list, gene_tag_type=gene_tag_type,
                        gene_tag_type_ano=gene_tag_type_ano, init_bcs=init_bcs, eval_bcs=eval_bcs)
        print("Writing filtered barcodes")
        fil_CB_list = writeBarcodeTsv(mtx_type='filtered', mtx_subpath=filtered_mtx_subpath, CB_list=CB_list, lanenum=lanenum, 
                                      init_bcs=init_bcs, eval_bcs=eval_bcs)
        print("Writing filtered mtx")
        writeMatrixMtx(mtx_type='filtered', mtx_subpath=filtered_mtx_subpath, gene_list=gene_list, CB_list=fil_CB_list, countdict=countdict, nlines=nlines_fil,
                       lanenum=lanenum, fil_pnum=fil_pnum, fil_tmp=fil_tmp, init_bcs=init_bcs, eval_bcs=eval_bcs)
           
    filtered_mtx_path = f"{mtxpath}/filtered_feature_bc_matrix"
    
    print("Done!")
    
    if bcplt:
        print("Plotting barcode rank plot...")
        # bc_rk_plt(raw_mtx=raw_mtx, init_bcs=init, eval_bcs=eval, fig_path=f"{workst}/{sampid}/03_count")
        bc_curve_plt(raw_mtx=raw_mtx, init_idx=init, eval_idx=eval, fig_path=f"{workst}/{sampid}/03_count")
        print("Plotting done!")
    return raw_mtx_path, filtered_mtx_path


def bc_rk_plt(raw_mtx:countMatrix, init_bcs, eval_bcs, fig_path:str='.'):

    my_bc_colors = {'blue_1': '#08336E',
                    'blue_2': '#105C9C',
                    'blue_3': '#3888C0',
                    'blue_3': '#68ACCD',
                    'blue_4': '#AAD7E5',
                    'blue_5': '#D2E3F3',
                    'blue_6': '#F4F9FE'}

    bc_colors = np.array([my_bc_colors['blue_5'] for _ in range(len(raw_mtx.bcs))])
    bc_colors[init_bcs] = my_bc_colors['blue_1']
    bc_colors[eval_bcs] = my_bc_colors['blue_2']
    ind_plt = np.argsort(raw_mtx.counts_per_bc)[::-1]                       

    fig = plt.figure()
    fig1 = fig.add_subplot(111)

    fig1.scatter(x=range(1, len(raw_mtx.counts_per_bc)+1), y=raw_mtx.counts_per_bc[ind_plt],
                 c=bc_colors[ind_plt], s=7, linewidths=0, zorder=10)
    fig1.loglog()
    fig1.set_xlim(0, len(raw_mtx.counts_per_bc)*1.25)
    fig1.set_xlabel('Barcode index')
    fig1.set_ylabel('Count')
    fig1.spines['bottom'].set_linewidth(0.5)
    fig1.spines['left'].set_linewidth(0.5)
    fig1.spines['top'].set_visible(False)
    fig1.spines['right'].set_visible(False)
    fig1.set_xticks([], minor=True)
    fig1.set_yticks([], minor=True)
    fig1.tick_params(bottom=False, top=False, left=False, right=False)
    fig1.grid(color="#F0F0F2", zorder=0)

    fig.savefig(f"{fig_path}/bcs_rk_plt.png", bbox_inches='tight', dpi=400, transparent=True)


def count_bam(bamfile:str):

    # error correct settings
    # umi_corr_flag = False
    # bc_corr_flag = False
    
    bam_file = pysam.AlignmentFile(bamfile, "rb")
    bam_iter = bam_file.fetch()

    count_dict = dict()
    count_dict_bc = dict()
    count_dict_only_bc = dict()
    gene_set = set()
    CB_set = set()
    
    global nreads
    nlines = 0
   
    for read in bam_iter:
        bl = bamLine(fname=bamfile, read=read)        
        nreads += 1
        ut.fileop.infoPasedReads(nreads=nreads, unit=1000000)

        if bl.XT == '':
            continue
        else:
            
            if bl.XT not in count_dict:
                count_dict[bl.XT] = {}
            if bl.CB not in count_dict[bl.XT]:
                nlines += 1                
                count_dict[bl.XT][bl.CB] = set()
                count_dict[bl.XT][bl.CB].add(bl.umi)
            else:
                count_dict[bl.XT][bl.CB].add(bl.umi)

            gene_set.add(bl.XT)
            CB_set.add(bl.CB)
                
    gene_list = list(gene_set)
    gene_list.sort()
    CB_list = list(CB_set)
    CB_list.sort()
    
    return nlines, count_dict, gene_list, CB_list, count_dict_bc


def getCountMatrix(workst:str, sampid:str, gtf_path:str, lanenum:int=1, init_method:str='ordmag', corenum:int=-1, min_umis=MIN_UMIS, max_adj_pvalue=MAX_ADJ_PVALUE, 
                   self_anno:bool=False, remain_tmp:bool=False, countbam:bool=True, get_raw:bool=True, get_filtered:bool=True, bcplt:bool=True) -> dict:
    
    bamfile = f"{workst}/{sampid}/02_mapping/assigned_sorted.bam"
    if self_anno:
        anno_prefix = generate_gtf_anno(workst=workst, sampid=sampid, f_gtf=gtf_path)
    else:
        anno_prefix = ''

    if countbam:
        nlines, count_dict, gene_list, CB_list, count_dict_bc = count_bam(bamfile=bamfile)
        # run_adj(count_dict_bc=count_dict_bc)
    else:
        nlines = 0
        count_dict = dict()
        gene_list = list()
        CB_list = list()
    
    outpath = f"{workst}/{sampid}/03_count"    
    ut.sysop.mkdir_p(dir=outpath)
    fil_pnum = ut.sysop.get_process_num(corenum)
    generate_legacy_matrix(workst=workst, sampid=sampid, annofile=f"{workst}/{sampid}/02_mapping/{anno_prefix}.self.anno.txt",
                           gene_list=gene_list, CB_list=CB_list, countdict=count_dict, nlines=nlines, lanenum=lanenum, fil_pnum=fil_pnum, fil_tmp=remain_tmp, init_method=init_method,
                           min_umis=min_umis, max_adj_pvalue=max_adj_pvalue, get_raw=get_raw, get_filtered=get_filtered, bcplt=bcplt)
    
    return count_dict

