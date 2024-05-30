import multiprocessing
from collections import Counter
from Bio.SeqRecord import SeqRecord
from .parsebc import *
from hellobc.fileprocess.fastq_process import *
import hellobc.utils as ut

samp_n_reads = 0
def slipFilesCallBack(arg_tup:tuple) -> dict:
    sbc_seq = arg_tup[0]
    rec = arg_tup[1]
    recoth = arg_tup[2]
    outind = arg_tup[3]
    outpath = arg_tup[4]
    if outind == '1':
        othind = '2'
    else:
        othind = '1'

    sampfname = f"{outpath}/split_fastqs/{sbc_seq}_R{outind}.fastq"
    sampfname_another = f"{outpath}/split_fastqs/{sbc_seq}_R{othind}.fastq"
    sbclistfname = f"{outpath}/sbc_list.txt"

    global sbc_stat_dict
    
    if sbc_seq in sbc_stat_dict:
        sbc_stat_dict[sbc_seq] += 1
    else:
        sbc_stat_dict[sbc_seq] = 1
        ut.fileop.addOneLine(fname=sbclistfname, theline=sbc_seq)
        
    outfqobj = fastqFile(fname=sampfname)
    outothfqobj = fastqFile(fname=sampfname_another)
    outfqobj.appendOneRecord(oprecord=rec)
    outothfqobj.appendOneRecord(oprecord=recoth)

    del outfqobj
    del outothfqobj

    global samp_n_reads
    samp_n_reads += 1
    if samp_n_reads%10000 == 0:        
        print(f"parsed {samp_n_reads} reads.")

    return sbc_stat_dict


def makeSplitPath(splitpath:str, pnum:int) -> str:
    if pnum > 1:
        split_subpath = splitpath + "/tmp"
    else:
        split_subpath = splitpath
    
    ut.sysop.mkdir_p(dir=split_subpath)

    return split_subpath


def extractSampBarcode(bcobj:BcPattern, fq1:fastqFile, fq2:fastqFile, 
                       pid:int, pnum:int, splitsubpath:str='.', swtlist:str=None) -> dict:
    
    sbc_stat_dict = dict()
    
    if fq1.type == fqprocess.fqType.FASTQ:
        
        for record1, record2 in zip(SeqIO.parse(fq1.filename, "fastq"), 
                                    SeqIO.parse(fq2.filename, "fastq")):
    
            recoth = SeqRecord(Seq(""), id="", name="")
            
            if bcobj.bc_pattern_dict.sread == 1:
                s_record = record1
                recoth = record2
                outind = '1'
            else:
                s_record = record2
                recoth = record1
                outind = '2'

            rec = SeqRecord(Seq(""), id="", name="")
            quality_list = list()
            all_regions = list()
            
            sbc_seq = ''
            for region in bcobj.bc_pattern_dict.sregions:
                sbc_seq += str(s_record.seq[region[0]: region[1]+1])

            for r in bcobj.bc_pattern_dict.sregions:
                all_regions.append(r)

            all_regions = sorted(all_regions, key=lambda x: x[0])

            for i in range(len(all_regions)):
                if i < len(all_regions)-1:
                    rec += s_record.seq[all_regions[i][1]: all_regions[i+1][0]-1]
                    quality_list += s_record.letter_annotations["phred_quality"][all_regions[i][1]: all_regions[i+1][0]-1]
                else:
                    rec += s_record.seq[all_regions[i][1]:]
                    quality_list += s_record.letter_annotations["phred_quality"][all_regions[i][1]:]

            rec.name = s_record.name
            rec.id = s_record.id
            rec.description = s_record.description
            rec.letter_annotations["phred_quality"] = quality_list

            if swtlist is not None:
                sbc_seq = ''
            
            if outind == '1':
                othind = '2'
            else:
                othind = '1'

            if pnum > 1:                
                sampfname = f"{splitsubpath}/{sbc_seq}_{pid + 1}_R{outind}.{fq1.fqstr}"
                sampfname_another = f"{splitsubpath}/{sbc_seq}_{pid + 1}_R{othind}.{fq2.fqstr}"
                sbclistfname = f"{splitsubpath}/sbc_list_{pid}.txt"

            else:
                sampfname = f"{splitsubpath}/{sbc_seq}_R{outind}.{fq1.fqstr}"
                sampfname_another = f"{splitsubpath}/{sbc_seq}_R{othind}.{fq2.fqstr}"
                sbclistfname = f"{splitsubpath}/sbc_list.txt"

            if sbc_seq in sbc_stat_dict:
                sbc_stat_dict[sbc_seq] += 1
            else:
                sbc_stat_dict[sbc_seq] = 1
                ut.fileop.addOneLine(fname=sbclistfname, theline=sbc_seq)
                
            outfqobj = fastqFile(fname=sampfname)
            outothfqobj = fastqFile(fname=sampfname_another)
            outfqobj.appendOneRecord(oprecord=rec)
            outothfqobj.appendOneRecord(oprecord=recoth)

            del outfqobj
            del outothfqobj

    else:
        
        for record1, record2 in zip(SeqIO.parse(gzip.open(fq1.filename, 'rt'), "fastq"), 
                                    SeqIO.parse(gzip.open(fq2.filename, 'rt'), "fastq")):
    
            recoth = SeqRecord(Seq(""), id="", name="")

            if bcobj.bc_pattern_dict.sread == 1:
                s_record = record1
                recoth = record2
                outind = '1'
            else:
                s_record = record2
                recoth = record1
                outind = '2'

            rec = SeqRecord(Seq(""), id="", name="")
            quality_list = list()
            all_regions = list()
            
            sbc_seq = ''
            for region in bcobj.bc_pattern_dict.sregions:
                sbc_seq += str(s_record.seq[region[0]: region[1]+1])

            for r in bcobj.bc_pattern_dict.sregions:
                all_regions.append(r)

            all_regions = sorted(all_regions, key=lambda x: x[0])

            for i in range(len(all_regions)):
                if i < len(all_regions)-1:
                    rec += s_record.seq[all_regions[i][1]: all_regions[i+1][0]-1]
                    quality_list += s_record.letter_annotations["phred_quality"][all_regions[i][1]: all_regions[i+1][0]-1]
                else:
                    rec += s_record.seq[all_regions[i][1]:]
                    quality_list += s_record.letter_annotations["phred_quality"][all_regions[i][1]:]

            rec.name = s_record.name
            rec.id = s_record.id
            rec.description = s_record.description
            rec.letter_annotations["phred_quality"] = quality_list

            if swtlist is not None:
                sbc_seq = ''
            
            if outind == '1':
                othind = '2'
            else:
                othind = '1'

            if pnum > 1:                
                sampfname = f"{splitsubpath}/{sbc_seq}_{pid + 1}_R{outind}.{fq1.fqstr}"
                sampfname_another = f"{splitsubpath}/{sbc_seq}_{pid + 1}_R{othind}.{fq2.fqstr}"
                sbclistfname = f"{splitsubpath}/sbc_list_{pid}.txt"

            else:
                sampfname = f"{splitsubpath}/{sbc_seq}_R{outind}.{fq1.fqstr}"
                sampfname_another = f"{splitsubpath}/{sbc_seq}_R{othind}.{fq2.fqstr}"
                sbclistfname = f"{splitsubpath}/sbc_list.txt"

            if sbc_seq in sbc_stat_dict:
                sbc_stat_dict[sbc_seq] += 1
            else:
                sbc_stat_dict[sbc_seq] = 1
                ut.fileop.addOneLine(fname=sbclistfname, theline=sbc_seq)
                
            outfqobj = fastqFile(fname=sampfname)
            outothfqobj = fastqFile(fname=sampfname_another)
            outfqobj.appendOneRecord(oprecord=rec)
            outothfqobj.appendOneRecord(oprecord=recoth)

            del outfqobj
            del outothfqobj

    return sbc_stat_dict


def slipIntoSampFiles(bcobj:BcPattern, fq1path:str, fq2path:str, outdir:str, corenum:int=-1, swtlist:str=None) -> None:
    
    fq1 = fastqFile(fname=fq1path)
    fq2 = fastqFile(fname=fq2path)

    process_num = ut.sysop.get_process_num(corenum=corenum)
    split_subpath = makeSplitPath(splitpath=f"{outdir}/splitfiles", pnum=process_num)

    ut.fileop.splitFastqFile(rawfq=fq1.filename, nfiles=process_num, outpath=split_subpath, threads=20)
    ut.fileop.splitFastqFile(rawfq=fq2.filename, nfiles=process_num, outpath=split_subpath, threads=20)

    presults = list()
    pool = multiprocessing.Pool(processes=process_num)
    for pid in range(process_num):
        sfq1 = fastqFile(fname=f"{split_subpath}/{fq1.prefix}.part_{ut.fileop.zeroFill3(pid + 1)}.{fq1.tail}")
        sfq2 = fastqFile(fname=f"{split_subpath}/{fq2.prefix}.part_{ut.fileop.zeroFill3(pid + 1)}.{fq2.tail}")                
        presults.append(pool.apply_async(func=extractSampBarcode, args=(bcobj, sfq1, sfq2, pid, process_num, split_subpath)))
        
    pool.close()
    pool.join()

    pres_get_list = list()
    collect_dict = dict()
    collect_counter = Counter(collect_dict)
    for i_res in presults:
        pres_get_list.append(i_res.get())
        tmp_counter = Counter(i_res.get())
        collect_counter += tmp_counter
    collect_dict = dict(collect_counter)

    for sbc_seq in collect_dict.keys():
        merge_name_list1 = list()
        merge_name_list2 = list()

        for pid in range(process_num):            
            if sbc_seq in pres_get_list[pid]:
                pid_fname1 = f"{split_subpath}/{sbc_seq}_{pid + 1}_R1.{fq1.fqstr}"
                pid_fname2 = f"{split_subpath}/{sbc_seq}_{pid + 1}_R2.{fq2.fqstr}"
                merge_name_list1.append(pid_fname1)
                merge_name_list2.append(pid_fname2)
        
        ut.fileop.mergeFastqFile(fqlist=merge_name_list1, inprefix=f"{sbc_seq}_R1", outpath=f"{outdir}/splitfiles", ftype="fastq")
        ut.fileop.mergeFastqFile(fqlist=merge_name_list2, inprefix=f"{sbc_seq}_R2", outpath=f"{outdir}/splitfiles", ftype="fastq")
    
