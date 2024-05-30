import shutil
import datetime
import multiprocessing
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hellobc.fileprocess import fastq_process
import hellobc.fileprocess.fastq_process as fqprocess
from hellobc.fileprocess.fastq_process import fastqFile
import hellobc.utils as ut

class BcPatternDict():
    '''
    Class for storing the parsed barcode patterns.
    '''
    def __init__(self, in_dict:dict) -> None:
        self.raw_dict = in_dict
        if 'S' in in_dict or 's' in in_dict:
            self.sread, self.slen, self.sregions = self._parseBcDict(bctype='S')
        self.bread, self.blen, self.bregions = self._parseBcDict(bctype='B')
        self.uread, self.ulen, self.uregions = self._parseBcDict(bctype='U')

    def _parseBcDict(self, bctype:str):
        p_tmp = self.raw_dict[bctype]
        type_ind = p_tmp.find(bctype)

        type_len = int(p_tmp[:type_ind])                            

        region_str = p_tmp.split('[')[-1].replace(']', '')
        region_str = region_str.split(';')
       
        region_num = p_tmp.count(';') + 1
        type_regions = list(range(region_num)) 

        for i in range(region_num):
            r_tmp = region_str[i].split(':')
            type_regions[i] = tuple([int(r_tmp[0]), int(r_tmp[-1])])    
        
        type_read = int(p_tmp[type_ind+1])                              
        
        return type_read, type_len, type_regions


class BcPattern():
    '''
    Classes for parsing and storing barcode patterns.
    '''    
    def __init__(self, in_pattern:str) -> None:
        self.raw_pattern = in_pattern                                                           
        self.bc_pattern = ''.join(self.raw_pattern.split())                                     
        self.depart_index = [i for i, char in enumerate(self.bc_pattern) if char == ']']        
        self.bc_pattern_dict, self.bc_pattern_norm_dict = self._getPatternDict()        

    def _getPatternDict(self) -> BcPatternDict:
        bcpt_dict = {}
        split_bcstr = list()
        ind2 = list(range(len(self.depart_index)))
        for i in range(len(self.depart_index)):
            if i == 0:
                ind2[i] = 0
            else:
                ind2[i] = self.depart_index[i-1] + 1

        for j, k in zip(ind2, self.depart_index):
            split_bcstr.append(self.bc_pattern[j:k+1])
        
        for sk in split_bcstr:
            if 'S' in sk or 's' in sk:
                bcpt_dict['S'] = sk
            elif 'B' in sk or 'b' in sk:
                bcpt_dict['B'] = sk
            elif 'U' in sk or 'u' in sk:
                bcpt_dict['U'] = sk
            else:
                print('Input -bc-pattern ERROR: need one of S B U')

        bc_dict_obj = BcPatternDict(in_dict=bcpt_dict)
        
        return bc_dict_obj, bcpt_dict


def getExtractRegions(bc_pattern:BcPattern, record1:SeqRecord, record2:SeqRecord):
    '''
    Extracts each type of regions of the barcode from the sequence, and returns the extracted sequence and the barcode.
    '''
    seq_name_tail = ''
    barcode_str = ''
    umi_str = ''
    all_regions_1 = list()
    all_regions_2 = list()
    quality_list1 = list()
    quality_list2 = list()
    rec1 = SeqRecord(Seq(""), id="", name="")
    rec2 = SeqRecord(Seq(""), id="", name="")
    
    if 'S' in bc_pattern.bc_pattern_norm_dict:
        if bc_pattern.bc_pattern_dict.sread == 1:
            s_record = record1
        else:
            s_record = record2
        for region in bc_pattern.bc_pattern_dict.sregions:
                seq_name_tail += s_record.seq[region[0]-1: region[1]]

    if bc_pattern.bc_pattern_dict.bread == 1:
        b_record = record1
    else:
        b_record = record2

    if bc_pattern.bc_pattern_dict.uread == 1:
        u_record = record1
    else:
        u_record = record2

    for region in bc_pattern.bc_pattern_dict.bregions:
        barcode_str += b_record.seq[region[0]-1: region[1]]    
    
    for region in bc_pattern.bc_pattern_dict.uregions:
        umi_str += u_record.seq[region[0]-1: region[1]]
    seq_name_tail = barcode_str + '_' + umi_str
    umi_str = str(umi_str)
    barcode_str = str(barcode_str)

    if 'S' in bc_pattern.bc_pattern_norm_dict:
        if bc_pattern.bc_pattern_dict.sread == 1:
            for r in bc_pattern.bc_pattern_dict.sregions:
                all_regions_1.append(r)
        else:
            for r in bc_pattern.bc_pattern_dict.sregions:
                all_regions_2.append(r)

    if bc_pattern.bc_pattern_dict.bread == 1:
        for r in bc_pattern.bc_pattern_dict.bregions:
            all_regions_1.append(r)
    else:
        for r in bc_pattern.bc_pattern_dict.bregions:
            all_regions_2.append(r)

    if bc_pattern.bc_pattern_dict.uread == 1:
        for r in bc_pattern.bc_pattern_dict.uregions:
            all_regions_1.append(r)
    else:
        for r in bc_pattern.bc_pattern_dict.uregions:
            all_regions_2.append(r)
    
    all_regions_1 = sorted(all_regions_1, key=lambda x: x[0])
    all_regions_2 = sorted(all_regions_2, key=lambda x: x[0])

    if len(all_regions_1) > 0:            
        for i in range(len(all_regions_1)):
            if i < len(all_regions_1)-1:
                rec1 += record1.seq[all_regions_1[i][1]: all_regions_1[i+1][0]-1]
                quality_list1 += record1.letter_annotations["phred_quality"][all_regions_1[i][1]: all_regions_1[i+1][0]-1]
            else:
                rec1 += record1.seq[all_regions_1[i][1]:]
                quality_list1 += record1.letter_annotations["phred_quality"][all_regions_1[i][1]:]
    else:
        rec1 += record1.seq
        quality_list1 += record1.letter_annotations["phred_quality"]

    if len(all_regions_2) > 0:
        for i in range(len(all_regions_2)):
            if i < len(all_regions_2)-1:
                rec2 += record2.seq[all_regions_2[i][1]: all_regions_2[i+1][0]-1]
                quality_list2 += record2.letter_annotations["phred_quality"][all_regions_2[i][1]: all_regions_2[i+1][0]-1]
            else:
                rec2 += record2.seq[all_regions_2[i][1]:]
                quality_list2 += record2.letter_annotations["phred_quality"][all_regions_2[i][1]:]
    else:
        rec2 += record2.seq
        quality_list2 += record2.letter_annotations["phred_quality"]

    return rec1, rec2, quality_list1, quality_list2, seq_name_tail, barcode_str, umi_str

    
def processSeqRecord(parse_subpath:str, bc_pattern:BcPattern, fq1:fastqFile, fq2:fastqFile) -> tuple:
    '''
    Barcode extraction of sequencing files.
    '''    
    outfq1 = fastqFile(fname=f"{parse_subpath}/{fq1.prefix}_parsedbc.{fq1.fqstr}")
    outfq2 = fastqFile(fname=f"{parse_subpath}/{fq2.prefix}_parsedbc.{fq2.fqstr}")
    unmapped_count_dict = dict()
    
    if fq1.type == fqprocess.fqType.FASTQ:    
        for record1, record2 in zip(SeqIO.parse(fq1.filename, "fastq"), 
                                    SeqIO.parse(fq2.filename, "fastq")):
                    
            rec1, rec2, quality_list1, quality_list2, seq_name_tail, barcode_str, umi_str = getExtractRegions(bc_pattern=bc_pattern, record1=record1, record2=record2) 
            # generate unmapped dict
            if barcode_str not in unmapped_count_dict:
                unmapped_count_dict[barcode_str] = set()
                unmapped_count_dict[barcode_str].add(umi_str)
            else:
                unmapped_count_dict[barcode_str].add(umi_str)
            # parse record
            rec1_discri_tail = ''
            if len(record1.description.split(' ')) == 1:
                rec1_discri_tail = ''
            else:
                rec1_discri_tail = str(record1.description.split(' ')[1])

            rec2_discri_tail = ''
            if len(record2.description.split(' ')) == 1:
                rec2_discri_tail = ''
            else:
                rec2_discri_tail = str(record2.description.split(' ')[1])
            
            rec1.name = str(record1.name) + '_' + str(seq_name_tail)
            rec1.id = str(record1.id) + '_' + str(seq_name_tail)
            rec1.description = str(rec1.id) + ' ' + rec1_discri_tail
            rec1.letter_annotations["phred_quality"] = quality_list1

            rec2.name = str(record2.name) + '_' + str(seq_name_tail)
            rec2.id = str(record2.id) + '_' + str(seq_name_tail)
            rec2.description = str(rec2.id) + ' ' + rec2_discri_tail
            rec2.letter_annotations["phred_quality"] = quality_list2
            
            outfq1.appendOneRecord(oprecord=rec1)
            outfq2.appendOneRecord(oprecord=rec2)

    else:
        for record1, record2 in zip(SeqIO.parse(gzip.open(fq1.filename, 'rt'), "fastq"), 
                                    SeqIO.parse(gzip.open(fq2.filename, 'rt'), "fastq")):                    

            rec1, rec2, quality_list1, quality_list2, seq_name_tail, barcode_str, umi_str = getExtractRegions(bc_pattern=bc_pattern, record1=record1, record2=record2)
            # generate unmapped dict
            if barcode_str not in unmapped_count_dict:
                unmapped_count_dict[barcode_str] = set()
                unmapped_count_dict[barcode_str].add(umi_str)
            else:
                unmapped_count_dict[barcode_str].add(umi_str)
            # parse record            
            rec1_discri_tail = ''
            if len(record1.description.split(' ')) == 1:
                rec1_discri_tail = ''
            else:
                rec1_discri_tail = str(record1.description.split(' ')[1])

            rec2_discri_tail = ''
            if len(record2.description.split(' ')) == 1:
                rec2_discri_tail = ''
            else:
                rec2_discri_tail = str(record2.description.split(' ')[1])
            
            rec1.name = str(record1.name) + '_' + str(seq_name_tail)
            rec1.id = str(record1.id) + '_' + str(seq_name_tail)
            rec1.description = str(rec1.id) + ' ' + rec1_discri_tail
            rec1.letter_annotations["phred_quality"] = quality_list1

            rec2.name = str(record2.name) + '_' + str(seq_name_tail)
            rec2.id = str(record2.id) + '_' + str(seq_name_tail)
            rec2.description = str(rec2.id) + ' ' + rec2_discri_tail
            rec2.letter_annotations["phred_quality"] = quality_list2
            
            outfq1.appendOneRecord(oprecord=rec1)
            outfq2.appendOneRecord(oprecord=rec2)
    
    return (outfq1, outfq2, unmapped_count_dict)
    
n_reads = 0


def makeParsePath(parsepath:str, pnum:int) -> str:
    if pnum > 1:
        parse_subpath = parsepath + "/tmp"
    else:
        parse_subpath = parsepath
    
    ut.sysop.mkdir_p(dir=parse_subpath)

    return parse_subpath


def parseBarcode(fq1path:str, fq2path:str, bcpattern:BcPattern, parsepath:str, corenum:int=-1) -> None:
    '''
    Parallel processing barcode extraction.
    '''
    fq1 = fastqFile(fname=fq1path)
    fq2 = fastqFile(fname=fq2path)
    
    process_num = ut.sysop.get_process_num(corenum=corenum)
    parse_subpath = makeParsePath(parsepath=parsepath, pnum=process_num)
    presults = list()
    ofq1_list = list()
    ofq2_list = list()
    unmapped_dic_list = list()
    
    if process_num > 1:
        ut.fileop.splitFastqFile(rawfq=fq1.filename, nfiles=process_num, outpath=parse_subpath, threads=process_num)
        ut.fileop.splitFastqFile(rawfq=fq2.filename, nfiles=process_num, outpath=parse_subpath, threads=process_num)

        pool = multiprocessing.Pool(processes=process_num)
        for pid in range(process_num):
            sfq1 = fastqFile(fname=f"{parse_subpath}/{fq1.prefix}.part_{ut.fileop.zeroFill3(pid + 1)}.{fq1.tail}")
            sfq2 = fastqFile(fname=f"{parse_subpath}/{fq2.prefix}.part_{ut.fileop.zeroFill3(pid + 1)}.{fq2.tail}")
            # print(sfq1.filename)
            # print(sfq2.filename)
            
            presults.append(pool.apply_async(func=processSeqRecord, args=(parse_subpath, bcpattern, sfq1, sfq2)))
        
        pool.close()
        pool.join()

        pres_get_list = list()
        for i_res in presults:
            pres_get_list.append(i_res.get())
            ofq1_list.append(i_res.get()[0].filename)
            ofq2_list.append(i_res.get()[1].filename)
            unmapped_dic_list.append(i_res.get()[2])
    
        unmapped_union_dict = unmapped_dic_list[0]
        for di in unmapped_dic_list[1:]:
            for k, v in di.items():
                if k not in unmapped_union_dict:
                    unmapped_union_dict[k] = v
                else:
                    unmapped_union_dict[k] |= v
        for k, v in unmapped_union_dict.items():
            unmapped_union_dict[k] = len(v)

        ut.fileop.mergeFastqFile(fqlist=ofq1_list, inprefix="bcparsed_R1", outpath=f"{parsepath}", ftype="fastq")    # 加缩进
        ut.fileop.mergeFastqFile(fqlist=ofq2_list, inprefix="bcparsed_R2", outpath=f"{parsepath}", ftype="fastq")

    else:
        processSeqRecord(parse_subpath=parse_subpath, bc_pattern=bcpattern, fq1=fq1, fq2=fq2)


def run_parseFastq(workst:str, sampid:str, fq1:str, fq2:str, bcpat:str, corenum:int=-1, remain_tmp=False):

    out_fq_path = f'{workst}/{sampid}/01_parsed'          
    ut.sysop.mkdir_p(dir=out_fq_path)
        
    test_bc_obj = BcPattern(in_pattern=bcpat)

    start_time = datetime.datetime.now()
    parseBarcode(fq1path=fq1, fq2path=fq2, bcpattern=test_bc_obj, parsepath=f"{out_fq_path}", corenum=corenum)
    if not remain_tmp:
        shutil.rmtree(path=f"{out_fq_path}/tmp")
    end_time = datetime.datetime.now()
    time_eclip = (end_time - start_time).total_seconds()
    print(f"Take time : {time_eclip}")
    