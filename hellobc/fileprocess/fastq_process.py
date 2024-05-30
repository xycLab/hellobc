import os
from enum import Enum
from Bio import SeqIO
import gzip


class fqType(Enum):
    FASTQ = 0
    FASTQ_GZ = 1


class fastqFile():

    def __init__(self, fname:str) -> None:
        self.filename = fname
        self.type, self.prefix, self.tail, self.fqstr = self._getType()
        
        
    def _getType(self) -> int:
        fpath_split = self.filename.split('/')
        fname_split = fpath_split[-1].split('.')
        #fprefix = fname_split[0]
        if fname_split[-1] == 'gz':
            if fname_split[-2] == 'fastq' or fname_split[-2] == 'fq':
                ftype = fqType.FASTQ_GZ
                fprefix = '.'.join(fname_split[0:-2])
                ftail = fname_split[-2] + '.' + fname_split[-1]
                ffqstr = fname_split[-2]
            else:
                raise ValueError("input file must be .fastq(fq) or .fastq(fq).gz")
        elif fname_split[-1] == 'fastq' or fname_split[-1] == 'fq':
            ftype = fqType.FASTQ
            fprefix = '.'.join(fname_split[0:-1])
            ftail = fname_split[-1]
            ffqstr = fname_split[-1]
        else:
            raise ValueError("input file must be .fastq(fq) or .fastq(fq).gz")
        return ftype, fprefix, ftail, ffqstr
    
    def gzipFastq(self) -> None:
        if self.type == fqType.FASTQ:
            os.system("gzip -c " + self.filename + " > " + self.filename + ".gz")
        else:
            print("The fastq file is already in gzip format!")

    def appendOneRecord(self, oprecord:SeqIO) -> None:
        if self.type == fqType.FASTQ:
            with open(self.filename, "a") as f:
                SeqIO.write(oprecord, f, "fastq")
            f.close()
        else:
            with open('.'.join((self.filename).split('.')[0:-1]), "a") as f:
                SeqIO.write(oprecord, f, "fastq")
            
            # with gzip.open(self.filename, "a") as f:
            #     SeqIO.write(oprecord, f, "fastq")
            f.close()
        

    def appendRecords(self, oprecords:list) -> None:
        n_rec = len(oprecords)
        if self.type == fqType.FASTQ:
            with open(self.filename, "a") as f:
                for i in range(n_rec):
                    SeqIO.write(oprecords[i], f, "fastq")
            f.close()
        else:
            with gzip.open(self.filename, "a") as f:
                for i in range(n_rec):
                    SeqIO.write(oprecords[i], f, "fastq")
            f.close()

