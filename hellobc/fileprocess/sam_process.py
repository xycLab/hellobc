import pysam
from pysam.libcalignedsegment import AlignedSegment


class bamLine():
    
    def __init__(self, fname:str, read:AlignedSegment) -> None:
        self.filename = fname
        self.read = read
        self.XS = self._getXS()
        self.XT = self._getXT()
        self.CB, self.umi = self._splitSeqID()
        

    def _getXS(self) -> str:
        asign_tag = dict(self.read.tags)
        return asign_tag['XS']
        

    def _getXT(self) -> str:
        asign_tag = dict(self.read.tags)
        if 'XT' in asign_tag:
            return asign_tag['XT']
        else:
            return ''
        

    def _splitSeqID(self):
        splitid = self.read.query_name.split('_')
        readCB = splitid[-2]
        readumi = splitid[-1]
        return readCB, readumi


