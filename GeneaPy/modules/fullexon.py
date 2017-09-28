from pyensembl.exon import Exon
import modules.custom_exceptions as ex
from  modules import get_exon_intron_number as gein

class FullExon(Exon):
    ''' Enhances the Exon object by adding exon/intron numbers'''
    def __init__(self, exon_id, contig, start, end, strand, gene_name, gene_id, position, transcript_id, hg_version):
        Exon.__init__(self, exon_id, contig, start, end, strand, gene_name, gene_id)
        self.position = position
        self.transcript_id = transcript_id
        self.hg_version = hg_version
        
    @property
    def exon_no(self):
        exon_id, intron_num, exon_num = gein.main(
            transcript=self.transcript_id, 
            hg_version=self.hg_version, 
            pos=self.position)
        if self.id != exon_id:
            raise ex.ExonMismatch(self.id, exon_id)
        return exon_num
    
    @property
    def intron_no(self):
        exon_id, exon_num, intron_num = gein.main(
            transcript=self.transcript_id, 
            hg_version=self.hg_version, 
            pos=self.position)
        return intron_num
    
    @classmethod
    def from_pyexon(cls, exon, position, transcript_id, hg_version):
        ''' Parse a pyensembl Exon object instead'''
        return cls(
            exon.id, exon.contig, exon.start, exon.end, 
            exon.strand, exon.gene_name, exon.gene_id, 
            position, transcript_id, hg_version)
    
    def __str__(self):
        return "Exon(exon_id=%s, gene_name=%s, contig=%s, start=%d, end=%s, position=%s, exon_no=%s, intron_no=%s)" % (
            self.id, self.gene_name, self.contig, self.start, 
            self.end, self.position, self.exon_no, self.intron_no)
