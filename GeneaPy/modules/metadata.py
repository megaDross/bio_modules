from pyensembl import EnsemblRelease
from GeneaPy.modules import pyensembl_wrappers
from GeneaPy.modules.common import correct_hg_version, get_ensembl_release
import GeneaPy.modules.custom_exceptions as ex
from GeneaPy import get_seq

# TODO: unit testing
# TODO: exception handeling
# TODO: logging

class LocusMetaData(object):
    ''' Store the gene, transcript and exon metadata of a given genomic position.

    Parameters:
        contig: chromosome nuumber
        position: position number
        hg_version: human genome version
        flank: desired number of bp to scrape from each side of the given genomic position (default=50)
        genome: path to human genome FASTA file (optional)
        gene_list: preffered gene(s) if position covers more than one gene
        seq: decide whether to scrape sequence from UCSC or not (Boolean)
    '''
    def __init__(self, contig, position, hg_version, flank=50, genome=None, gene_list=[], seq=True):
        self.contig = contig
        self.position = position
        self.hg_version = correct_hg_version(hg_version)
        self.flank = flank
        self.genome = genome
        self.gene_list = gene_list
        self.seq = seq
        self.ensembl = EnsemblRelease(get_ensembl_release(self.hg_version))
        self.gene = self._get_gene()
        self._transcript = None
        self.sequence = self._sequence()

    def _get_gene(self):
        return pyensembl_wrappers.get_gene_locus(
            data=self.ensembl, 
            contig=self.contig, 
            position=self.position,
            gene_list=self.gene_list
        ) 
 
    @property
    def transcript(self):
        transcript = pyensembl_wrappers.get_transcript(
            data=self.ensembl, 
            contig=self.contig, 
            position=self.position, 
            gene_list=self.gene_list
        )
        self._transcript = transcript
        return transcript

    @transcript.setter
    def transcript(self, transcript_id):
        new_transcript = self.ensembl.transcript_by_id(transcript_id)
        self._transcript = new_transcript

    @property
    def exon(self):
        if not self._transcript:
            self.transcript
        return pyensembl_wrappers.get_exon(self.position, self._transcript)

    def _sequence(self):
        if not self.seq:
            return None
        query = '{}:{}'.format(self.contig, self.position)
        seq = get_seq.get_seq(query, self.hg_version, self.genome, self.flank, self.flank, header=False)
        return seq

    def _get_seq_range(self):
        start = self.position - self.flank
        end = self.position + self.flank
        return '{}:{}-{}'.format(self.contig, start, end)

    @classmethod
    def from_position(cls, genomic_position, hg_version, flank=50, genome=None, gene_list=[]):
        contig, position = genomic_position.lower().replace('chr', '').split(':')
        return cls(contig=int(contig), position=int(position),
                   hg_version=hg_version, flank=flank, genome=genome, gene_list=gene_list)

    def __str__(self):
        s = '-'*50
        query = '{}\nchr{}:{} in {}({})\n{}\n\n'.format(
                s, self.contig, self.position, self.hg_version, self.ensembl.release, s)
        gene_location = '{}:{}-{}'.format(self.gene.contig, self.gene.start, self.gene.end)
        gene = 'Gene\nname: {}\nid: {}\nloc: {}\ntype: {}\n\n'.format(
               self.gene.name, self.gene.id, gene_location, self.gene.biotype)
        transcript = 'Transcript\nname: {}\nid: {}\ntype: {}\ncanon: {}\n\n'.format(
                     self.transcript.name, self.transcript.id, self.transcript.biotype, 
                     self.transcript.canonical)
        name = 'exon' if self.exon.exon else 'intron'
        exon = '{}\nid: {}\nno: {}\nstrand: {}\n\n'.format(
               name.capitalize(), self.exon.id, self.exon.number, self.exon.strand)
        scrapped_seq = '{}\n>{}\n{}\n{}'.format(
                       s, self._get_seq_range(), self.sequence, s)
        all_metadata = (query+gene+transcript+exon+scrapped_seq)        
        return all_metadata
