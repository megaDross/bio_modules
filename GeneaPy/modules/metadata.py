from pyensembl import EnsemblRelease
from modules.fullexon import FullExon
from modules import pyensembl_wrappers
from modules import common
import get_seq

# TODO: unit testing
# TODO: exception handeling
# TODO: logging

# inherit from Sequence
class LocusMetaData(object):
    ''' Store the gene, transcript and exon metadata of a given genomic position.

    Parameters:
        contig: chromosome nuumber
        position: position number
        hg_version: human genome version
        flank: desired number of bp to scrape from each side of the given genomic position (default=50)
        genome: path to human genome FASTA file (optional)
    '''
    def __init__(self, contig, position, hg_version, flank=50, genome=None):
        self.contig = contig
        self.position = position
        self.hg_version = common.correct_hg_version(hg_version)
        self.flank = flank
        self.genome = genome
        self._transcript = None

    @property
    def ensembl(self):
        ensembl_hg = common.get_ensembl_release(self.hg_version)
        return EnsemblRelease(ensembl_hg)

    @property
    def gene(self):
        return pyensembl_wrappers.get_gene_locus(
            data=self.ensembl, 
            contig=self.contig, 
            position=self.position
        ) 
       
    @property
    def transcript(self):
        if not self._transcript:
            canonical = pyensembl_wrappers.get_canonical_transcript(
                data=self.ensembl, 
                contig=self.contig, 
                position=self.position
            )
            self._transcript = canonical
        return self._transcript
   
    @transcript.setter
    def transcript(self, transcript_id):
        new_transcript = self.ensembl.transcript_by_id(transcript_id)
        self._transcript = new_transcript

    @property
    def exon(self):
        pyexon = pyensembl_wrappers.get_exon_at_locus_of_transcript(
            data=self.ensembl,
            contig=self.contig,
            position=self.position,
            transcript=self._transcript
        )
        genomic_position = 'chr{}:{}'.format(self.contig, self.position)
        fullexon = FullExon.from_pyexon(pyexon, genomic_position, 
                                        self._transcript.id, self.hg_version)
        return fullexon

    @property
    def sequence(self):
        query = '{}:{}'.format(self.contig, self.position)
        seq = get_seq.main(query, self.hg_version, self.genome, self.flank, self.flank, header=False)
        return seq

    @property
    def seq_range(self):
        start = self.position - self.flank
        end = self.position + self.flank
        return '{}:{}-{}'.format(self.contig, start, end)

    @classmethod
    def from_position(cls, genomic_position, hg_version, flank=50, genome=None):
        contig, position = genomic_position.lower().replace('chr', '').split(':')
        return cls(contig=int(contig), position=int(position),
                   hg_version=hg_version, flank=flank, genome=genome)

    def __str__(self):
        s = '-'*50
        query = '{}\nchr{}:{} in {}\n{}\n\n'.format(
                 s, self.contig, self.position, self.hg_version,  s)
        gene_location = '{}:{}-{}'.format(self.gene.contig, self.gene.start, self.gene.end)
        gene = 'Gene\nname: {}\nid: {}\nloc: {}\ntype: {}\n\n'.format(
                self.gene.name, self.gene.id, gene_location, self.gene.biotype)
        transcript = 'Transcript\nid: {}\n\n'.format(
                      self.transcript.id)
        exon = 'Exon\nid: {}\nexon: {}\nintron: {}\n\n'.format(
                self.exon.id, self.exon.exon_no, self.exon.intron_no)
        scrapped_seq = '{}\n\n>{}\n{}\n\n{}'.format(
                        s, self.seq_range, self.sequence, s)
        all_metadata = (query+gene+transcript+exon+scrapped_seq)        
        return all_metadata
