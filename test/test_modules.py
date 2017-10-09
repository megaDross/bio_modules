from GeneaPy.modules import pyensembl_wrappers
from GeneaPy.modules.fullexon import FullExon
from GeneaPy.modules.metadata import LocusMetaData
from GeneaPy.modules import common
from pyensembl import EnsemblRelease, Exon, Gene, Transcript
import unittest

DATA = EnsemblRelease(75)

class TestFullExon(unittest.TestCase):
    def test_full_exon(self):
        ''' Tests the objects attributes are the same'''
        exon = Exon('ENSE00003605533', 15, 48752443, 48752514, '-',
                    'FBN1', 'ENSG00000166147')
        exon.position = 48752450
        exon.number = '43/66'
        exon.exon = True
        full_exon = FullExon('ENSE00003605533', 15, 48752443, 48752514, '-',
                             'FBN1', 'ENSG00000166147', 48752450, '43/66', True)
        self.assertEqual(full_exon.__dict__, exon.__dict__)


class TestPyensemblWrappers(unittest.TestCase):
    def test_get_canonical_transcript(self):
        correct = DATA.transcript_by_id('ENST00000458208')
        canon = pyensembl_wrappers.get_canonical_transcript(DATA, 10, 90751147)
        self.assertEqual(canon, correct)

    def test_get_gene_locus(self):
        correct = DATA.gene_by_id('ENSG00000026103')
        gene = pyensembl_wrappers.get_gene_locus(DATA, 10, 90752100)
        self.assertEqual(gene, correct)

    def test_get_exon(self):
        transcript = DATA.transcript_by_id('ENST00000316623')
        exon = pyensembl_wrappers.get_exon(48752450, transcript)
        correct = FullExon('ENSE00003605533', 15, 48752443, 48752514, '-',
                           'FBN1', 'ENSG00000166147', 48752450, '43/66', True)
        self.assertEqual(exon.__dict__, correct.__dict__)

    def test_get_intron_positive(self):
        ''' Expects a FullExon object of an Intron on positive strand '''
        transcript = DATA.transcript_by_id('ENST00000487314')
        intron = pyensembl_wrappers.get_exon(90752100, transcript)
        correct = FullExon('N/A', 10, 90751384, 90762785, '+', 'FAS', 
                           'ENSG00000026103', 90752100, '1/7', False) 
        self.assertEqual(intron.__dict__, correct.__dict__)
        
    def test_get_intron_negative(self):
        ''' Expects a FullExon object of an Intron on negative strand '''
        transcript = DATA.transcript_by_id('ENST00000316623')
        intron = pyensembl_wrappers.get_exon(48778271, transcript)
        correct = FullExon('N/A', 15, 48776141, 48777570, '-', 'FBN1',
                           'ENSG00000166147', 48778271, '29/65', False)
        self.assertEqual(intron.__dict__, correct.__dict__)


class TestCommon(unittest.TestCase):
    def test_correct_hg19_version(self):
        hg_version = common.correct_hg_version('GrCh37')
        self.assertEqual(hg_version, 'hg19')
    
    def test_correct_hg38_version(self):
        hg_version = common.correct_hg_version('GrCh38')
        self.assertEqual(hg_version, 'hg38')
   
    def test_get_ensembl_release_hg19(self):
        ensembl_release = common.get_ensembl_release('hg19')
        self.assertEqual(ensembl_release, 75)

    def test_get_ensembl_release_hg38(self):
        ensembl_release = common.get_ensembl_release('hg38')
        self.assertEqual(ensembl_release, 83)


class TestMetaData(unittest.TestCase):
    correct = {'genome': None, 
               'ensembl': DATA, 
               'flank': 50, 
               'sequence': 'atagtgaatgggacagacacaatctttgacttc'
                           'aaaatgattacactgtg\nGccaggagacagat'
                           'gaacaattaattgcaccatgcatgatgtgccat'
                           'ttg\nc', 
               '_transcript': None, 
               'gene': Gene(gene_id='ENSG00000166147', gene_name='FBN1', biotype='protein_coding', 
                            contig='15', start=48700503, end=48938046, strand='-', genome=DATA), 
               'position': 48778271, 
               'contig': 15, 
               'gene_list': [],
               'hg_version': 'hg19'}
    metadata = LocusMetaData(15, 48778271, 'hg19') 

    def test_metadata(self):
        self.assertEqual(self.metadata.__dict__, self.correct)

    def test_metadata_exon(self):
        intron = FullExon(exon_id='N/A', gene_name='FBN1', contig=15, start=48776141, end=48777570, 
                          position=48778271, number='29/65', strand='-', gene_id='ENSG00000166147',
                          exon=False)
        self.assertEqual(self.metadata.exon, intron)

    def test_metadata_transcript(self):
        transcript = Transcript(transcript_id='ENST00000316623', transcript_name='FBN1-001', strand='-', 
                                biotype='protein_coding', contig=15, start=48700503, end=48938046,
                                genome=DATA, gene_id='ENSG00000166147')
        self.assertEqual(self.metadata.transcript, transcript)


class TestMetaData2(unittest.TestCase):
    ''' Tests metadata object if position overlaps two genes and the 
        automatically select gene is as expected.
    '''
    correct = {'genome': None, 
               'ensembl': DATA, 
               'flank': 50, 
               'sequence': 'actaaaagtagttcctggttggtgaaaataaatcattaatgcgttttaaa\n'
                           'Tgaaaaagaaatgcatgcgtcttgtaaaaaatgtgaaataaaagaggcat\na',
               '_transcript': None, 
               'gene': Gene(gene_id='ENSG00000267699', gene_name='RP11-729L2.2',
                            biotype='protein_coding', contig='18', start=48494389,
                            end=48584514, strand='+', genome=DATA),
               'position': 48555816, 
               'contig': 18, 
               'gene_list': [],
               'hg_version': 'hg19'}
    metadata = LocusMetaData(18, 48555816, 'hg19') 

    def test_metadata(self):
        self.assertEqual(self.metadata.__dict__, self.correct)

    def test_metadata_exon(self):
        intron = FullExon(exon_id='N/A', gene_name='RP11-729L2.2', 
                          contig=18, start=48500932, end=48573289, 
                          position=48555816, number='2/8', strand='+', 
                          gene_id='ENSG00000267699',
                          exon=False)
        self.assertEqual(self.metadata.exon, intron)

    def test_metadata_transcript(self):
        transcript = Transcript(transcript_id='ENST00000590722', 
                                transcript_name='RP11-729L2.2-001', strand='+', 
                                biotype='nonsense_mediated_decay', contig=18, start=48494389, end=48584514,
                                genome=DATA, gene_id='ENSG00000267699')
        self.assertEqual(self.metadata.transcript, transcript)


class TestMetaData3(unittest.TestCase):
    ''' Tests metadata object if position overlaps two genes and one uses
        gene_list to select one explicitly
    '''
    correct = {'genome': None, 
               'ensembl': DATA, 
               'flank': 50, 
               'sequence': 'actaaaagtagttcctggttggtgaaaataaatcattaatgcgttttaaa\n'
                           'Tgaaaaagaaatgcatgcgtcttgtaaaaaatgtgaaataaaagaggcat\na',
               '_transcript': None, 
               'gene': Gene(gene_id='ENSG00000141646', gene_name='SMAD4',
                            biotype='protein_coding', contig='18', start=48494410,
                            end=48611415, strand='+', genome=DATA),
               'position': 48555816, 
               'contig': 18, 
               'gene_list': ['SMAD4'],
               'hg_version': 'hg19'}
    metadata = LocusMetaData(18, 48555816, 'hg19', gene_list=['SMAD4']) 

    def test_metadata(self):
        self.assertEqual(self.metadata.__dict__, self.correct)

    def test_metadata_exon(self):
        intron = FullExon(exon_id='N/A', gene_name='SMAD4', 
                          contig=18, start=48500932, end=48573289, 
                          position=48555816, number='2/8', strand='+', 
                          gene_id='ENSG00000141646',
                          exon=False)
        self.assertEqual(self.metadata.exon, intron)

    def test_metadata_transcript(self):
        transcript = Transcript(transcript_id='ENST00000452201', 
                                transcript_name='SMAD4-202', strand='+', 
                                biotype='protein_coding', contig=18, start=48494410,
                                end=48584514,
                                genome=DATA, gene_id='ENSG00000141646')
        self.assertEqual(self.metadata.transcript, transcript)



if __name__ == '__main__':
    unittest.main()
