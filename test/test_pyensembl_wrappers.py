from GeneaPy.modules import pyensembl_wrappers
from GeneaPy.modules.fullexon import FullExon
from pyensembl import EnsemblRelease
import unittest

DATA = EnsemblRelease(75)

class TestWrappers(unittest.TestCase):
    def test_get_canonical_transcript(self):
        correct = DATA.transcript_by_id('ENST00000458208')
        canon = pyensembl_wrappers.get_canonical_transcript(DATA, 10, 90751147)
        self.assertEqual(canon, correct)

    def test_get_gene_locus(self):
        correct = DATA.gene_by_id('ENSG00000026103')
        gene = pyensembl_wrappers.get_gene_locus(DATA, 10, 90752100)
        self.assertEqual(gene, correct)

    def test_get_exon(self):
        ''' Test retrieving FullExon object.'''
        transcript = DATA.transcript_by_id('ENST00000316623')
        exon = pyensembl_wrappers.get_exon(48752450, transcript)
        correct = FullExon('ENSE00003605533', 15, 48752443, 48752514, '-',
                           'FBN1', 'ENSG00000166147', 48752450, '43/66', True)
        self.assertEqual(exon.__dict__, correct.__dict__)

    def test_get_intron_positive(self):
        transcript = DATA.transcript_by_id('ENST00000487314')
        intron = pyensembl_wrappers.get_exon(90752100, transcript)
        correct = FullExon('N/A', 10, 90751384, 90762785, '+', 'FAS', 
                           'ENSG00000026103', 90752100, '1/7', False) 
        self.assertEqual(intron.__dict__, correct.__dict__)
        
    def test_get_intron_negative(self):
        transcript = DATA.transcript_by_id('ENST00000316623')
        intron = pyensembl_wrappers.get_exon(48778271, transcript)
        correct = FullExon('N/A', 15, 48776141, 48777570, '-', 'FBN1',
                           'ENSG00000166147', 48778271, '29/65', False)
        self.assertEqual(intron.__dict__, correct.__dict__)
