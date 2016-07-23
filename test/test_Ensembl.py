from GeneaPy.Ensembl import *
from GeneaPy.get_gene_exon_info import *
import unittest
# python3 -m unittest test_Ensembl.py


class TestEnsembl(unittest.TestCase):
    ''' Methods that test the ScarpeEnsembl 
        class found in Ensembl.py
    '''

    def test_gene_info_hg19(self):
        ''' Ensure the expected gene information is returned when
            parsing a genomic position (hg19)
        '''
        self.assertEqual(ScrapeEnsembl("15:48762884", "hg19").get_gene_info(),(
            'FBN1', 'ENSG00000166147', 'protein_coding', '15:48700503-48938046'))
        self.assertEqual(ScrapeEnsembl('16:15812194', "hg19").get_gene_info(),(
            'MYH11', 'ENSG00000133392', 'protein_coding', '16:15797029-15950890'))
    
    def test_gene_info_hg38(self):        
        ''' Ensure the expected gene information is returned when
            parsing a genomic position (hg38)
        '''
        self.assertEqual(ScrapeEnsembl("15:48762884", "hg38").get_gene_info(),(
            'CEP152', 'ENSG00000103995', 'protein_coding', '15:48712928-48811146'))
        self.assertEqual(ScrapeEnsembl("16:15812194", "hg38").get_gene_info(),(
            'MYH11', 'ENSG00000133392', 'protein_coding', '16:15703135-15857033'))

    def test_no_gene_found(self):
        ''' Test an error message is recieved when a genome position 
            associated with no gene is parsed
        '''
        self.assertEqual(ScrapeEnsembl("17:2415500", "hg19").get_gene_info(),
                         "No gene found at 17:2415500 for genome version 75")

    def test_canonical_transcript(self):
        ''' Ensure the correct canonical transcript is returned when 
            parsing a genomic position
        '''
        self.assertEqual(ScrapeEnsembl("15:48762884", "hg38").
                         get_canonical_transcript('CEP152'), 
                         "ENST00000380950")
        self.assertEqual(ScrapeEnsembl("21:33036150", "hg19").
                         get_canonical_transcript('SOD1'),
                         "ENST00000270142")

        
class TestExon(unittest.TestCase):
    ''' Methods that test the get_exon_number function 
        in get_gene_exon_info.py
    '''

    def test_get_exon_number_exon(self):
        ''' Test the expected exon information is returned when a
            transcript, human genome version and genomic position 
            is parsed
        '''
        self.assertEqual(get_exon_number("ENST00000316623", "hg19", "15:48762884"),
                         ('ENSE00003582511', '-', '36/66'))
        self.assertEqual(get_exon_number("ENST00000617185", "hg38", "17:7675200"),
                         ('ENSE00003518480', '-', '5/12'))

    def test_get_exon_number_intron(self):
        ''' As above but test for expected intron information is returned
        '''
        self.assertEqual(get_exon_number("ENST00000420246", "hg19", "17:7576660"),
                         ('-', '9/11', '-'))
        self.assertEqual(get_exon_number("ENST00000520751", "hg38", "8:127736724"),
                         ('-', '1/1', '-'))


class TestGeneTranscriptExon(unittest.TestCase):
    ''' Methods that test the gene_transcript_exon and create_position
        methods in get_exon_gene_info.py
    '''

    def test_function(self):
        '''
        '''
        self.assertEqual(gene_transcript_exon("15:48762884", "hg19"),
                         (('FBN1', 'ENSG00000166147', 'protein_coding', 
                           '15:48700503-48938046'), 'ENST00000316623',
                          ('ENSE00003582511', '-', '36/66')))

        self.assertEqual(gene_transcript_exon("16:15812194", "hg19"),
                         (('MYH11', 'ENSG00000133392', 'protein_coding', 
                           '16:15797029-15950890'), 'ENST00000452625',
                          ('ENSE00003528672', '-', '38/43')))
    
    def test_hg38(self):
        '''
        '''
        self.assertEqual(gene_transcript_exon("15:48762884", "hg38"),
                         (('CEP152', 'ENSG00000103995', 'protein_coding', 
                           '15:48712928-48811146'), 'ENST00000380950', 
                          ('-', '17/26', '-')))
        self.assertEqual(gene_transcript_exon("3:123418990", "hg38"),
                        (('ADCY5', 'ENSG00000173175', 'protein_coding', 
                          '3:123282296-123449758'), 'ENST00000462833', 
                         ('-', '1/20', '-')))

    #def test_exon_error(self):
    #    '''
    #    '''
    #    self.assertEqual(gene_transcript_exon("16:15812194", "hg38"),
    #                     "ERROR: No exon information found for 16:15812194 in hg38")


class TestCreatePosition(unittest.TestCase):
    

    def test_it(self):
        self.assertEqual(create_position("1:10000900-10001000"),
                         "1:10000950")
    
if __name__ == '__main__':
    unittest.main()


