from get_seq.Ensembl import *
import unittest
# python3 -m unittest test_Ensembl.py


class TestGeneTranscript(unittest.TestCase):
    def test_gene_info_hg19(self):
        self.assertEqual(ScrapeEnsembl("15:48762884", "hg19").get_gene_info(),(
            'FBN1', 'ENSG00000166147', 'protein_coding', '15:48700503-48938046'))
        self.assertEqual(ScrapeEnsembl('16:15812194', "hg19").get_gene_info(),(
            'MYH11', 'ENSG00000133392', 'protein_coding', '16:15797029-15950890'))
    def test_gene_info_hg38(self):        
        self.assertEqual(ScrapeEnsembl("15:48762884", "hg38").get_gene_info(),(
            'CEP152', 'ENSG00000103995', 'protein_coding', '15:48712928-48811146'))
        self.assertEqual(ScrapeEnsembl("16:15812194", "hg38").get_gene_info(),(
            'MYH11', 'ENSG00000133392', 'protein_coding', '16:15703135-15857033'))



if __name__ == '__main__':
    unittest.main()

        
