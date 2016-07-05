# python3 -m pytest test_Ensembl.py
from get_seq.Ensembl import *


def check_gene_info():
    assert ScrapeEnsembl("15:48762884", "hg19").get_gene_info() == (
        'FBN1', 'ENSG00000166147', 'protein_coding', '15:48700503-48938046')
    assert ScrapeEnsembl('16:15812194', "hg19").get_gene_info() == (
        'MYH11', 'ENSG00000133392', 'protein_coding', '16:15797029-15950890')
    
