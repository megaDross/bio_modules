from GeneaPy import get_seq, unknown_primer
import unittest
import os

class GetSeq(unittest.TestCase):
    def test_get_seq(self):
        correct = 'agacacttacCttggcacctt'
        seq = get_seq.get_seq('chr15:48733918', upstream=10, 
                              downstream=10, header=False,
                              hg_version='hg19')
        self.assertEqual(seq, correct)
    
    def test_header(self):
        correct = '>15:48733908,48733928 hg19\nagacacttacCttggcacctt'
        seq = get_seq.get_seq('chr15:48733918', upstream=10,
                              downstream=10, header=True,
                              hg_version='hg19')
        self.assertEqual(seq, correct)

    def test_hg38(self):
        correct = 'actgtccttcAtatataaatt'
        seq = get_seq.get_seq('chr15:48733918', upstream=10,
                              downstream=10, header=False,
                              hg_version='hg38')
        self.assertEqual(seq, correct)

    def test_range(self):
        correct = 'cttggcaccttcttccactggaggacaaggaaaacccttctggacacaga\n' \
                  'catttgaagctgccttcagtgttactgcatgt'
        seq = get_seq.get_seq('chr15:48733918-48733999', header=False)
        self.assertEqual(seq, correct)


class UnknownPrimer(unittest.TestCase):
    def test_unknown_primer(self):
        correct = ("query","CTGTTCACAGGGCTTGTTCC","CTGGGCAGAGAGTCATTTAAAGT",
                   "hg19", "FBN1", "-", "41/65", "421bp","chr15:48755298-48755718", 39.2)
        primer_info = unknown_primer.unknown_primer(f_primer='CTGTTCACAGGGCTTGTTCC',
                                                    r_primer='CTGGGCAGAGAGTCATTTAAAGT',
                                                    hg_version='hg19',
                                                    primer_name='query',
                                                    max_size=4000,
                                                    min_perfect=15,
                                                    min_good=15)
        self.assertEqual(primer_info, correct)

    def check_new_feature(self):
        # this primer pair lies within 2 genes. Use a new option for preferable genes
        primer_info = unknown_primer.unknown_primer(f_primer='CCAGGTCCAGATTCAGAGCC', 
                                                    r_primer='GACATGGCGCGGTTACCT',
                                                    hg_version='hg19',
                                                    primer_name='query',
                                                    max_size=4000,
                                                    min_perfect=15,
                                                    min_good=15)


if __name__ == '__main__':
    unittest.main()
