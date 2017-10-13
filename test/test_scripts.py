from GeneaPy import get_seq, unknown_primer, primer_finder
import logging
import unittest
import os

HERE = os.path.dirname(os.path.realpath(__file__))
DATABASE = HERE+'/expected_output/primer_database.txt'

log = logging.getLogger()
log.disabled = True

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
                   "hg19", "FBN1", "FBN1-001", "-", "41/65", "421bp",
                   "chr15:48755298-48755718", 39.2)
        primer_info = unknown_primer.unknown_primer(f_primer='CTGTTCACAGGGCTTGTTCC',
                                                    r_primer='CTGGGCAGAGAGTCATTTAAAGT',
                                                    hg_version='hg19',
                                                    primer_name='query',
                                                    max_size=4000,
                                                    min_perfect=15,
                                                    min_good=15)
        self.assertEqual(primer_info, correct)

    def test_IO(self):
        args = {'input': HERE+'/expected_output/unknown_primer_in.txt',
                'genome_version': 'hg19',
                'max_size': 4000,
                'min_perfect': 15,
                'min_good': 15,
                'output': 'temp.txt'}
        header = '\t'.join(('Primer', 'F_Primer', 'R_Primer', 'Genome', 'Gene', 'Transcript', 'Exon',
                            'Intron', 'Product_Size', 'Primer_Range', 'GC%'))
        unknown_primer.parse2output(args, header)
        correct = open(HERE+'/expected_output/unknown_primer_out.txt').read()
        output = open('temp.txt').read()
        self.assertEqual(output, correct)

    def tearDown(self):
        try:
            os.remove('temp.txt')
        except FileNotFoundError:
            pass


class PrimerFinder(unittest.TestCase):
    def test_intron_filter(self):
        primer_finder.primer_finder(DATABASE, intron=21, output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_intron.txt').read()
        intron = open('temp.txt').read()
        self.assertEqual(intron, correct)

    def test_exon_filter(self):
        primer_finder.primer_finder(DATABASE, exon=9, output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_exon.txt').read()
        exon = open('temp.txt').read()
        self.assertEqual(exon, correct)

    def test_gene_filter(self):
        primer_finder.primer_finder(DATABASE, gene='FBN1', output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_gene.txt').read()
        fbn1 = open('temp.txt').read()
        self.assertEqual(fbn1, correct)

    def test_gc_filter(self):
        primer_finder.primer_finder(DATABASE, gc=60, output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_gc.txt').read()
        gc = open('temp.txt').read()
        self.assertEqual(gc, correct)

    def test_size_filter(self):
        primer_finder.primer_finder(DATABASE, size=350, output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_size.txt').read()
        size = open('temp.txt').read()
        self.assertEqual(size, correct)

    def test_variant_filter(self):
        primer_finder.primer_finder(DATABASE, variant='1:2160444', output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_variant.txt').read()
        variant = open('temp.txt').read()
        self.assertEqual(variant, correct)

    def test_distance_filter(self):
        primer_finder.primer_finder(DATABASE, variant='15:48787400', distance=50, output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_distance.txt').read()
        distance = open('temp.txt').read()
        self.assertEqual(distance, correct)

    def test_hg_filter(self):
        primer_finder.primer_finder(DATABASE, hg='hg19', output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_hg.txt').read()
        hg = open('temp.txt').read()
        self.assertEqual(hg, correct)

    def test_input_file(self):
        primer_finder.primer_finder(DATABASE, input_file=HERE+'/expected_output/primer_finder_input.txt', output='temp.txt')
        correct = open(HERE+'/expected_output/primer_finder_test_input.txt').read()
        input_file = open('temp.txt').read()
        self.assertEqual(input_file, correct)

    def tearDown(self):
        try:
            os.remove('temp.txt')
        except FileNotFoundError:
            pass


if __name__ == '__main__':
    unittest.main()
