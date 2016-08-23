from GeneaPy.check_sanger import *
from GeneaPy import useful
import unittest, warnings


test_dir = useful.cwd_file_path(__file__)+"test_files/"

class TestCheckSanger(unittest.TestCase):
    
    
    def test_get_matching_seq_file(self):
        ''' Ensure the query matches the expected file in he selected dir
        '''
        self.assertEqual(CompareSeqs.get_matching_seq_file("24HC1541-WYN1", 
                         test_dir)[0], 
                         test_dir + "F07_24HC1541_WYN1_011.ab1")

    def test_compare_nucleotides(self):
        ''' Test both conditional outputs of the compare_nucleotides method
        '''
        self.assertEqual(CompareSeqs.compare_nucleotides("A","G"),
                         ("the nucleotides given are DIFFERENT",2))
        self.assertEqual(CompareSeqs.compare_nucleotides("C","C"),
                         ("the nucleotides given are the SAME",1))

    def test_get_start_end_indexes_postseq(self):
        ''' Test that the preseq returns the expected index numbers 
        '''
        warnings.filterwarnings("ignore")
        self.assertEqual(CompareSeqs.get_start_end_indexes("TTCCTCCTTCAAACTTCGCA",   CompareSeqs.handle_seq_file(test_dir+"F01_LX18_SP_UKB2_LX18_F_011.ab1", test_dir)), (337, 357, "TTCCTCCTTCAAACTTCGCA"))

    def test_get_start_and_indexes_preseq(self):
        ''' Test that the postseq returns the expected index numbers
        '''
        warnings.filterwarnings("ignore")
        self.assertEqual(CompareSeqs.get_start_end_indexes("agcctatctcacactcacag".upper()
, open(test_dir+"B02_CX1_AD_UKB2_CX1_F_004.seq", "r").read().replace("\n", "")), (270, 290, "AGCCTATCTCACACTCACAG"))

    def test_match_with_seq_file_preseq(self):
        ''' Test that the presequence is found within the seq file
        '''
        warnings.filterwarnings("ignore")
        self.assertEqual(CompareSeqs(20, 20, test_dir+"B02_CX1_AD_UKB2_CX1_F_004.seq").match_with_seq_file("agcctatctcacactcacagCggaacaggccagggaggttg"),('agcctatctcacactcacag', 'ggaacaggccagggaggttg', 'C', 'G/C', 290))

    def test_match_with_seq_file_postseq(self):
        ''' Test that the postsequence is found within the seq file
        '''
        warnings.filterwarnings("ignore")
        self.assertEqual(CompareSeqs(60, 60, test_dir+"C02_HX15_RD_GS_U_HX15_RD_F_006.seq").match_with_seq_file("agtaggcagcgtgactgtggtgtccaggcggccctcacctgctgtgtggctttgcggaccCggtcgctcatggcctccatgttgccctgctcctcctccagctcctcctccagctgggcga"),('agtaggcagcgtgnctgtggtgtccaggcggccctcacctgctgtgtggctttgcggacc', 'ggtcgctcatggcctccatgttgccctgctcctcctccagctcctcctccagctgggcga', 'C', 'N', 119))


    def test_ab1_file_name_conversion(self):
        self.assertEqual(CompareSeqs.convert_ab1_to_seq(test_dir+
                         "A09_29XX1917_WYN13_F_001.ab1"), test_dir+
                         "A09_29XX1917_WYN13_F_001.seq")
    
    def test_converted_ab1_sequence(self):
        ''' Tests SeqIO in BioPython is converting AB1 files correctly
        '''
        warnings.filterwarnings("ignore")
        self.assertEqual(open(test_dir+"A09_29XX1917_WYN13_F_001.seq").
                         read().replace("\n", ""), "NNNNNNNNNNNCNNNCNNNNGANGGTGAGACCTGGTCCCTGGGCCACTTGCNNAGCCCCTTCCACGCTGCCCTCACCTTAGCACCATYGTTGCCGGGAGCACCGTTGGCCCCTCGGGGACCAGCAGGACCAGGGGGACCTTGCACACCACGCTCGCCAGGGAAACCTCTCTCGCCCTAGAAGGGAAGGACAGGGCATGTGAAGGCTGCTCTGGAGATAGGGCCAAGTACAACGCACCTTGACGGATGCAGCGAGAGAGGCCTACTTACTCTTGCTCCAGAGGGGCCAGGGGCGCCAAGGTCTCCAGGAACACCCTGAGGGGGAGGGAGAGAGGAACAGACAGTGAGCAAAACCTACCTGGGGCCACTCTGCTGGAAAGCACGGTCCTTCCTCCTGGGGTCTGGCCGTGATTAGAGAGGAACCCCTTCTCAGCACTGAATTGAGATTATCCCAAACAGCCCCTCTCTTCCTCCTAGGGATGTGTCAAAGGCTTCCCCCCTCCCCAAACTATGGACCAAGATTTATCAATAGAAGGGTTGAGGGAAAGTCACAGNNNNNA")

       
class TtunerHetCalls(unittest.TestCase):

    def test_handle_seq_file(self):
        ''' ensure .seq, .abi files are processed s expected
        '''
        self.maxDiff = None
        ab1 = CompareSeqs.handle_seq_file(test_dir+"F01_LX18_SP_UKB2_LX18_F_011.ab1", 
                                          test_dir)
        seq = CompareSeqs.handle_seq_file(test_dir+"F01_LX18_SP_UKB2_LX18_F_011.ab1.seq", test_dir) 
        self.assertEqual(ab1, 'AAAAAAACTTAGCAGTTGTGTCTACTCCTTGGGCATCTGTGATTCATGACATAATTAAAAAGAAAATAACTACCAAGCTCCTGAGATAAAATAATTTTTAATGAACATTCTAATTAGACTTACAATTAAATTAAATTGTGTTACTGTCTTTAAGGCCTACAGTCTTACTTACATCATGGCCAGTCTGCACCCTGCATGGCCCAGAGAGAAATGCAGATGACAGACATACCATCAGGTTCCGTGGGGCAGAGCTCGCAGGGGTCTCCCCAGCCTTCTCCCTTCAAGGCACAGCAGCATTCCTGCTTGGAGTGATTTCTGGATTTGGGTGATGAACACTTTCCTCCTTCAAACTTCGCATAACCATAACCCCATCCCCAATCCGCCACAAAAATTTATGACCCCCCTTCATTGGCTTCCCAACGGATTCCTAATTGCCCAAAAAA')
        self.assertEqual(ab1,seq)

    

if __name__ == '__main__':
    unittest.main()


