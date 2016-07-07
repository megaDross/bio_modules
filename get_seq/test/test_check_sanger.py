from get_seq.check_sanger import *
from useful_tools import useful
import unittest
        
class TestStaticMethods(unittest.TestCase):
    
    test_dir = useful.cwd_file_path(__file__)+"test_files/"
    
    def test_get_matching_seq_file(self):
        self.assertEqual(CompareSeqs.get_matching_seq_file("test2", 
                         TestStaticMethods.test_dir), 
                         TestStaticMethods.test_dir+"test2_R.seq")

    def test_compare_nucleotides(self):
        self.assertEqual(CompareSeqs.compare_nucleotides("A","G"),
                         ("the nucleotides given are DIFFERENT",2))
        self.assertEqual(CompareSeqs.compare_nucleotides("C","C"),
                         ("the nucleotides given are the SAME",1))

    def test_get_start_end_indexes(self):
        self.assertEqual(CompareSeqs.get_start_end_indexes("TTCCTCCTTCAAACTTCGCA",   open(TestStaticMethods.test_dir+"F01_LX18_SP_UKB2_LX18_F_011.seq", "r").read().replace("\n", "")), (337, 357, "TTCCTCCTTCAAACTTCGCA"))


if __name__ == '__main__':
    unittest.main()


