import subprocess, warnings
import unittest
from useful_tools import useful


# use and manipulate file_apth to get path to get_seq.py
file_path = useful.cwd_file_path(__file__)
get_seq = file_path[:-5]+"get_seq.py"



class TestGetSeqPrint(unittest.TestCase):
    ''' Test the printed output of the get_seq.py click application
    '''
    def test_genomic_position(self):
        ''' Test the printed output from parsing a genomic position is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\nReference Sequence:\taggagcgatgtctcctct'
                                 'ttCattcctgattttgataattt\n\n\n')
    
    def test_piping(self):
        ''' Check that piping a genomic position into get_seq.py works
        '''
        # stop ResourceWarning printing to screen
        warnings.filterwarnings("ignore")

        # use subprocess to test piping a position into the get_seq.py click application,
        # grab the stdout from the resulting tuple and then convert to string
        echo = subprocess.Popen(["echo", "chr1:169314424"], stdout=subprocess.PIPE)
        ref_seq_popen = subprocess.Popen(["python3", get_seq], stdin=echo.stdout,
                                                               stdout=subprocess.PIPE)
        ref_seq_bytes = ref_seq_popen.communicate()[0]
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\nReference Sequence:'
                         '\taggagcgatgtctcctctttCattcctgattttgataattt\n\n\n')
 
    def test_genomic_range(self):
        ''' Test the printed output from parsing a genomic range is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, 
                                                 "chr1:169314404,169314444"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\nReference Sequence:\taggagcgatgtctcctct'
                                 'ttCattcctgattttgataattt\n\n\n')
 
    def test_hg_version(self):
        ''' Test the printed output from parsing a genomic position and hg38
            is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr15:48762884",
                                                "--hg_version", "hg38"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\nReference Sequence:\tgaacttaactatatgacaaaA'
                         'atcacatgaaagatttaagt\n\n\n')
    
   
    def test_header(self):
        ''' Test the printed output from parsing a genomic position and header
            option is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424",
                                                "--header"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,
                         '> query chr1:169314424 chr1:169314404,169314444'
                         '\nReference Sequence:\taggagcgatgtctcctct'
                                 'ttCattcctgattttgataattt\n\n\n')

    def test_upstream_downstream(self):
        ''' Test the printed output from parsing a genomic position, header option
            , upstream and downstream options is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424",
                                                "--header", "--upstream", "40",
                                                "--downstream", "70"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,
                         '> query chr1:169314424 chr1:169314384,169314494'
                         '\nReference Sequence:\tcttttaatttctgcagggatag'
                         'gagcgatgtctcctctttCattcctgattttgataatttgtatcttcagtcttttttct'
                         'tggtcagtctaaccaaagttttgcctatgt\n\n\n')


    def test_seq_file(self):
        ''' Test the printed output from parsing a genomic position, header option
            , upstream and downstream options is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr15:48762884",
                                                "--header", "--seq_file",
                                                 file_path[:-5]+'test/test_files/'
                                                'B02_CX1_AD_UKB2_CX1_F_004.seq'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,
                         '> query chr15:48762884 chr15:48762864,48762904'
                         '\nReference Sequence:\t'
                         'agcctatctcacactcacagCggaacaggccagggaggttg'                    
                         '\nSanger Sequence:\t'
                         'agcctatctcacactcacagG/Cggaacaggccaggg'
                         'aggttg\nthe nucleotides given are '
                         'DIFFERENT\n\n\n'                       
                        )

    def test_input_file(self):
        ''' Test the printed output from parsing an input file is as expected and that
            the correct seq file is automatically selected for comparing with the 
            reference sequence
        
            WARNING: the below assertEqual function is painfully ugly
        '''
        ref_seq_bytes = subprocess.check_output(["python3", get_seq, file_path[:-5]+
                                                 'test/test_in.txt', "--seq_dir",
                                                 file_path[:-5]+'test/test_files/',
                                                 '--header'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')

        self.assertEqual(ref_seq,'> CX1_AD 15:48762884 15:48762864,48762904'
                         '\nReference Sequence:\tagcctatctcacactcacagCggaacaggccagg'
                         'gaggttg\nSanger Sequence:\tagcctatctcacactcacag'
                         'G/Cggaacaggccagggaggttg\nthe nucleotides given are DIFFERENT'
                         
                         '\n\n\n> 24AS15 15:48776065 15:48776045,48776085\n'
                         'Reference Sequence:\ttgaatccatcataacacaagCacctgtactctccaggga'
                         'ta\nSanger Sequence:\ttgaatccatcataacacaagCacctgtactctccaggg'
                         'ata\nthe nucleotides given are the SAME\n\n\n'

                         '> 24AS15 15:48762884 15:48762864,48762904\n'
                         'Reference Sequence:\tagcctatctcacactcacagCggaacaggccagggagg'
                         'ttg\n\n\n'

                         '> 24AS15 15:48776065 15:48776045,48776085\n'
                         'Reference Sequence:\ttgaatccatcataacacaagCacctgtactctccaggga'
                         'ta\nSanger Sequence:\ttgaatccatcataacacaagCacctgtactctccaggg'
                         'ata\nthe nucleotides given are the SAME\n\n\n'

                         '> HX15_RD_GS 16:15812194 16:15812174,15812214\n'
                         'Reference Sequence:\tgctgtgtggctttgcggaccCggtcgctcatggcctcc'
                         'atg\nSanger Sequence:\tgctgtgtggctttgcggaccNggtcgctcatggcctc'
                         'catg\nthe nucleotides given are DIFFERENT\n\n\n'

                         '> LX15_NY 15:48707763 15:48707743,48707783\n'
                         'Reference Sequence:\ttgcggaagtaaccaggtggaCagccacacaggtaacc'
                         'gccc\nSanger Sequence:\ttgcggaagtaaccaggtggaNagccacacaggtaacc'
                         'gccc\nthe nucleotides given are DIFFERENT\n\n\n'

                         '> LX16_SB 15:48712915 15:48712895,48712935\n'
                         'Reference Sequence:\tttccactggtagtgctggagGtagccctgggggcag'
                         'ctgca\nSanger Sequence:\tttccactggtagtgctggagG/Ttagcc'
                         'ctgggggcagctgca\nthe nucleotides given are DIFFERENT\n\n\n'

                         '> FUK27_MM 15:48892343 15:48892323,48892363\n'
                         'Reference Sequence:\tacagtgtacttacgttgtccAcagtgagtccc'
                         'tatgtatcc\nSanger Sequence:\tacagtgtacttacgttgtccAcagtg'
                         'agtccctatgtatcc\nthe nucleotides given are the SAME\n\n\n'

                         '> FUK26_SH 15:48812913 15:48812893,48812933\n'
                         'Reference Sequence:\tgacccctggagaccagcatcGgccggcatcacagcagc'
                         'act\nSanger Sequence:\tgacccctggagaccagcatcNgccggcatcacagc'
                         'agcact\nthe nucleotides given are DIFFERENT\n\n\n'

                         '> LX18_SP 15:48730093 15:48730073,48730113\n'
                         'Reference Sequence:\tttcctccttcaaacttcgcaTaacagtagctcattcg'
                         'caaa\nSanger Sequence:\tttcctccttcaaacttcgcaTaammrnnacnnmwt'  
                         'nccnaa\nthe nucleotides given are the SAME\n\n\n'

                         '> test_R 15:48812913 15:48812893,48812933\n'
                         'Reference Sequence:\tgacccctggagaccagcatcGgccggcatcacagcag'
                         'cact\nSanger Sequence:\tgacccctggagaccagttttNgccggcatcacag'
                         'cagcact\nthe nucleotides given are DIFFERENT\n\n\n'

                        '> test2_R 15:48892343 15:48892323,48892363\n'
                        'Reference Sequence:\tacagtgtacttacgttgtccAcagtgagtccctatg'
                        'tatcc\nSanger Sequence:\tacagtgtacttatgttggggAcagtgagtccctat'
                        'gtatcc\nthe nucleotides given are the SAME\n\n\n')






if __name__ == '__main__':
    unittest.main()
