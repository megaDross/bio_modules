import subprocess, warnings, os, unittest
import GeneaPy.useful


# use and manipulate file_apth to get path to get_seq.py
file_path = GeneaPy.useful.cwd_file_path(__file__)
get_seq = file_path[:-5]+"GeneaPy/get_seq.py"


class TestGetSeqPrint(unittest.TestCase):
    ''' Test the printed output of the get_seq.py click application
    '''
    def test_genomic_position(self):
        ''' Test the printed output from parsing a genomic position is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\naggagcgatgtctcctctttCattcctgattttgataattt\n\n\n')
    
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
        
        self.assertEqual(ref_seq,'\naggagcgatgtctcctctttCattcctgattttgataattt\n\n\n')
 
    def test_genomic_range(self):
        ''' Test the printed output from parsing a genomic range is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, 
                                                 "chr1:169314404,169314444"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\naggagcgatgtctcctctttCattcctgattttgataattt\n\n\n')
 
    def test_hg_version(self):
        ''' Test the printed output from parsing a genomic position and hg38
            is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr15:48762884",
                                                "--hg_version", "hg38"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\ngaacttaactatatgacaaaAatcacatgaaagatttaagt\n\n\n')
    
   
    def test_header(self):
        ''' Test the printed output from parsing a genomic position and header
            option is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424",
                                                "--header"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,
                         '> query chr1:169314424 chr1:169314404,169314444 -'
                         '\naggagcgatgtctcctctttCattcctgattttgataattt\n\n\n')

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
                         '> query chr1:169314424 chr1:169314384,169314494 -'
                         '\ncttttaatttctgcagggatag'
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
                         '> query chr15:48762884 chr15:48762864,48762904 -'
                         '\nReference Sequence:\t'
                         'agcctatctcacactcacagCggaacaggccagggaggttg'                    
                         '\nSanger Sequence:\t'
                         'agcctatctcacactcacagG/Cggaacaggccaggg'
                         'aggttg\nthe nucleotides given are '
                         'DIFFERENT\n\n\n')


    def test_input_file(self):
        ''' Test the printed output from parsing an input file is as expected and that
            the correct seq file is automatically selected for comparing with the 
            reference sequence
        '''
        # ignore Resource Mangemenet error
        warnings.filterwarnings("ignore")
        ref_seq_bytes = subprocess.check_output(["python3", get_seq, file_path[:-5]+
                                                 'test/test_in.txt', "--seq_dir",
                                                 file_path[:-5]+'test/test_files/',
                                                 '--header', '--ensembl'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')

        self.assertEqual(ref_seq, open(file_path[:-5]+'test/test_out_print.txt').read())


    def test_ab1_conversion(self):
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "17:48269160",
                                                "--header", "--seq_file",
                                                 file_path[:-5]+'test/test_files/'
                                                'A09_29XX1917_WYN13_F_001.ab1'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,
                         '> query 17:48269160 17:48269140,48269180 -'
                         '\nReference Sequence:\t'
                         'tgccctcaccttagcaccatCgttgccgggagcaccgttgg'                    
                         '\nSanger Sequence:\t'
                         'tgccctcaccttagcaccatC/Tgttgccgggagcaccgttgg'
                         '\nthe nucleotides given are DIFFERENT\n\n\n')



class GetSeqFileOut(unittest.TestCase):
    ''' Test the written output file of the get_seq.py click application
    '''
    def test_output_file(self):
        '''  Test the written ouput from parsing an input file is as expected.
        '''
        # ignore Resource Management error
        warnings.filterwarnings("ignore")
        ref_seq_bytes = subprocess.check_output(["python3", get_seq, file_path[:-5]+
                                                 'test/test_in.txt', "--seq_dir",
                                                 file_path[:-5]+'test/test_files/',
                                                 '--header', '--output_file',
                                                  file_path[:-5]+'test/yyy.txt',
                                                 '--ensembl'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')

        self.assertEqual(open(file_path[:-5]+'test/yyy.txt').read(), 
                          open(file_path[:-5]+'test/test_file_out.txt').read())
        

    def tearDown(self):
        ''' Remove unwanted test files
        '''
        os.remove(file_path[:-5]+'test/yyy.txt')



class TestErrors(unittest.TestCase):
    ''' Test the expected error messages are printed
    '''
    def test_gene_info_error(self):
        ''' esnsure gene information error is printed
        ''' 
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "16:15812194",
                                                 "--header", "--ensembl", "--hg_version",
                                                 "hg38"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')        
        self.assertEqual(ref_seq,'ERROR: No gene information found for query'
                         ' at 16:15812194 in hg38\n')


if __name__ == '__main__':
    unittest.main()
