import subprocess, warnings, os, unittest
import GeneaPy.useful


# use and manipulate file_apth to get path to get_seq.py
file_path = GeneaPy.useful.cwd_file_path(__file__)
get_seq = file_path[:-5]+"GeneaPy/get_seq.py"


#class TestGetSeqPrint(unittest.TestCase):
#    ''' Test the printed output of the get_seq.py click application
#    '''
#
#    def test_seq_tab_file_creation(self):
#        ''' This test the same as test_input_file() except delete the .tab and .seq
#            files associated with the query names in the test_in.txt. This will only pass
#            if the .seq and .tab files are created successfully and checks whether the handle_seq_file method in CompareSeqs works as expected.
#        '''
#        subprocess.Popen("rm /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*CX1*AD*tab /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*24AS*15*tab /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*LX15*NY*tab /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*LX16*SB*tab /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*FUK27*MM*tab /home/david/scripts-x14.04/python/modules/GeneaPy/test/test_files/*LX18*SP*tab", shell=True)
#        TestGetSeqPrint.test_input_file(self)
#
#
#
#    def test_ab1_conversion(self):
#        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "17:48269160",
#                                                "--header", "--seq_file",
#                                                 file_path[:-5]+'test/test_files/'
#                                                'A09_29XX1917_WYN13_F_001.ab1'])
#        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
#        
#        self.assertEqual(ref_seq,
#                         '> query 17:48269160 17:48269140,48269180 -'
#                         '\nReference Sequence:\t'
#                         'tgccctcaccttagcaccatCgttgccgggagcaccgttgg'                    
#                         '\nSanger Sequence:\t'
#                         'tgccctcaccttagcaccatT/Cgttgccgggagcaccgttgg'
#                         '\nthe nucleotides given are DIFFERENT\n\n\n')
#


class GetSeqFileOut(unittest.TestCase):
    ''' Test the written output file of the get_seq.py click application
    '''
    def test_output_file(self):
        '''  Test the written ouput from parsing an input file is as expected.
        '''
        self.maxDiff = None
        # ignore Resource Management error
        warnings.filterwarnings("ignore")
        ref_seq_bytes = subprocess.check_output(["python3", get_seq, file_path[:-5]+
                                                 'test/test_in.txt', "--seq_dir",
                                                 file_path[:-5]+'test/test_files/',
                                                 '--output_file',
                                                  file_path[:-5]+'test/yyy.txt'])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')

        self.assertEqual(open(file_path[:-5]+'test/yyy.txt').read(), 
                          open(file_path[:-5]+'test/test_file_out.txt').read())
        

    def tearDown(self):
        ''' Remove unwanted test files
        '''
        os.remove(file_path[:-5]+'test/yyy.txt')



#class TestErrors(unittest.TestCase):
#    ''' Test the expected error messages are printed
#    '''
#    def test_gene_info_error(self):
#        ''' esnsure gene information error is printed
#        ''' 
#        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "16:15812194",
#                                                 "--header", "--ensembl", "--hg_version",
#                                                 "hg38"])
#        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')        
#        self.assertEqual(ref_seq,'ERROR: No exon information found for query'
#                         ' at 16:15812194 in hg38\nNone\n')


if __name__ == '__main__':
    unittest.main()
