import subprocess
import unittest
from useful_tools import useful


# use and manipulate file_apth to get path to get_seq.py
file_path = useful.cwd_file_path(__file__)
get_seq = file_path[:-5]+"get_seq.py"



class TestGetSeqPrint(unittest.TestCase):
    ''' Test the printed output of the get_seq.py click application
    '''
    def test_simple_get_reference(self):
        ''' Test the printed output from parsing a genomic position is as expected
        '''
        # use subprocess to test get_seq.py click application and then convert to string
        ref_seq_bytes = subprocess.check_output(["python3",get_seq, "chr1:169314424"])
        ref_seq = ref_seq_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(ref_seq,'\nReference Sequence:\taggagcgatgtctcctct'
                                 'ttCattcctgattttgataattt\n\n\n')
    
    def test_simple_get_reference_header(self):
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




if __name__ == '__main__':
    unittest.main()
