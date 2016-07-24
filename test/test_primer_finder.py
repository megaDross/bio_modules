import subprocess, warnings, os, unittest
import GeneaPy.useful


# use and manipulate file_apth to get path to get_seq.py
file_path = GeneaPy.useful.cwd_file_path(__file__)
primer_finder = file_path[:-5]+"GeneaPy/primer_finder.py"


class TestPrimerFinderPrint(unittest.TestCase):
    ''' Test the printed output of the primer_finder.py click application
    '''
    def test_genomic_position(self):
        ''' Test the printed output from parsing a genomic position is as expected
        '''
        warnings.filterwarnings("ignore") 
        # use subprocess to test get_seq.py click application and then convert to string
        primer_bytes = subprocess.check_output(["python3", primer_finder, 
                                                 "chr15:48762884"])
        primer = primer_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(primer, open("primer_finder_test_print.txt").read())
    
    def test_over_700bp(self):
        ''' Test the --size option works as expected
        '''
        primer_bytes = subprocess.check_output(["python3", primer_finder, 
                                                "chr2:189910700", "--size", "700",
                                                "--primer_database", 
                                                "test_primer_database.txt"])
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(primer, 'Variant\tPrimer\tPosition\tGene_Name\tAmplicon_Size\t'
                         'GC%\tAmplicon_Number\tDist_from_F\tDist_from_R\t\n'
                         'query\t-\t-\t-\t-\t-\t-\t-\n')
        
    def test_over_70_GC(self):
        ''' Test the --gc option works as expected
        '''
        primer_bytes = subprocess.check_output(["python3", primer_finder, 
                                                "chr18:48556600", "--gc", "70",
                                                "--primer_database", 
                                                "test_primer_database.txt"])
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(primer, 'Variant\tPrimer\tPosition\tGene_Name\tAmplicon_Size\t'
                         'GC%\tAmplicon_Number\tDist_from_F\tDist_from_R\t\n'
                         'query\t-\t-\t-\t-\t-\t-\t-\n')
 
    def test_over_100_from_primer(self):
        ''' Test the --gc option works as expected
        '''
        primer_bytes = subprocess.check_output(["python3", primer_finder, 
                                                "chr18:48556350", "--distance", "100",
                                                "--primer_database", 
                                                "test_primer_database.txt"])
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(primer, 'Variant\tPrimer\tPosition\tGene_Name\tAmplicon_Size\t'
                         'GC%\tAmplicon_Number\tDist_from_F\tDist_from_R\t\n'
                         'query\t-\t-\t-\t-\t-\t-\t-\n')
 
    def test_piping(self):
        ''' Check that piping a genomic position into get_seq.py works
        '''
        # stop ResourceWarning printing to screen
        warnings.filterwarnings("ignore")

        # use subprocess to test piping a position into the get_seq.py click application,
        # grab the stdout from the resulting tuple and then convert to string
        echo = subprocess.Popen(["echo", "chr15:48762884"], stdout=subprocess.PIPE)
        primer_popen = subprocess.Popen(["python3", primer_finder], stdin=echo.stdout,
                                                               stdout=subprocess.PIPE)
        primer_bytes = primer_popen.communicate()[0]
        primer = primer_bytes.decode(encoding='UTF-8')
        
        self.assertEqual(primer, open("primer_finder_test_print.txt").read())


class TestPrimerFinderFile(unittest.TestCase):
    ''' Test IO works as expected
    '''
    def test_file_in(self):
        ''' Test the expected input files produce the expected output
        '''
        warnings.filterwarnings("ignore")
        primer_bytes = subprocess.check_output(["python3", primer_finder,
                                                "test_in_primer_finder.txt",
                                                "--primer_database", 
                                                "test_primer_database.txt"])
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(primer, open("test_out_primer_finder.txt").read())
    
    def test_file_out(self):
        ''' Test the expected input files produce the expected output
        '''
        warnings.filterwarnings("ignore")
        primer_bytes = subprocess.check_output(["python3", primer_finder,
                                                "test_in_primer_finder.txt",
                                                "--primer_database", 
                                                "test_primer_database.txt", 
                                                "--output_file", "oot.txt"])
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(open("oot.txt").read(), 
                         open("test_out_primer_finder.txt").read())
   
        os.remove(file_path[:-5]+"test/oot.txt")



if __name__ == '__main__':
    unittest.main()
