import subprocess, warnings, os, unittest
import GeneaPy.useful


# use and manipulate file_apth to get path to get_seq.py
file_path = GeneaPy.useful.cwd_file_path(__file__)
unknown_primer = file_path[:-5]+"GeneaPy/unknown_primer.py"
header = "\t".join(("Primer", "F_Primer","R_Primer", "Gene", "Product_Size",
                            "Primer_Range","GC%","Number_Amplicons", "\n"))



class TestUnknownPrimerPrint(unittest.TestCase):
    ''' Test the printed output of the unknown_primer.py click application
    '''
    def test_primer_pair(self):
        ''' Test the printed output from parsing a primer pair is as expected
        '''
        warnings.filterwarnings("ignore") 
        primer_bytes = subprocess.check_output(["python3", unknown_primer, 
                                                 "--primers","CCAGGTCCAGATTCAGAGCC",
                                                "GACATGGCGCGGTTACCT"])
        primer = primer_bytes.decode(encoding='UTF-8')
        answer = "\t".join(("query","CCAGGTCCAGATTCAGAGCC","GACATGGCGCGGTTACCT",
                            "RP11-729L2.2","597bp","18:48556413-48557009", "73.23%",
                            "1\n"))
        self.assertEqual(primer, header + answer)
    
    def test_primer_pair_another(self):
        ''' Test the printed output from parsing a primer pair is as expected
        '''
        warnings.filterwarnings("ignore") 
        primer_bytes = subprocess.check_output(["python3", unknown_primer, 
                                                 "--primers","CTGTTCACAGGGCTTGTTCC",
                                                "CTGGGCAGAGAGTCATTTAAAGT"])
        primer = primer_bytes.decode(encoding='UTF-8')
        answer = "\t".join(("query","CTGTTCACAGGGCTTGTTCC","CTGGGCAGAGAGTCATTTAAAGT",
                            "FBN1","421bp","15:48755298-48755718", "38.37%",
                            "1\n"))
        self.assertEqual(primer, header + answer)

    def test_primer_pair_hg38(self):
        ''' Test the printed output is as expected for hg38 option
        '''
        warnings.filterwarnings("ignore") 
        primer_bytes = subprocess.check_output(["python3", unknown_primer, 
                                                 "--primers","CTGTTCACAGGGCTTGTTCC",
                                                "CTGGGCAGAGAGTCATTTAAAGT", 
                                                "--hg_version", "hg38"])
        primer = primer_bytes.decode(encoding='UTF-8')
        answer = "\t".join(("query","CTGTTCACAGGGCTTGTTCC","CTGGGCAGAGAGTCATTTAAAGT",
                            "FBN1","421bp","15:48463101-48463521", "38.37%",
                            "1\n"))
        self.assertEqual(primer, header + answer)


class TestUnknownPrimerIO(unittest.TestCase):
    ''' Test IO processing in unknwon primer
    '''
    def test_in_file(self):
        ''' Test input file prints the expected output
        '''
        warnings.filterwarnings("ignore") 
        primer_bytes = subprocess.check_output(["python3", unknown_primer, 
                                                "--input_file", 
                                                "test_primer_database.txt"]),
        primer = primer_bytes.decode(encoding='UTF-8')
        self.assertEqual(primer, open("test_print_out_unknown_primer.txt").read())


class TestUnknownPrimerExceptions(unittest.TestCase):
    ''' Test the exceptions are raised as expected
    '''
    def test_multiple_amplicons(self):
        ''' Test that the MultipleAmplicons exception is raised as expected
        '''
        warnings.filterwarnings("ignore") 
        primer_bytes = subprocess.check_output(["python3", unknown_primer, 
                                                 "--primers","TACTGACTACGGTGACCAGC",
                                                "TGAAGAACAGCTGGAGGAGG", 
                                                "--hg_version", "hg38"])
        primer = primer_bytes.decode(encoding='UTF-8')
        answer = "The following primers generate more than one amplicon: query\nNone\n"
        self.assertEqual(primer, header + answer)


    

if __name__ == '__main__':
    unittest.main()
        
