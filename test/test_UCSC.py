import unittest 
from GeneaPy.UCSC import *


class TestUCSC(unittest.TestCase):
    ''' Test all methods in UCSC.py under the Default
        arguments
    '''
    Default = ScrapeSeq("Test.Default", 20, 20, "hg19", "Y", "Y")
    UpDownstream = ScrapeSeq("Test.UpDownstream", 35, 70, "hg19", "Y")
    HG = ScrapeSeq("Test.HG", 20, 20, "hg38", "Y")
    NoHeader = ScrapeSeq("Test.NoHeader", 20, 20, "hg19")
    var_pos = "15:48762884"

    def test_region(self):
        ''' Ensure the parsed arguments are outputting the expected genomic range
        '''
        self.assertEqual(TestUCSC.Default.create_region(TestUCSC.var_pos),
                         "15:48762864,48762904")
        self.assertEqual(TestUCSC.UpDownstream.create_region(TestUCSC.var_pos),
                         "15:48762849,48762954")
        self.assertEqual(TestUCSC.Default.create_region("15:48762849,48762954"),
                         "15:48762849,48762954")


    def test_region_info(self):
        ''' Check that the expected sequence is being scrapped from UCSC DAS server
        '''
        self.assertEqual(TestUCSC.Default.get_region_info("15:48762864,48762904"),
                                'agcctatctcacactcacagCggaacaggccagggaggttg')
        self.assertEqual(TestUCSC.UpDownstream.get_region_info("15:48762849,48762954"),
                        'ttctgtccagttcgtagcctatctcacactcacagCggaacaggccagggaggttgtgg'
                         'caagttccaaagacacagatgttcggaagggagcactcatcaatatc')
        self.assertEqual(TestUCSC.HG.get_region_info("15:48762864,48762904"),
                        'gaacttaactatatgacaaaAatcacatgaaagatttaagt')


    def test_header(self):
        ''' Ensure the header option works as expected
        '''
        self.assertEqual(TestUCSC.Default.header_option(
            "Test", TestUCSC.var_pos, "15:48762864,48762904",
            'agcctatctcacactcacagCggaacaggccagggaggttg'),
            '> Test 15:48762884 15:48762864,48762904')
        self.assertEqual(TestUCSC.NoHeader.header_option(
            "Test", TestUCSC.var_pos, "15:48762864,48762904",
            'agcctatctcacactcacagCggaacaggccagggaggttg'),
            '')
        
            

class TestErrorHandeling(unittest.TestCase):
    ''' Ensure custom exceptions work as anticipated
    '''
    def test_all_exceptions(self):
        ''' Ensure custom exceptions are raised upon parsing 
            invalid arguments
        '''
        with self.assertRaises(WrongHGversion):
            ScrapeSeq("Test", 20, 20, "hg1", "Y").handle_argument_exception(
            "15:48762864")
        with self.assertRaises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "151671617")
        with self.assertRaises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15::1671617,")
        with self.assertRaises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15:1671617,,1671680")
        with self.assertRaises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15:1671--617")

if __name__ == '__main__':
    unittest.main()
