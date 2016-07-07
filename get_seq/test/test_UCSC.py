import unittest 
from get_seq.UCSC import *


class TestUSCS(unittest.TestCase):
    ''' Test all methods in UCSC.py under the Default
        arguments
    '''
    Default = ScrapeSeq("Test.Default", 20, 20, "hg19", "Y")
    UpDownstream = ScrapeSeq("Test.UpDownstream", 35, 70, "hg19", "Y")
    HG = ScrapeSeq("Test.HG", 20, 20, "hg38", "Y")
    NoHeader = ScrapeSeq("Test.NoHeader", 20, 20, "hg19")
    var_pos = "15:48762884"

    def test_region(self):
        ''' Ensure the parsed arguments are outputting the expected genomic range
        '''
        self.assertEqual(TestUSCS.Default.create_region(TestUSCS.var_pos),
                         "15:48762864,48762904")
        self.assertEqual(TestUSCS.UpDownstream.create_region(TestUSCS.var_pos),
                         "15:48762849,48762954")
        self.assertEqual(TestUSCS.Default.create_region("15:48762849,48762954"),
                         "15:48762849,48762954")


    def test_region_info(self):
        ''' Check that the expected sequence is being scrapped from UCSC DAS server
        '''
        self.assertEqual(TestUSCS.Default.get_region_info("15:48762864,48762904"),
                                'agcctatctcacactcacagCggaacaggccagggaggttg')
        self.assertEqual(TestUSCS.UpDownstream.get_region_info("15:48762849,48762954"),
                        'ttctgtccagttcgtagcctatctcacactcacagCggaacaggccagggaggttgtgg'
                         'caagttccaaagacacagatgttcggaagggagcactcatcaatatc')
        self.assertEqual(TestUSCS.HG.get_region_info("15:48762864,48762904"),
                        'gaacttaactatatgacaaaAatcacatgaaagatttaagt')


    def test_header(self):
        ''' Ensure the header option works as expected
        '''
        self.assertEqual(TestUSCS.Default.header_option(
            "Test", TestUSCS.var_pos, "15:48762864,48762904",
            'agcctatctcacactcacagCggaacaggccagggaggttg'),
            '> Test 15:48762884 15:48762864,48762904')
        self.assertEqual(TestUSCS.NoHeader.header_option(
            "Test", TestUSCS.var_pos, "15:48762864,48762904",
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
