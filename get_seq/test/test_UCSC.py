from get_seq.UCSC import *
import pytest 

class TestNormal:
    ''' Test default arguments generally given to the UCSC class
    '''
    
    # init doesnt work for some reason
    Test = ScrapeSeq("Test", 20, 20, "hg19", "Y")
    var_pos = "15:48762884"

    def test_region(self):
        assert TestNormal.Test.create_region(TestNormal.var_pos) == \
                "15:48762864,48762904"


    def test_region_info(self):
        assert TestNormal.Test.get_region_info("15:48762864,48762904") == \
                ('agcctatctcacactcacagCggaacaggccagggaggttg')

    def test_header(self):
        assert TestNormal.Test.header_option(
            "Test", TestNormal.var_pos, "15:48762864,48762904",
            'agcctatctcacactcacagCggaacaggccagggaggttg') == \
            '> Test 15:48762884 15:48762864,48762904'


class TestErrorHandeling:
    ''' Ensure custom exceptions work as anticipated
    '''
    def test_all_exceptions(self):
        with pytest.raises(WrongHGversion):
            ScrapeSeq("Test", 20, 20, "hg1", "Y").handle_argument_exception(
            "15:48762864")
        with pytest.raises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "151671617")
        with pytest.raises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15::1671617,")
        with pytest.raises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15:1671617,,1671680")
        with pytest.raises(TypographyError):
            ScrapeSeq("Test", 20, 20, "hg19", "Y").handle_argument_exception(
                "15:1671--617")
