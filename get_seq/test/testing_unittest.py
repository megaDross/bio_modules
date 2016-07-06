import unittest 
from get_seq.UCSC import *


class TestNormal(unittest.TestCase):
    Test = ScrapeSeq("Test", 20, 20, "hg19", "Y")
    var_pos = "15:48762884"

    def test_region(self):
        self.assertEqual(TestNormal.Test.create_region(TestNormal.var_pos),
                         "15:48762864,48762904")

if __name__ == '__main__':
    unittest.main()
