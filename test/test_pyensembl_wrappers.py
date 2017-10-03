import GeneaPy
from GeneaPy.GeneaPy.modules import pyensembl_wrappers
from pyensembl import EnsemblRelease
import unittest

class TestWrappers(unittest.TestCase):
    def test_get_canonical_transcript(self):
        data = EnsemblRelease(75)
        correct = data.transcript_by_id('ENST00000458208')
        canon = pyensembl_wrappers.get_canonical_transcript(data, 10, 88888)
        self.assertEqual(canon, correct)
