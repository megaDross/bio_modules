import unittest 
from get_seq.transcription_translation import *


######### variables shared across methods #########################
rna_seqs = ["atgcgatcgatcagtttgctgc","atgcgatcgatcagtttgctg",
            "atgcgatcgatcagtttgct", "atgcgtcagttttgggccgcat",
            "cccccggggttttccccggggcccccccc"]

all_options = ProteinRNA("transcribe", "translate", "rc")
no_rc = ProteinRNA("transcribei", "translate")
###################################################################


class TestTranscriptionTranslation(unittest.TestCase):
    ''' Test the transcription() and translation() functions
    '''
    def test_transcription_func(self):
        ''' Ensure the parsed arguments are returning the expceted
            sequence in the presence and then the absence of a start 
            codon
        '''
        self.assertEqual(transcription(rna_seqs[0]), Seq('augcgaucgaucaguuugcugc'))
        self.assertEqual(transcription(rna_seqs[3], "reverse_complement"), 
                                       Seq('augcggcccaaaacugacgcau'))
        self.assertEqual(transcription(rna_seqs[4]), "NOTH")

    def test_translation_func(self):
        ''' Ensure the parsed RNA sequences are returning the expected
            protein sequence and that the recursive recursive element
            of the function works 
        '''
        self.assertEqual(translation(transcription(rna_seqs[0])), Seq('MRSISLL'))
        self.assertEqual(translation(transcription(rna_seqs[1])), Seq('MRSISLL'))
        self.assertEqual(translation(transcription(rna_seqs[2])), Seq('MRSISL'))


class TestProteinRNA(unittest.TestCase):
    ''' Test the ProteinRNA() class
    '''
    def test_get_rna_seq(self):
        ''' Test the options parsed to get_rna_seq() method will return
            the expected sequence
        '''
        self.assertEqual(no_rc.get_rna_seq(rna_seqs[0]),
                         Seq('augcgaucgaucaguuugcugc'))
        self.assertEqual(all_options.get_rna_seq(rna_seqs[3]),
                        Seq('augcggcccaaaacugacgcau'))
        self.assertEqual(no_rc.get_rna_seq(rna_seqs[4]), "NOTH")

    def test_get_protein_seq(self):
        ''' Test the options parsed to get_protein_seq() method will
            return the expected sequence
        '''
        self.assertEqual(no_rc.get_protein_seq(transcription(rna_seqs[0])),
                         Seq('MRSISLL'))
        self.assertEqual(no_rc.get_protein_seq(transcription(rna_seqs[1])),
                         Seq('MRSISLL'))
        self.assertEqual(no_rc.get_protein_seq(transcription(rna_seqs[2])),
                         Seq('MRSISL'))


if __name__ == '__main__':
    unittest.main()
