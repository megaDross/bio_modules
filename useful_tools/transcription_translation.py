import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re


def transcription(dna, rc=None):
    ''' Transcribes a given DNA sequence into an mRNA sequence from the first intiation codon.
    '''
    # ensure the DNA has no ambigous bases
    dna_seq = Seq(dna, IUPAC.unambiguous_dna)

    # option as whether to perform reverse complement or not before transcribing
    if not rc:
        mRNA= dna_seq.transcribe()
    if rc:
        reverse_comp = dna_seq.reverse_complement()
        mRNA = reverse_comp.transcribe()
    
    # find intiation index and start transcription from this site
    if "aug" in mRNA:
        intiation = str(mRNA).index("aug")
        mRNA = mRNA[intiation:]
    else:
        mRNA = "NOTH"

    return mRNA


def translation(*args):
    ''' A decorator that translates a given RNA sequence into a protein sequence.
    '''
    # the first element in the args tuple is assumed to be an RNA sequence
    rna = args[0]
    
    # if the rna is devisible by 3 then translate the rna
    if (len(rna)/3).is_integer():
        protein = rna.translate()
        return protein
    
    # else, take off the last character of the string and return to inner()
    else:
        rna = rna[:-1]
        return rna








