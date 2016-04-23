import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import re

# need a get_ATG decorator?

def transcription(func):
    ''' A decorator that transcribes a given DNA sequence into an mRNA sequence 
    '''
    def inner(dna,rc=None):
        # get the dna sequence by returning the decorated func with the dna variable
        dna = func(dna)
        
        # ignore any header lines
        if dna.startswith(">"):
            dna = dna.split("\n")[1]
        
        # ensure the DNA has no ambigous bases
        dna_seq = Seq(dna, IUPAC.unambiguous_dna)

        # option as whether to perform reverse complement or not before transcribing
        if not rc:
            mRNA= dna_seq.transcribe()
        if rc:
            reverse_comp = dna_seq.reverse_complement()
            mRNA = reverse_comp.transcribe()

        # find intiation index and start transcription from this site
        intiation = str(mRNA).index("aug")
        mRNA = mRNA[intiation:]
        
        return mRNA
    return inner


def translation(func):
    ''' A decorator that translates a given RNA sequence into a protein sequence.
    '''
    def inner(*args):
        # the first element in the args tuple is assumed to be an RNA sequence
        rna = args[0]
        
        # if the rna is devisible by 3 then translate the rna
        if (len(rna)/3).is_integer():
            protein = func(rna).translate()
            return protein
        
        # else, take off the last character of the string and return to inner()
        else:
            rna = rna[:-1]
            return inner(rna)
    
    return inner


#@translation
#@transcription
#def add_polya_tail(dna):
#    ''' Adds a poly A tail to your sequence
#    '''
#    return dna+"AAAAAA"
#
##print(add_polya_tail(">this seq\nCTAGGG"))
#mm = add_polya_tail("CTTATGGGG")
#print(mm)

