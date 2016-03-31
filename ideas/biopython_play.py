from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature


my_seq = Seq("GATCCGATGGGCCTATATAGGATCGAAAATCGC", IUPAC.unambiguous_dna)
# reverse complement is required prioir to transcription
# as the template originates from the negative strand 
mRNA = my_seq.reverse_complement().transcribe()
protein = mRNA.translate()

seq_record = SeqRecord(my_seq)
seq_record.id = "Test1"
seq_record.description = "Testing SeqRecord Class"
seq_record.name = "crazy mad test"
# annotations allows one to add a custom made annotation
seq_record.annotations["evidence "] = " None, it's a silly little test you fool"


start_pos = SeqFeature.AfterPosition(5)
end_pos = SeqFeature.BetweenPosition(9, left=8, right=9)
my_location = SeqFeature.FeatureLocation(start_pos, end_pos)
print my_location