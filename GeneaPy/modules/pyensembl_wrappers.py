''' Scrape ensembl, via pyensembl, using a genomic position as input'''
from modules.fullexon import FullExon
import modules.custom_exceptions as ex

def get_canonical_transcript(data, contig, position):
    ''' Get the canonical transcript at a given position. '''
    all_transcripts = get_transcripts_by_length(data, contig, position)
    canonical_transcript = extract_canonical_transcript(all_transcripts)
    return canonical_transcript

def get_transcripts_by_length(data, contig, position):
    ''' Return the all transcripts in order of length from a given position'''
    transcripts = data.transcripts_at_locus(contig=contig, position=position)
    order_by_length = []
    for transcript in transcripts:
        length = transcript.end - transcript.start
        order_by_length.append((length, transcript))

    transcripts_by_length = [x[1] for x in sorted(order_by_length, reverse=True, key=lambda x: x[0])]
    return transcripts_by_length

def extract_canonical_transcript(transcripts_by_length):
    ''' Return the longest protein coding transcript from a list of transcripts.'''
    protein_coding_by_length = [x for x in transcripts_by_length if x.biotype == 'protein_coding']
    if not protein_coding_by_length:
        raise ex.NoProteinCodingTranscript(transcripts_by_length[0])
    canonical_transcript = protein_coding_by_length[0]
    return canonical_transcript

def get_gene_locus(data, contig, position):
    ''' Get Gene object at a given genomic position'''
    gene_name = data.gene_names_at_locus(contig=contig, position=position)
    if not gene_name:
        raise ex.NoGene(contig, position)
    gene = data.genes_by_name(gene_name[0])[0]
    return gene

def get_exon(pos, transcript):
    ''' Return a FullExon object from position and Transcript object
 
    Args:
        pos: position e.g. 48733600 or 15:48733600
        transcript: pyensembl Transcript object
    '''
    if isinstance(pos, str):
        pos = int(pos.split(':')[1])
    for number, exon in enumerate(transcript.exons, 1):
        if number > len(transcript.exons):
            raise ex.NoExon(transcript.genome.release, transcript.contig, pos)
        next_exon = transcript.exons[number] # need to break if number > t.exons length
        if exon.start <= pos <= exon.end:
            exon = FullExon.from_pyexon(Exon=exon, position=pos, number=number, exon=True)
            return exon
        # The below operators work with negative strand, they may need to be reversed for positive strand
        elif next_exon.start < pos > exon.end:
            intron = FullExon(exon_id='N/A', contig=exon.contig, start=next_exon.end+1, end=exon.start-1,
                              strand=exon.strand, gene_name=exon.gene_name, gene_id=exon.gene_id,
                              position=pos, number=number-1, exon=False)
            return intron
