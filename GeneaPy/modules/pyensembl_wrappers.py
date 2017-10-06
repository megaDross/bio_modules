''' Scrape ensembl, via pyensembl, using a genomic position as input'''
from GeneaPy.modules.fullexon import FullExon
import GeneaPy.modules.custom_exceptions as ex

def get_transcript(data, contig, position, gene_list=[]):
    ''' Get the canonical or largest Transcript object associated with a given position. '''
    try:
        transcripts = data.transcripts_at_locus(contig=contig, position=position)
        all_transcripts = get_transcripts_by_length(transcripts, gene_list)
        canonical_transcript = get_canonical_transcript(data, contig, position, gene_list)
        if not canonical_transcript in all_transcripts:
            raise ex.NoProteinCodingTranscript(all_transcripts[0])
        canonical_transcript.canonical = True
        return canonical_transcript
    except ex.NoProteinCodingTranscript:
        largest_transcript = all_transcripts[0]
        largest_transcript.canonical = False
        return largest_transcript

def get_transcripts_by_length(transcripts, gene_list=[]):
    ''' Return the all transcripts in order of length from a given position'''
    if gene_list:
        transcripts = [x for x in transcripts if x.gene.name in gene_list]
    transcripts_by_length =  sorted(transcripts, reverse=True, key=lambda x: len(x))
    return transcripts_by_length

def get_canonical_transcript(data, contig, position, gene_list=[]):
    ''' Get the canonical transcript of a gene at a given position'''
    gene = get_gene_locus(data, contig, position, gene_list)
    transcripts_by_length = get_transcripts_by_length(gene.transcripts, gene_list)
    protein_coding_by_length = [x for x in transcripts_by_length if x.biotype == 'protein_coding']
    if not protein_coding_by_length:
        raise ex.NoProteinCodingTranscript(transcripts_by_length[0])
    canonical = protein_coding_by_length[0]
    return canonical

def get_gene_locus(data, contig, position, gene_list=[]):
    ''' Get Gene object at a given genomic position'''
    gene_names = data.gene_names_at_locus(contig=contig, position=position)
    if not gene_names:
        raise ex.NoGene(contig, position)
    try:
        if len(gene_names) > 1:
            raise ex.MultipleGenes(contig, position, gene_names)
    except ex.MultipleGenes as e:
        # should be loggining
        gene_intersect = set(gene_list).intersection(set(gene_names))
        if gene_intersect:
            gene_names = list(gene_intersect)
        else:
            print('ERROR: {}\nINFO: continuing with {}'.format(e, gene_names[0]))
    finally:
        gene = data.genes_by_name(gene_names[0])[0]
        return gene

def get_exon(pos, transcript):
    ''' Return a FullExon object from position and Transcript object
 
    Args:
        pos: position e.g. 48733600 or 15:48733600
        transcript: pyensembl Transcript object
    '''
    total_exons = len(transcript.exons)
    if isinstance(pos, str):
        pos = int(pos.split(':')[1])
    for number, exon in enumerate(transcript.exons, 1):
        if number > len(transcript.exons):
            raise ex.NoExon(transcript.genome.release, transcript.contig, pos)
        next_exon = transcript.exons[number] # need to break if number > t.exons length
        if exon.start <= pos <= exon.end:
            number = '{}/{}'.format(number, total_exons)
            exon = FullExon.from_pyexon(Exon=exon, position=pos, number=number, exon=True)
            return exon
        elif next_exon.start < pos > exon.end and exon.strand == '-':
            number = '{}/{}'.format(number-1, total_exons-1)
            intron = FullExon(exon_id='N/A', contig=exon.contig, start=next_exon.end+1, end=exon.start-1,
                              strand=exon.strand, gene_name=exon.gene_name, gene_id=exon.gene_id,
                              position=pos, number=number, exon=False)
            return intron
        elif exon.end < pos < next_exon.start and exon.strand == '+':
            number = '{}/{}'.format(number, total_exons-1)
            intron = FullExon(exon_id='N/A', contig=exon.contig, start=exon.end+1, end=next_exon.start-1,
                              strand=exon.strand, gene_name=exon.gene_name, gene_id=exon.gene_id,
                              position=pos, number=number, exon=False)
            return intron
