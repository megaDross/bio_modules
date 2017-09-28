''' Scrape ensembl, via pyensembl, using a genomic position as input'''
import modules.custom_exceptions as ex

def get_canonical_locus_exon(data, contig, position):
    ''' Get the Exon object at a locus of the canonical transcript'''
    canonical_transcript = get_canonical_transcript(data, contig, position)
    exon_at_locus = get_exon_at_locus_of_transcript(data, contig, position,
                                                    canonical_transcript)
    return exon_at_locus

def get_canonical_transcript(data, contig, position):
    ''' Return the canonical transcript at a given position'''
    transcripts = data.transcripts_at_locus(contig=contig, position=position)
    order_by_length = []
    for transcript in transcripts:
        length = transcript.end - transcript.start
        order_by_length.append((length, transcript))

    transcripts_by_length = [x[1] for x in sorted(order_by_length, reverse=True)]
    protein_coding_by_length = [x for x in transcripts_by_length if x.biotype == 'protein_coding']
    if not protein_coding_by_length:
        raise ex.NoProteinCodingTranscript(position, position)
    canonical_transcript = protein_coding_by_length[0]
    return canonical_transcript

def get_exon_at_locus_of_transcript(data, contig, position, transcript):
    ''' Get the Exon object at a locus of a transcript'''
    locus_exons = data.exons_at_locus(contig=contig, position=position)
    transcript_exons = data.exon_ids_of_transcript_id(transcript.id)
    # find overlap between transcript exon list and locus exons
    for locus_exon in locus_exons:
        for transcript_exon in transcript_exons:
            if locus_exon.id == transcript_exon:
                return locus_exon

def get_gene_locus(data, contig, position):
    ''' Get Gene object at a given genomic position'''
    gene_name = data.gene_names_at_locus(contig=contig, position=position)
    if not gene_name:
        raise ex.NoGene(contig, position)
    gene = data.genes_by_name(gene_name[0])[0]
    return gene
