from pyensembl import EnsemblRelease




def get_canonical_transcript(data, contig, position):
    ''' Return the canonical transcript at a given position'''
    transcripts = data.transcripts_at_locus(contig=6, position=29945884)
    order_by_length = []
    for transcript in transcripts:
        length = transcript.end - transcript.start
        order_by_length.append((length, transcript))

    transcripts_by_length = [x[1] for x in sorted(order_by_length, reverse=True)]
    protein_coding_by_length = [x for x in transcripts_by_length if x.biotype == 'protein_coding']
    canonical_transcript = protein_coding_by_length[0]
    return canonical_transcript


def get_exon_at_locus_of_transcript(data, contig, position, transcript):
    ''' Get the Exon object at a locus of a transcript'''
    locus_exons = data.exons_at_locus(contig=contig, position=position)
    transcript_exons = data.exon_ids_of_transcript_id(transcript.id)
    for locus_exon in locus_exons:
        for transcript_exon in transcript_exons:
            if locus_exon.id == transcript_exon:
                return locus_exon


data = EnsemblRelease(83)
canonical_transcript = get_canonical_transcript(data, 6, 29945884)
exon_at_locus = get_exon_at_locus_of_transcript(data, 6, 29945884, canonical_transcript)
print(exon_at_locus)
