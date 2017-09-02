import re, ensembl_exon
from pyensembl import EnsemblRelease


class GeneMetaData(object):
    ''' Store gene, transcript and exon metadata of a given genomic position.
    
    Attributes:
        self.query: genomic position
        self.hg_version = human genome version
    '''
    genome = {"hg19": 75, "hg38": 83, 'GRCh37': 75, 'GRCh38': 83}

    def __init__(self, query, hg_version):
        self.query = query.replace("chr","")
        self.hg_version = hg_version 
        self.hg = EnsemblRelease(self.genome.get(hg_version)) # convert to ensembl release object
        self.gene, self.id, self.type, self.location = self.get_gene_info()
        self.transcript = self.get_canonical_transcript()
        self.exon_id, self.intron, self.exon = self.get_exon_info()


    def get_gene_info(self):
        ''' Get the gene information at a given genomic position.

        Returns:
            (Gene Name, Gene Ensembl ID, Gene Type, Gene Location)
        '''
        try: 
            # check if the input is a genomic position or genomic range
            if re.search(r"[-:]", self.query) and self.query.replace(":","").isdigit():
                chrom, pos = [int(x) for x in self.query.split(":")]
                gene_name = self.hg.gene_names_at_locus(contig=chrom, position=pos)
                gene_obj = self.hg.genes_by_name(gene_name[0])[0]
                gene_info = (gene_obj.name, gene_obj.id, gene_obj.biotype, 
                             '{}:{}-{}'.format(gene_obj.contig, gene_obj.start, gene_obj.end))
                return(gene_info)

        except IndexError:
                print("No gene found at {} for genome version {}".format(self.query, str(self.hg_version)))


    def get_canonical_transcript(self):
        ''' Determine and return the canonical transcript of the given gene.
        '''
        all_transcripts = self.hg.transcript_ids_of_gene_name(self.gene)
        all_transcripts = [self.hg.transcript_by_id(x) for x in all_transcripts]
        protein_coding_transcripts = []

        for transcript in all_transcripts:
            size = transcript.end - transcript.start
            if transcript.biotype == "protein_coding":
                protein_coding_transcripts.append((size, transcript.id, transcript.biotype)) 
        # sort by size and return the largest protein coding transcript
        if protein_coding_transcripts:    
            canonical_transcript = sorted(protein_coding_transcripts)[-1][1]
            return canonical_transcript


    def get_exon_info(self):
        ''' Get the exon/intron numbers and IDs.
        '''
        return ensembl_exon.main(self.transcript, self.hg_version, self.query)    


    def update_transcript(self, transcript):
        ''' Update transcript ID and Exon info.

        NOTE: returns No Intron/Exon Matched error. 
        '''
        self.transcript = transcript
        self.exon_id, self.intron, self.exon = self.get_exon_info()










