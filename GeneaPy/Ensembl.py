import re, ensembl_exon
from pyensembl import EnsemblRelease


class GeneInfo(object):
    ''' Store gene, transcript and exon metadata of a given genomic position.
    
    Attributes:
        self.query: genomic position
        self.hg_version = human genome version
    '''
    genome = {"hg19": 75, "hg38": 83, 'GRCh37': 75, 'GRCh38': 83}

    def __init__(self, query, hg_version):
        self.query = query.replace("chr","")
        self.hg_version = hg_version 
        self.hg = EnsemblRelease(GeneInfo.genome.get(hg_version)) # convert to ensembl release object
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
                print(" ".join(("No gene found at",self.query,"for genome version", str(self.hg_version))))


    def get_canonical_transcript(self):
        ''' Determine and return the canonical transcript of the given gene.
        '''
        all_transcripts = self.hg.transcript_ids_of_gene_name(self.gene)
        all_transcript_details = [self.hg.transcript_by_id(x) for x in all_transcripts]
        protein_coding_transcripts = []
        for x in all_transcript_details:
            split_transcript_info = re.split(r"[=,]",str(x))
            transcript = split_transcript_info[1]
            transcript_type = split_transcript_info[9]
            location = split_transcript_info[-1][:-1]
            start = re.split(r"[:-]", location)[1]
            stop = re.split(r"[:-]", location)[2]
            size = int(stop) - int(start)
            if transcript_type == "protein_coding":
                protein_coding_transcripts.append((size,transcript,transcript_type)) 
        
        # sort by size and return the largest protein coding transcript
        if protein_coding_transcripts:    
            canonical_transcript = sorted(protein_coding_transcripts)[-1][1]
            return canonical_transcript


    def get_exon_info(self):
        ''' From a transcript and genomic position, determine the exon number
            using the ExonInfo class within the Ensmebl module
        '''
        return ensembl_exon.main(self.transcript, self.hg_version, self.query)    














