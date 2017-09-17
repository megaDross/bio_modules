import re
import argparse
from pyensembl import EnsemblRelease
import ensembl_exon
import UCSC
import custom_exceptions as ex

# TODO: function which iterates through a file of pos that write all metadata to csv
# TODO: unit testing
# TODO: exception handeling
# TODO: logging

class GeneMetaData(object):
    ''' Store gene, transcript and exon metadata of a given genomic position.

    Args:
        query: genomic position
        hg_version: human genome version
        flank: desired number of bp to scrape from each side of the given genomic position (default=50)
        genome: path to human genome FASTA file (optional)

    Attributes:
        self.query: genomic position
        self.hg_version: human genome version
        self.hg: Ensembl human genome release
        self.gene: gene in which the query resides within
        self.id: Ensembl gene ID
        self.type: gene type
        self.location: genomic co-ordinates of the gene
        self.exon_id: Ensembl exon ID
        self.exon: exon number in which the query resides within
        self.intron: intron number in which the query resides within
        self.seq: scrapped FASTA sequence of the query + and - the flank
    '''
    hg_dict = {"hg19": 75, "hg38": 83, 'grch37': 75, 'grch38': 83}

    def __init__(self, query, hg_version, flank=50, genome=None):
        self.query = query.replace("chr", "")
        self.hg_version = hg_version 
        self.hg = EnsemblRelease(self.hg_dict.get(hg_version.lower())) # convert to ensembl release object
        self.gene, self.id, self.type, self.location = self._get_gene_info()
        self.transcript = self._get_canonical_transcript()
        self.exon_id, self.intron, self.exon = self._get_exon_info()
        self.seq = UCSC.main(query, hg_version, genome, flank, flank)


    def _get_gene_info(self):
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


    def _get_canonical_transcript(self):
        ''' Determine and return the canonical transcript of the given gene.
        '''
        try:
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
            else:
                raise ex.NoProteinCodingTranscript(self.query, self.gene)

        except ex.NoProteinCodingTranscript as e:
            return '-'


    def _get_exon_info(self):
        ''' Get the exon/intron numbers and IDs.
        '''
        if self.transcript == '-':
            return ('-', '-', '-')
        else:
            return ensembl_exon.main(self.transcript, self.hg_version, self.query)    


    def update_transcript(self, transcript):
        ''' Update transcript ID and Exon info.
        '''
        self.transcript = transcript
        self.exon_id, self.intron, self.exon = self.get_exon_info()

    
    def print_data(self):
        ''' print all scraped data associated with the query
        '''
        s = '-'*50
        query = '{}\nchr{} in {}\n{}\n\n'.format(
                 s, self.query, self.hg_version,  s)
        gene = 'Gene\nname: {}\nid: {}\nloc: {}\ntype: {}\n\n'.format(
                self.gene, self.id, self.location, self.type)
        transcript = 'Transcript\nid: {}\n\n'.format(
                      self.transcript)
        exon = 'Exon\nid: {}\nexon: {}\nintron: {}\n\n'.format(
                self.exon_id, self.exon, self.intron)
        scrapped_seq = '{}\n\n{}\n\n{}'.format(
                        s, self.seq, s)
        msg = (query+gene+transcript+exon+scrapped_seq)        
        print(msg)



def get_parser():
    parser = argparse.ArgumentParser(description='Scrape a genomic positions meta-data from Ensembl')
    parser.add_argument('query', type=str, help='genomic position')
    parser.add_argument('-hg', '--genome_version', type=str, help='human genome verison (default=hg38)', default='hg38')
    parser.add_argument('-f', '--flank', type=int, help='desired number of bp to scrape from each side of the given genomic position (default=50)', default=50)
    parser.add_argument('-g', '--genome', type=str, help='path to genome FASTA file', default=None) 
    return parser


def cli():
    parser = get_parser()
    args = vars(parser.parse_args())
    metadata = GeneMetaData(args['query'], args['genome_version'], 
                            args['flank'], args['genome'])
    metadata.print_data()


if __name__ == '__main__':
    cli()

    # 'chr15:48733918'
